#!/usr/bin/env python3

import os
import sys
import gzip
import mmap
import struct
import datetime
import json
import socket
import hashlib
from array import array

import shutil
import subprocess
import time
# import pandas as pd

from typing import Callable, Generator, Tuple, Union, List

# import numpy as np

# TODO
# - https://numpy.org/doc/stable/reference/generated/numpy.memmap.html
# - block compress
# - json header
#
# Profiling
# time pypy -m cProfile -s cumulative ./indexer.py
#
# Compress
# for BIN in *.bin; do if [[ ! -f "$BIN.gz" ]]; then gzip -v $BIN -c > $BIN.gz; fi; done
#
#
# K=13; time (gunzip -c -k example-$K.fasta.gz | pypy ./indexer.py )
#
# time pypy ./indexer.py - K=13
#
#
#
# 4.0K example-03.fasta.gz  2.1K example-03.fasta.gz.03.08.bin
#   8K example-05.fasta.gz  3.1K example-05.fasta.gz.05.08.bin
#  80K example-07.fasta.gz   19K example-07.fasta.gz.07.08.bin
# 1.3M example-09.fasta.gz  259K example-09.fasta.gz.09.08.bin
#  19M example-11.fasta.gz  4.1M example-11.fasta.gz.11.08.bin
# 291M example-13.fasta.gz   65M example-13.fasta.gz.13.08.bin
# 4.7G example-15.fasta.gz  1.1G example-15.fasta.gz.15.08.bin
#  74G example-17.fasta.gz   17G example-17.fasta.gz.17.08.bin
#
# 225M S_lycopersicum_chromosomes.4.00.fa.gz
#
#
#
# pypy K=15
# pygzip  C=8 PIPE    C=8 pygz    C=8 STDIN   C=8 Popen
# real    56m11.378s  33m48.443s  35m50.502s  66m 5.744s
# user    41m56.097s  33m36.810s  39m46.295s  44m20.704s
# sys     06m45.105s  00m07.941s  06m21.575s   7m47.551s
#
#
#
# Example
#       real        user        sys        speed
# K=11   0m 6.213s   0m 6.172s   0m0.040s   8,554,287 bp/s
# K=13   2m19.482s   2m16.338s   0m1.341s  11,796,623 bp/s
# K=15  27m21.105s  26m23.109s  0m14.690s  11,800,204 bp/s
# K=17 
# K=19 
#
#
#time pypy ./indexer.py S_lycopersicum_chromosomes.4.00.fa.gz 15 8
#      real         speed
#K= 3   4m41.935s   2,754,109 bp/s
#K= 7   
#K= 9   5m 7.420s   2,576,543 bp/s
#K=11   
#K=13   5m11.103s   2,538,100 bp/s
#K=15   
#K=17   8m24.544s   2,184,808 bp/s
#K=19  
#K=21  


DEBUG      = False
USE_PYTHON = True

ALFA       = 'ACGT'
CONV       = [None] * 255
ALFA_SET   = set(ALFA)
for pos, nuc in enumerate(ALFA):
    CONV[ord(nuc.upper())] = pos
    CONV[ord(nuc.lower())] = pos

# IND         0     1            2            3     4            5     6     7     8
# VALD              1            3                  15                             255
# VALB              1            11                 1111                           11111111
# MAX               0-1          0-3                0-128                          0-255
BIT_MASKS  = [None, 0b0000_0001, 0b0000_0011, None, 0b0000_1111, None, None, None, 0b1111_1111]

# IND          0     1   2  3     4   5     6     7     8
# VALB               1   11       1111                  11111111
# MAX                0-1 0-3      0-128                 0-255
# BIT_MASKS = [None, 1,  3, None, 15, None, None, None, 255]

# HEADER_FMT = '<' + 'B'*2 + 'Q'*4 + 'Q'*(2**8)
# HEADER_LEN = struct.calcsize(HEADER_FMT)


class Timer:
    def __init__(self):
        self.time_begin   = datetime.datetime.now()
        self.time_last    = self.time_begin

        self.val_last     = 0
        self.val_delta    = 0
        self.time_ela     = datetime.timedelta(seconds=0)
        self.time_delta   = datetime.timedelta(seconds=0)
        self.time_ela_s   = 'none'
        self.time_delta_s = 'none'
        self.speed_ela    = 0
        self.speed_delta  = 0

    @property
    def time_delta_seconds(self):
        return (datetime.datetime.now() - self.time_last).total_seconds()

    def update(self, val):
        time_now          = datetime.datetime.now()

        self.time_ela     = time_now       -  self.time_begin
        self.time_delta   = time_now       -  self.time_last

        self.time_ela_s   = str(self.time_ela  ).split('.', 2)[0]
        self.time_delta_s = str(self.time_delta).split('.', 2)[0]

        self.val_delta    = val                -  self.val_last

        self.speed_ela    = int(val            // self.time_ela  .total_seconds())
        self.speed_delta  = int(self.val_delta // self.time_delta.total_seconds())

        self.time_last    = time_now
        self.val_last     = val

    def __str__(self):
        rep = (
            f"ela   time {self.time_ela_s  } val {self.val_last :12,d} speed {self.speed_ela  :12,d}\n"
            f"delta time {self.time_delta_s} val {self.val_delta:12,d} speed {self.speed_delta:12,d}"
        )
        return rep

class Header:
    HEADER     = b'KMER001'
    HEADER_LEN = len(HEADER) + 2*8

    def __init__(self, project_name: str, fasta_file: str, kmer_len: int, field_len: int):
        self.project_name     = project_name
        self.fasta_file       = fasta_file
        self.time_start       = datetime.datetime.now()
        self.kmer_len         = kmer_len
        self.field_len        = field_len

        self.hist             = [0 for _ in range(2**8)]
        self.num_kmers        = 0
        self.num_unique_kmers = 0

        # self.kmer_size        = 4 ** self.kmer_len
        # self.word_len         = (8 // self.field_len)
        # self.data_size        = self.kmer_size * self.field_len // 8
        # self.max_size         = self.HEADER_LEN + self.data_size
        # self.bit_mask         = BIT_MASKS[self.field_len]

    @property
    def hist_sum(self): return sum(self.hist)
    @property
    def hist_uniq(self): return len([h for h in self.hist if h])
    @property
    def kmer_size(self): return 4 ** self.kmer_len
    @property
    def word_len(self): return (8 // self.field_len)
    @property
    def data_size(self): return self.kmer_size * self.field_len // 8
    @property
    def max_size(self): return self.HEADER_LEN + self.data_size
    @property
    def bit_mask(self): return BIT_MASKS[self.field_len]


    def write(self, fhd):
        fhd.seek(0)
        fhd.write(self.HEADER)
        
        ba = bytearray(2*8)
        struct.pack_into('<QQ', ba, 0, self.kmer_len, self.field_len)
        fhd.write(ba)
        time_end   = datetime.datetime.now()

        header_data = {
            "kmer_len"           : self.kmer_len        ,
            "field_len"          : self.field_len       ,
            "num_kmers"          : self.num_kmers       ,
            "num_unique_kmers"   : self.num_unique_kmers,
            "hist_sum"           : self.hist_sum        ,
            "hist_uniq"          : self.hist_uniq       ,
            "hist"               : self.hist            ,
            "project_name"       : self.project_name    , 
            "input_file_name"    : self.fasta_file      ,
            "input_file_size"    : os.path.getsize(self.fasta_file),
            "input_file_ctime"   : os.path.getctime(self.fasta_file),
            "input_file_cheksum" : gen_checksum(self.fasta_file),
            "creation_time_start": str(self.time_start) ,
            "creation_time_end"  : str(time_end)        ,
            "creation_duration"  : str(time_end - self.time_start),
            "hostname"           : socket.gethostname(),
            "checksum_script"    : gen_checksum(os.path.basename(__file__))
        }
        header_json = json.dumps(header_data, indent=None, sort_keys=1, separators=(',', ':'))
        # json_size   = len(header_json)

        fhd.seek(self.max_size + 1)
        # ba = bytearray(8)
        # struct.pack_into('<Q', ba, 0, json_size)
        # f.write(ba)
        fhd.write(header_json.encode())

    def read(self, fhd):
        file_ver = fhd.read(len(self.HEADER))
        assert file_ver == self.HEADER
        print(f"file_ver    : {file_ver.decode()}")

        self.kmer_len, self.field_len = struct.unpack_from("<QQ", fhd.read(2*8))
        print(f"kmer_len    : {self.kmer_len}")
        print(f"field_len   : {self.field_len}")
        print(f"kmer_size   : {self.kmer_size}")
        print(f"data_size   : {self.data_size}")

        fhd.seek(self.HEADER_LEN + self.data_size + 1)
        header_json = fhd.read().decode()
        # print(header_json)

        header_data = json.loads(header_json)
        # print(header_data)

        (
            self.kmer_len           , self.field_len          ,
            self.num_kmers          , self.num_unique_kmers ,
            self.hist,
            self.project_name       ,
            self.input_file_name    , self.input_file_size  , self.input_file_ctime , self.input_file_cheksum,
            self.creation_time_start, self.creation_time_end, self.creation_duration,
            self.hostname           , self.checksum_script
        )= [header_data[k] for k in [
            "kmer_len"           , "field_len"        ,
            "num_kmers"          , "num_unique_kmers" ,
            "hist",
            "project_name"       ,
            "input_file_name"    , "input_file_size"  , "input_file_ctime" , "input_file_cheksum",
            "creation_time_start", "creation_time_end", "creation_duration",
            "hostname"           , "checksum_script"
        ]]

        fhd.seek(0)


def gen_checksum(filename, chunk_size=65536):
    file_hash = hashlib.sha256()
    with open(filename, "rb") as f:
        chunk = f.read(chunk_size)
        while len(chunk) > 0:
            file_hash.update(chunk)
            chunk = f.read(chunk_size)

    return file_hash.hexdigest()

def parse_fasta(fhd, print_every: int=25_000_000) -> Generator[Tuple[str, str, int], None, None]:
    seq_name : str       = None
    seq      : List[str] = []
    seq_num  : int       = 0
    line_num : int       = 0
    bp_num   : int       = 0
    bp_last_e: int       = 0

    timer    :Timer      = Timer()
    conv     :List[Tuple[int,None]] = CONV

    for line in fhd:
        line = line.strip()

        if len(line) == 0:
            continue

        line_num += 1

        if isinstance(line, bytes):
            line = line.decode()

        if line[0] == ">":
            if seq_name is not None:
                if bp_num // print_every > bp_last_e:
                    timer.update(bp_num)
                    bp_last_e   = bp_num // print_every
                    print(f"seq_name: {seq_name:25s} seq_num   : {seq_num        :14,d} line_num  : {line_num          :14,d}")
                    print(f"          {''      :25s} bp_num    : {timer.val_last :14,d} time_ela  : {timer.time_ela_s  :>14s} speed_ela  : {timer.speed_ela  :14,d} bp/s")
                    print(f"          {''      :25s} bp_delta  : {timer.val_delta:14,d} time_delta: {timer.time_delta_s:>14s} speed_delta: {timer.speed_delta:14,d} bp/s")
                seq_str  = (ord(s) for b in seq for s in b)
                seq_str  = [conv[s] for s in seq_str]
                seq_len  = len(seq_str)
                bp_num  += seq_len
                yield seq_name, seq_str, seq_len
            seq_name  = line[1:]
            seq_num  += 1
            seq.clear()
        else:
            seq.append(line)
    
    if seq_name is not None:
        print(f"seq_name: {seq_name:25s} seq_num   : {seq_num        :14,d} line_num  : {line_num          :14,d}")
        print(f"          {''      :25s} bp_num    : {timer.val_last :14,d} time_ela  : {timer.time_ela_s  :>14s} speed_ela  : {timer.speed_ela  :14,d} bp/s")
        print(f"          {''      :25s} bp_delta  : {timer.val_delta:14,d} time_delta: {timer.time_delta_s:>14s} speed_delta: {timer.speed_delta:14,d} bp/s")
        seq_str  = (ord(s) for b in seq for s in b)
        seq_str  = [conv[s] for s in seq_str]
        seq_len  = len(seq_str)
        bp_num  += seq_len
        yield seq_name, seq_str, seq_len

    print(f"seq_name: {'DONE'  :25s} seq_num   : {seq_num        :14,d} line_num  : {line_num          :14,d}")
    print(f"          {''      :25s} bp_num    : {timer.val_last :14,d} time_ela  : {timer.time_ela_s  :>14s} speed_ela  : {timer.speed_ela  :14,d} bp/s")
    print(f"          {''      :25s} bp_delta  : {timer.val_delta:14,d} time_delta: {timer.time_delta_s:>14s} speed_delta: {timer.speed_delta:14,d} bp/s")

def read_fasta(fasta_file: Union[str, None]) -> Generator[Tuple[str, str, int], None, None]:
    filehandler = None

    if fasta_file is None:
        print("READING FASTA FROM STDIN")
        filehandler = lambda: sys.stdin

    else:
        print(f"READING FASTA FROM {fasta_file}")
        filehandler = lambda: open(fasta_file, 'rt')

        if fasta_file.endswith(".gz"):
            if USE_PYTHON:
                print(f"READING FASTA FROM PYGZ {fasta_file}")
                filehandler = lambda: gzip.open(fasta_file, 'rt')

            else:
                print(f"READING FASTA FROM PIPE {fasta_file}")
                filehandler = lambda: os.popen(f"gunzip -k -c {fasta_file}", mode="r")

                # p1     = subprocess.Popen(f"gunzip -k -c {fasta_file}", stdout=subprocess.PIPE, shell=True)
                # for row in parse_fasta(p1.stdout):
                #     yield row
                # print("loop done")

        with filehandler() as fhd:
            for row in parse_fasta(fhd):
                yield row

def gen_kmers(fasta_file: str, kmer_len: int, opener: Callable) -> Generator[Tuple[int, str, int, int, mmap.mmap], None, None]:
    pos_val: List[int] = [4**(kmer_len-p-1) for p in range(kmer_len)]

    for num, (name, seq, seq_len) in enumerate(read_fasta(fasta_file)):
        mm = opener()
        # print(f"{num+1:11,d} {name}")
        # print(f"{num+1:03d} {name} {seq}")

        for i in range(0, seq_len - kmer_len + 1):
            ints:List[Union[int,None]] = seq[i:i+kmer_len]

            if None in ints: continue

            fwd:int = 0
            rev:int = 0
            eints:List[Tuple[int,int]] = enumerate(ints)
            for p, i in eints: fwd += pos_val[         p  ]*   i
            for p, i in eints: rev += pos_val[kmer_len-p-1]*(3-i) 

            # for p, i in enumerate(ints):
            #     fwd += pos_val[         p  ]*   i
            #     rev += pos_val[kmer_len-p-1]*(3-i) 

            # fwd = sum((pos_val[         p  ]*   i  for p, i in enumerate(ints)))
            # rev = sum((pos_val[kmer_len-p-1]*(3-i) for p, i in enumerate(ints)))

            # val = ((pos_val[         p  ]*   i, pos_val[kmer_len-p-1]*(3-i)) for p, i in enumerate(ints))
            # fwd, rev = (sum(f) for f in zip(*val))

            # print(f"{num+1:03d} {name} {i} {kmer} {ints} {ints_p} {fwd:03d} {stni} {stni_p} {rev:03d}")
            # print(f"{num+1:03d} {name} {seq} {fwd:03d} {rev:03d}")
            # print(" ints   ", ints)
            # print(" stni   ", stni)
            # print(" pos_val", pos_val)
            # print(" ints_p ", ints_p, fwd)
            # print(" stni_p ", stni_p, rev)

            yield num, name, fwd, rev, mm

        mm.flush()
        mm.close()

def create_fasta_index(project_name: str, fasta_file: str, index_file: str, kmer_len: int, field_len: int, overwrite: bool, debug: bool = False) -> None:
    header = Header(project_name, fasta_file, kmer_len, field_len)

    print(f"project_name {header.project_name} kmer_len {header.kmer_len:14,d} field_len {header.field_len:14,d} kmer_size {header.kmer_size:14,d} word_len {header.word_len:14,d} max_size {header.max_size:14,d} bytes {header.max_size//1024:14,d} Kb {header.max_size//1024//1024:14,d} Mb {header.max_size//1024//1024//1024:14,d} Gb")
    assert header.bit_mask

    if os.path.exists(index_file) and overwrite:
        os.remove(index_file)

    index_file_tmp = f"{index_file}.tmp"
    if os.path.exists(index_file_tmp):
        os.remove(index_file_tmp)

    with open(index_file_tmp, "w+b") as fhd:
        fhd.seek(header.max_size - 1)
        fhd.write('\0'.encode())
        # print(f"{f.tell():,d} bytes {f.tell()//1024:,d} Kb {f.tell()//1024//1024:,d} Mb {f.tell()//1024//1024//1024:,d} Gb")
        fhd.seek(0)

        # https://docs.python.org/3/library/mmap.html
        # print("initializing database")
        # sys.stdout.flush()
        # mm = [0 for _ in range(header.data_size)]
        # mm = array('B', (0 for _ in range(header.data_size)))
        # print("  done")
        sys.stdout.flush()

        opener = lambda : mmap.mmap(fhd.fileno(), header.max_size, access=mmap.ACCESS_WRITE)

        for num, name, fwd, rev, mm in gen_kmers(fasta_file, kmer_len, opener):
            pos               = fwd if fwd < rev else rev
            header.num_kmers += 1
            # print(num, name, fwd, rev, pos)

            byte_pos        = pos // header.word_len
            bit_pos         = pos %  header.word_len
            # print(f"num {num} name {name} fwd {fwd} rev {rev} pos {pos} byte_pos {byte_pos}")

            bit_shift       = (bit_pos*field_len)                   # number of shifts
            byte_val        = int(mm[header.HEADER_LEN+byte_pos])   # get value
            # byte_val        = 0   # get value
            # byte_val        = mm[byte_pos]   # get value
            bit_val_mask    = int(header.bit_mask<<bit_shift)       # bit mask
            bit_val         = byte_val & bit_val_mask               # bit value
            bit_val_s       = bit_val >> bit_shift                  # bit value shifted

            if bit_val_s == 0:
                header.num_unique_kmers += 1                               # add first appearance of a kmer

            if bit_val_s < (2**field_len-1):
                bit_res         = bit_val_s + 1                         # add to value
                bit_res_s       = bit_res << bit_shift                  # shift back to position
                byte_res        = byte_val & ~ bit_val_mask | bit_res_s # replace old value

                mm[header.HEADER_LEN+byte_pos] = byte_res                      # save
                # mm[byte_pos] = byte_res                      # save

                header.hist[byte_val] -= 1 if byte_val else 0                  # update histogram
                header.hist[byte_res] += 1                                     # update histogram

                if (kmer_len <= 5 and DEBUG) or debug:
                    print(f"  fwd          {fwd         :3d} rev      {rev     :3d} pos       {pos      :3d} word_len     {header.word_len    :3d} ")
                    print(f"  byte_pos     {byte_pos    :3d} bit_pos  {bit_pos :3d} bit_shift {bit_shift:3d}")
                    print(f"  byte_val     {byte_val    :3d} {byte_val    :08b}")
                    print(f"  bit_val_mask {bit_val_mask:3d} {bit_val_mask:08b}")
                    print(f"  bit_val      {bit_val     :3d} {bit_val     :08b}")
                    print(f"  bit_val_s    {bit_val_s   :3d} {bit_val_s   :08b}")
                    print(f"  bit_res      {bit_res     :3d} {bit_res     :08b}")
                    print(f"  bit_res_s    {bit_res_s   :3d} {bit_res_s   :08b}")
                    print(f"  byte_res     {byte_res    :3d} {byte_res    :08b}\n")

        fhd.flush()

        # if (kmer_len <= 5 and DEBUG) or debug:
        #     fhd.seek(0)
        #     mm = mmap.mmap(fhd.fileno(), 0, access=mmap.ACCESS_READ)
        #     for byte_pos, byte_val in enumerate(memoryview(mm[header.HEADER_LEN:])):
        #         print(byte_val, end=" ")
        #     print()

        # print(hist_uniq, hist_sum, hist)
        print(f"project_name {header.project_name} kmer_len {header.kmer_len:14,d} field_len {header.field_len:14,d} num_kmers {header.num_kmers:14,d} num_unique_kmers {header.num_unique_kmers:14,d} kmer_size {header.kmer_size:14,d} word_len {header.word_len:14,d} max_size {header.max_size:14,d} hist_sum {header.hist_sum:14,d} hist_uniq {header.hist_uniq:14,d}")
        if (kmer_len <= 5 and DEBUG) or debug:
            print(f" hist {header.hist}")
        # save header

        print("  done indexing")

        print("  saving data")
        fhd.seek(header.HEADER_LEN+1)
        # for c in mm:
        #     fhd.write(chr(c).encode())
        # mm.tofile(fhd)
        fhd.flush()
        fhd.seek(0)
        print("    done")

        print("  creating header")
        header.write(fhd)
        print("  header saved")
        print("  DONE")

    print("renaming")
    os.rename(index_file_tmp, index_file)

    # print("compressing")
    # with open(index_file_tmp, 'rb') as f_in:
    #     with gzip.open(f'{index_file}.gz', compresslevel=9, mode='wb') as f_out:
    #         shutil.copyfileobj(f_in, f_out)

    # print("cleaning up")
    # os.remove(index_file_tmp)

    print("done")

def read_fasta_index(project_name: str, index_file: str, debug: bool = False) -> None:
    opener = open
    if os.path.exists(f"{index_file}.gz"):
        index_file = index_file + ".gz"
        opener = gzip.open

    header = Header(project_name, None, None, None)
    header_t = Header(project_name, None, None, None)

    with opener(index_file, "r+b") as fhd:
        # read header
        # kmer_len, field_len, f_num_kmers, f_num_unique_kmers, f_hist_sum, f_hist_uniq, *f_hist = struct.unpack_from(HEADER_FMT, f.read(HEADER_LEN))

        header.read(fhd)

        print(f"project_name       : {header.project_name}")
        print(f"input_file name    : {header.input_file_name}")
        print(f"input_file_size    : {header.input_file_size}")
        print(f"input_file_ctime   : {header.input_file_ctime}")
        print(f"input_file_cheksum : {header.input_file_cheksum}")
        print(f"creation_time_start: {header.creation_time_start}")
        print(f"creation_time_end  : {header.creation_time_end}")
        print(f"creation_duration  : {header.creation_duration}")
        print(f"hostname           : {header.hostname}")
        print(f"checksum_script    : {header.checksum_script}")

        print(f"project_name {header.project_name} kmer_len {header.kmer_len:14,d} bytes     {header.max_size//1024:14,d} Kb {header.max_size//1024//1024:14,d} Mb {header.max_size//1024//1024//1024:14,d} Gb")
        print(f"project_name {header.project_name} kmer_len {header.kmer_len:14,d} field_len {header.field_len     :14,d} num_kmers {header.num_kmers:14,d} num_unique_kmers {header.num_unique_kmers:14,d} kmer_size {header.kmer_size:14,d} word_len {header.word_len:14,d} max_size {header.max_size:14,d} hist_sum {header.hist_sum:14,d} hist_uniq {header.hist_uniq:14,d}")

        assert header.bit_mask

        print(f"project_name {project_name} kmer_len {header.kmer_len:14,d} field_len {header.field_len:14,d} hist_sum {header.hist_sum:14,d} hist_uniq {header.hist_uniq:14,d}")
        if (header.kmer_len <= 5 and DEBUG) or debug:
            print(f"  hist {header.hist}")

        fhd.seek(0)

        # https://docs.python.org/3/library/mmap.html
        mm = mmap.mmap(fhd.fileno(), 0, access=mmap.ACCESS_READ)

        if (header.kmer_len <= 5 and DEBUG) or debug:
            fhd.seek(0)
            mm = mmap.mmap(fhd.fileno(), 0, access=mmap.ACCESS_READ)
            for byte_pos, byte_val in enumerate(memoryview(mm[header.HEADER_LEN:header.HEADER_LEN+header.data_size])):
                print(byte_val, end=" ")
            print()

        # https://www.programiz.com/python-programming/methods/built-in/memoryview
        for byte_pos, byte_val in enumerate(memoryview(mm[header.HEADER_LEN:header.HEADER_LEN+header.data_size])):
            # update histogram
            header_t.hist[byte_val] += 1 if byte_val else 0

            for bit_pos in range(header.word_len):
                bit_shift       = (bit_pos*header.field_len)           # number of shifts
                bit_val_mask    = int(header.bit_mask<<bit_shift)      # bit mask
                bit_val         = byte_val & bit_val_mask       # bit value
                bit_val_s       = bit_val >> bit_shift          # bit value shifted
                kmer_pos        = byte_pos * header.word_len + bit_pos # kmer id

                if (header.kmer_len <= 5 and DEBUG) or debug:
                    print(f"  word_len     {header.word_len:3d} ")
                    print(f"  byte_pos     {byte_pos    :3d} bit_pos  {bit_pos :3d} bit_shift {bit_shift:3d}")
                    print(f"  byte_val     {byte_val    :3d} {byte_val    :08b}")
                    print(f"  bit_val_mask {bit_val_mask:3d} {bit_val_mask:08b}")
                    print(f"  bit_val      {bit_val     :3d} {bit_val     :08b}")
                    print(f"  bit_val_s    {bit_val_s   :3d} {bit_val_s   :08b}\n")
        mm.close()

    print(f"project_name {project_name} kmer_len {header.kmer_len:14,d} field_len {header.field_len:14,d} hist_sum {header_t.hist_sum:14,d} hist_uniq {header_t.hist_uniq:14,d}")
    if (header.kmer_len <= 5 and DEBUG) or debug:
        print(f"  hist {header_t.hist}")
    
    assert header.hist_sum  == header_t.hist_sum
    assert header.hist_uniq == header_t.hist_uniq
    assert header.hist      == header_t.hist

def run_test(overwrite=False):
    project_name = "example"
    # kmer_lens =  [3, 5, 7, 9, 11, 13, 15, 17]
    # kmer_lens =  [3]
    # kmer_lens =  [5]
    # kmer_lens =  [7]
    # kmer_lens =  [9]
    # kmer_lens =  [11]
    # kmer_lens =  [13]
    # kmer_lens =  [15]
    kmer_lens =  [17]

    # field_lens = [1, 2, 4, 8]
    # field_lens = [1]
    # field_lens = [2]
    # field_lens = [4]
    field_lens = [8]

    for kmer_len in kmer_lens:
        print(f"project_name {project_name} kmer_len {kmer_len:14,d}")
        for field_len in field_lens:
            fasta_file = f"{project_name}-{kmer_len:02d}.fasta.gz"
            # fasta_file = None # None for stdin

            index_file = f"{project_name}-{kmer_len:02d}.fasta.gz.{kmer_len:02d}.{field_len:02d}.bin"

            print(f"project_name {project_name} kmer_len {kmer_len:14,d} field_len {field_len:14,d}")
            create_fasta_index(project_name, fasta_file, index_file, kmer_len, field_len, overwrite=overwrite)
            read_fasta_index(project_name, index_file)

        print()

def main() -> None:
    if len(sys.argv) == 1:
        run_test(overwrite=False)

    else:
        fasta_file =     sys.argv[1]
        kmer_len   = int(sys.argv[2])
        field_len  = int(sys.argv[3])

        project_name = fasta_file
        index_file   = f"{fasta_file}.{kmer_len:02d}.{field_len:02d}.bin"

        print(f"project_name {project_name:s} fasta_file {fasta_file:s} index_file {index_file:s} kmer_len {kmer_len:14,d} field_len {field_len:14,d}")

        create_fasta_index(project_name, fasta_file, index_file, kmer_len, field_len, overwrite=True)
        
        read_fasta_index(project_name, index_file)

    print()


if __name__ == "__main__":
    main()
