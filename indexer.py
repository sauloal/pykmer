#!/usr/bin/env python3

import os
import sys
import io
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
import tempfile
# import pandas as pd

from typing import Callable, Generator, Tuple, Union, List
import typing

import gc

import numpy as np

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
#K= 3   1,303,590 bp/s
#K= 7
#K= 9
#K=11
#K=13
#K=15
#K=17
#K=19
#K=21  2,263,620 bp/s


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
# BIT_MASKS  = [None, 0b0000_0001, 0b0000_0011, None, 0b0000_1111, None, None, None, 0b1111_1111]

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
    HEADER_VER       = b'KMER001'
    HEADER_VAL_FMT   = '<Q'
    HEADER_VAL_SIZE  = struct.calcsize(HEADER_VAL_FMT)
    HEADER_VAL_NAMES = ['kmer_len']
    HEADER_LEN       = len(HEADER_VER) + HEADER_VAL_SIZE
    HEADER_DATA      = [
        "project_name"       ,
        "kmer_len"           ,
        "input_file_name"    , "input_file_size"  , "input_file_ctime" , "input_file_cheksum",
        "num_kmers"          ,
        "creation_time_start", "creation_time_end", "creation_duration",
        "hostname"           , "checksum_script"  ,
        "hist"               ,
        "hist_sum"           , "hist_count"       , "hist_min"         , "hist_max",
        "vals_sum"           , "vals_count"       , "vals_min"         , "vals_max"

    ]

    def __init__(self, project_name: str, input_file_name: str, kmer_len: int):
        self.project_name           = project_name
        self.input_file_name        = input_file_name
        self.kmer_len               = kmer_len

        self.input_file_size        = None
        self.input_file_ctime       = None
        self.input_file_cheksum     = None

        self.num_kmers              = None

        self._creation_time_start_t = datetime.datetime.now()
        self.creation_time_start    = None
        self.creation_time_end      = None
        self.creation_duration      = None

        self.hostname               = None
        self.checksum_script        = None

        self.hist                   = None
        self.hist_sum               = None
        self.hist_count             = None
        self.hist_min               = None
        self.hist_max               = None

        self.vals_sum               = None
        self.vals_count             = None
        self.vals_min               = None
        self.vals_max               = None

    @property
    def kmer_size(self): return 4 ** self.kmer_len

    @property
    def data_size(self): return self.kmer_size

    @property
    def max_size(self): return self.HEADER_LEN + self.data_size

    @property
    def file_ver(self): return self.HEADER_VER.decode()

    @property
    def header_len(self): return self.HEADER_LEN

    @property
    def max_val(self): return np.iinfo(np.uint8).max

    def _get_mmap(self, fhd, offset=0, mode="r+"):
        return np.memmap(fhd, dtype=np.uint8, mode=mode, offset=offset, shape=(self.data_size,))

    def as_array(self, fhd, mode="r+"):
        npmm = self._get_mmap(fhd, offset=self.header_len, mode=mode)
        return npmm

    def get_empty_array(self, fhd, mode="r+"):
        npmm = self._get_mmap(fhd, offset=0, mode=mode)
        return npmm

    def calc_stats(self, fhd):
        arr             = self.as_array(fhd)
        hist_v, _       = np.histogram(arr, bins=self.max_val, range=(0,self.max_val))

        self.hist       = hist_v.tolist()

        self.hist_sum   = np.sum(hist_v).item()
        self.hist_count = np.count_nonzero(hist_v)
        self.hist_min   = np.min(hist_v).item()
        self.hist_max   = np.max(hist_v).item()

        self.vals_sum   = np.sum(arr).item()
        self.vals_count = np.count_nonzero(arr)
        self.vals_min   = np.min(arr).item()
        self.vals_max   = np.max(arr).item()

    def update(self):
        time_end                 = datetime.datetime.now()

        self.input_file_size     = os.path.getsize(self.input_file_name)
        self.input_file_ctime    = os.path.getctime(self.input_file_name)
        self.input_file_cheksum  = gen_checksum(self.input_file_name)

        self.creation_time_start = str(self._creation_time_start_t)
        self.creation_time_end   = str(time_end)
        self.creation_duration   = str(time_end - self._creation_time_start_t)

        self.hostname            = socket.gethostname()
        self.checksum_script     = gen_checksum(os.path.basename(__file__))

    def write(self, fhd, buffer_size : int  = io.DEFAULT_BUFFER_SIZE):
        _is_filename = False
        if isinstance(fhd, str):
            _is_filename = True
            fhd = open(fhd, "w+b", buffering=buffer_size)

        self.update()
        self.calc_stats(fhd)

        header_data = {k: getattr(self,k) for k in self.HEADER_DATA}

        header_json = json.dumps(header_data, indent=None, sort_keys=1, separators=(',', ':'))

        fhd.seek(0)
        fhd.write(self.HEADER_VER)

        ba = bytearray(self.HEADER_VAL_SIZE)
        struct.pack_into(self.HEADER_VAL_FMT, ba, 0, *(getattr(self, k) for k in self.HEADER_VAL_NAMES))
        fhd.write(ba)

        fhd.seek(self.max_size + 1)
        fhd.write(header_json.encode())

        if _is_filename:
            fhd.flush()
            fhd.close()

    def read(self, fhd):
        file_ver = fhd.read(len(self.HEADER_VER))
        assert file_ver.decode() == self.file_ver

        # print(f"{'file_ver':20s}: {file_ver.decode()}")

        vals = struct.unpack_from(self.HEADER_VAL_FMT, fhd.read(self.HEADER_VAL_SIZE))
        for n, k in enumerate(self.HEADER_VAL_NAMES):
            setattr(self, k, vals[n]) 

        # print(f"{'kmer_len ':20s}: {self.kmer_len :15,d}")
        # print(f"{'kmer_size':20s}: {self.kmer_size:15,d}")
        # print(f"{'data_size':20s}: {self.data_size:15,d}")
        # print(f"{'max_size ':20s}: {self.max_size :15,d}")

        fhd.seek(self.HEADER_LEN + self.data_size + 1)
        header_json = fhd.read().decode()
        # print(header_json)

        header_data = json.loads(header_json)
        # print(header_data)

        for k in self.HEADER_DATA:
            v = header_data[k]
            # print(f"{k:20s}: {str(v)[:50]}")
            setattr(self, k, v)

        fhd.seek(0)

    def check(self, fhd):
        self.read(fhd)

        other = self.__class__(self.project_name, self.input_file_name, self.kmer_len)
        other.read(fhd)
        assert self.project_name     == other.project_name
        assert self.input_file_name  == other.input_file_name
        assert self.kmer_len         == other.kmer_len
        assert self.num_kmers        == other.num_kmers

        other.calc_stats(fhd)
        assert self.hist             == other.hist
        assert self.hist_sum         == other.hist_sum
        assert self.hist_count       == other.hist_count
        assert self.hist_min         == other.hist_min
        assert self.hist_max         == other.hist_max
        assert self.vals_sum         == other.vals_sum
        assert self.vals_count       == other.vals_count
        assert self.vals_min         == other.vals_min
        assert self.vals_max         == other.vals_max

        del other

    def __str__(self):
        res = []
        for k in ["file_ver", "kmer_len", "kmer_size", "data_size", "max_size"] + self.HEADER_DATA:
            v = getattr(self, k)
            if isinstance(v, int):
                res.append(f"{k:20s}: {v:16,d}")
            else:
                res.append(f"{k:20s}: {str(v)[:50]}")
        return "\n".join(res) + "\n"



def gen_checksum(filename, chunk_size=65536):
    file_hash = hashlib.sha256()
    with open(filename, "rb") as f:
        chunk = f.read(chunk_size)
        while len(chunk) > 0:
            file_hash.update(chunk)
            chunk = f.read(chunk_size)

    return file_hash.hexdigest()

def test_np(kmer_len, seq):
    pos_val: List[int] = [4**(kmer_len-p-1) for p in range(kmer_len)]

    seq = tuple(4 if s is None else s for s in seq)

    w = np.array(pos_val, dtype='int64')
    f = np.array(seq, dtype='int8')
    f = np.nan_to_num(f, nan=4, copy=False)
    # print(np.where(f == 4))
    # f = f[438600:438800]
    # print(f == 4)

    print("w.shape", w.shape)
    print("w.dtype", w.dtype)
    print("f.shape", f.shape)
    print("f.dtype", f.dtype)

    for fc, (unique, counts) in enumerate(test_np_ord(w, f)):
        print(" fc", fc+1)
        print("  unique.shape", unique.shape)
        print("  unique.dtype", unique.dtype)
        print("  counts.shape", counts.shape)
        print("  counts.dtype", counts.dtype)
        yield unique, counts
        # for k in ls:
        #     print(" k", k)

def test_np_ord(w, l):
    ws = w.shape[0]
    # print("w", w)

    kmins = []

    for y in range(ws):
        m = l[y:]
        ms = m.shape[0]
        # print("m.shape", ms)

        lt = m[:ms//ws*ws]
        lts = lt.shape[0]
        # print("lt.shape", lts)

        lm = lt.reshape(lts//ws, ws)
        # print("lm.shape", lm.shape)

        lq = np.any(lm == 4, axis=1)

        lr = lm[~lq]
        lrs = lr.shape
        # print("lr", lr)
        # print("lr.shape", lrs)

        rc_lr = lr[:,::-1]
        rc    = 3 - rc_lr
        # print("rc_lr", rc_lr)
        # print("rc   ", rc)

        lv = lr * w
        # print("lv      ", lv)
        # print("lv.shape", lv.shape)

        rc_lv = rc * w
        # print("rc_lv      ", rc_lv)
        # print("rc_lv.shape", rc_lv.shape)

        ls = lv.sum(axis=1)
        # print("ls      ", ls)
        # print("ls.shape", ls.shape)
        # print("ls.dtype", ls.dtype)

        rc_ls = rc_lv.sum(axis=1)
        # print("rc_ls      ", rc_ls)
        # print("rc_ls.shape", rc_ls.shape)
        # print("rc_ls.dtype", rc_ls.dtype)

        # ls.sort()
        # rc_ls.sort()

        kmin = np.minimum(ls, rc_ls)
        # print("kmin      ", kmin)
        # print("kmin.shape", kmin.shape)
        # print("kmin.dtype", kmin.dtype)

        kmins.append(kmin)
        del kmin
        del ls
        del rc_ls
        del lv
        del rc_lv
        del lr
        del rc_lr
        del lm
        del lq
        del lt
        del m

        # unique, counts = np.unique(kmin, return_counts=True)
        # yield unique, counts

    kmins_c = np.concatenate(kmins)
    del kmins
    kmins_c.sort()

    print("summarizing")
    unique, counts = np.unique(kmins_c, return_counts=True)
    print("unique      ", unique)
    print("unique.shape", unique.shape)
    print("unique.dtype", unique.dtype)
    print("counts.shape", counts.shape)
    print("counts.dtype", counts.dtype)
    yield unique, counts
    del unique
    del counts
    del kmins_c

def test_np_example():
    w = np.array([4, 5])
    w
    # array([4, 5])

    l = np.array([3, 3, 2, 2, 0, 4, 5])
    l
    # array([2, 2, 3, 3, 0, 4, 5])

    for y in range(w.shape[0]):
        print(y)

        m = l[y:]
        m
        # 0
        # array([2, 2, 3, 3, 0, 4, 5])
        # 1
        # array([2, 3, 3, 0, 4, 5])

        lt = m[:m.shape[0]//2*2]
        lt
        # array([2, 2, 3, 3, 0, 4])

        lm = lt.reshape(lt.shape[0]//2, 2)
        lm
        # array([
        #     [2, 2],
        #     [3, 3],
        #     [0, 4]])

        lq = np.any(lm == 0, axis=1)
        lq
        # array([False, False,  True])

        lr = lm[~lq]
        lr
        # array([
        #     [2, 2],
        #     [3, 3]])

        lv = lr * w
        lv
        # array([
        #     [12, 15],
        #     [ 8, 10]])

        ls = lv.sum(axis=1)
        ls
        # array([27, 18])

        ls.sort()
        ls
        # array([18, 27])

        # unique, counts = np.unique(ls, return_counts=True)
        # unique
        # # array([18, 27])
        # counts
        # # array([1, 1])

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
                    print(f"seq_name: {seq_name:25s} seq_num   : {seq_num        :15,d} line_num  : {line_num          :15,d}")
                    print(f"          {''      :25s} bp_num    : {timer.val_last :15,d} time_ela  : {timer.time_ela_s  :>14s} speed_ela  : {timer.speed_ela  :15,d} bp/s")
                    print(f"          {''      :25s} bp_delta  : {timer.val_delta:15,d} time_delta: {timer.time_delta_s:>14s} speed_delta: {timer.speed_delta:15,d} bp/s")
                seq_str  = (ord(s) for b in seq for s in b)
                seq_str  = tuple(conv[s] for s in seq_str)
                seq_len  = len(seq_str)
                bp_num  += seq_len
                yield seq_name, seq_str, seq_len
            seq_name  = line[1:]
            seq_num  += 1
            seq.clear()
        else:
            seq.append(line)

    if seq_name is not None:
        print(f"seq_name: {seq_name:25s} seq_num   : {seq_num        :15,d} line_num  : {line_num          :15,d}")
        print(f"          {''      :25s} bp_num    : {timer.val_last :15,d} time_ela  : {timer.time_ela_s  :>14s} speed_ela  : {timer.speed_ela  :15,d} bp/s")
        print(f"          {''      :25s} bp_delta  : {timer.val_delta:15,d} time_delta: {timer.time_delta_s:>14s} speed_delta: {timer.speed_delta:15,d} bp/s")
        seq_str  = (ord(s) for b in seq for s in b)
        seq_str  = tuple(conv[s] for s in seq_str)
        seq_len  = len(seq_str)
        bp_num  += seq_len
        yield seq_name, seq_str, seq_len

    print(f"seq_name: {'DONE'  :25s} seq_num   : {seq_num        :15,d} line_num  : {line_num          :15,d}")
    print(f"          {''      :25s} bp_num    : {timer.val_last :15,d} time_ela  : {timer.time_ela_s  :>14s} speed_ela  : {timer.speed_ela  :15,d} bp/s")
    print(f"          {''      :25s} bp_delta  : {timer.val_delta:15,d} time_delta: {timer.time_delta_s:>14s} speed_delta: {timer.speed_delta:15,d} bp/s")

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

def gen_kmers(fasta_file: str, kmer_len: int) -> Generator[Tuple[int, str, int, int], None, None]:
    pos_val: List[int] = [4**(kmer_len-p-1) for p in range(kmer_len)]

    for num, (name, seq, seq_len) in enumerate(read_fasta(fasta_file)):
        # mm = None
        # print(f"{num+1:11,d} {name}")
        print(f"{num+1:03d} {name} {seq_len:15,d}")

        ints:List[Union[int,None]] = []
        fwd:int = 0
        rev:int = 0
        for i in range(0, seq_len - kmer_len + 1):
            ints:List[Union[int,None]] = seq[i:i+kmer_len]

            if None in ints: continue

            fwd:int = 0
            rev:int = 0
            for p, i in enumerate(ints):
                fwd += pos_val[         p  ]*   i
                rev += pos_val[kmer_len-p-1]*(3-i)

            # print(f"{num+1:03d} {name} {i} {kmer} {ints} {ints_p} {fwd:03d} {stni} {stni_p} {rev:03d}")
            # print(f"{num+1:03d} {name} {fwd:03d} {rev:03d} {ints}")
            # print(" ints   ", ints)
            # print(" stni   ", stni)
            # print(" pos_val", pos_val)
            # print(" ints_p ", ints_p, fwd)
            # print(" stni_p ", stni_p, rev)

            yield num, name, fwd, rev

def process_kmers(
        kmers      : np.ndarray,
        header     : Header,
        index_file : str,
        buffer_size: int             = io.DEFAULT_BUFFER_SIZE,
        frag_size  : int             = 1_000_000_000,
        debug      : bool            = False):

    print(f"    summarizing    {len(kmers):15,d}")
    # print(f"      kmers.shape  {kmers.shape}")
    unique, counts = np.unique(kmers, return_counts=True)

    # print(f"      unique.shape {unique.shape}")
    # print(f"      unique.min   {np.min(unique):15,d}")
    # print(f"      unique.max   {np.max(unique):15,d}")
    # print(f"      counts.shape {counts.shape}")
    # print(f"      counts.min   {np.min(counts):15,d}")
    # print(f"      counts.max   {np.max(counts):15,d}")
    # print(f"      counts.sum   {np.sum(counts):15,d}")
    if debug:
        print( "      unique      ", unique)
        print( "      counts      ", counts)
    # del kmers_a
    gc.collect()
    sys.stdout.flush()

    # https://stackoverflow.com/questions/29611185/avoid-overflow-when-adding-numpy-arrays

    if frag_size > header.data_size: frag_size = header.data_size

    print(f"    creating zeros")
    sys.stdout.flush()

    for frag_from in range(0, header.data_size, frag_size):
        frag_to   = frag_from + frag_size
        if frag_to > header.data_size: frag_to = header.data_size
        if frag_from == frag_to: break

        print(f"        adding [{frag_from:16,d}:{frag_to:16,d}] out of {header.data_size:16,d} [{frag_to-frag_from:16,d}] - {frag_to / header.data_size * 100.0:6.2f} %")

        unique_w = np.where((unique >= frag_from) & (unique < frag_to))[0]
        unique_f = unique[unique_w]
        counts_f = counts[unique_w]

        # print("        unique_w       ", unique_w)
        # print("        unique_w.shape ", unique_w.shape)
        # print("        unique_w.min   ", np.min(unique_w))
        # print("        unique_w.max   ", np.max(unique_w))

        # print("        unique_f       ", unique_f)
        # print("        unique_f.shape ", unique_f.shape)
        # print("        unique_f.min   ", np.min(unique_f))
        # print("        unique_f.max   ", np.max(unique_f))

        # print("        counts_f       ", counts_f)
        # print("        counts_f.shape ", counts_f.shape)
        # print("        counts_f.min   ", np.min(counts_f))
        # print("        counts_f.max   ", np.max(counts_f))
        # print("        counts_f.sum   ", np.sum(counts_f))

        vals_sum  = np.sum(counts_f)

        if  vals_sum == 0:
            print("          zeroes")
            sys.stdout.flush()
        else:
            # print("          saving")

            counts_f[counts_f > header.max_val] = header.max_val

            # print("        counts_f       ", counts_f)
            # print("        counts_f.shape ", counts_f.shape)
            # print("        counts_f.min   ", np.min(counts_f))
            # print("        counts_f.max   ", np.max(counts_f))
            # print("        counts_f.sum   ", np.sum(counts_f))

            vals_frag                       = np.zeros(frag_to-frag_from, dtype=np.uint8)
            vals_frag[unique_f - frag_from] = counts_f.astype(np.uint8)

            # print("        vals_frag      ", vals_frag)
            # print("        vals_frag.shape", vals_frag.shape)
            # print("        vals_frag.min  ", np.min(vals_frag))
            # print("        vals_frag.max  ", np.max(vals_frag))
            # print("        vals_frag.sum  ", np.sum(vals_frag))

            sys.stdout.flush()
            with open(index_file, "w+b", buffering=buffer_size) as fhd:
                npmm      = header.as_array(fhd)

                npmm_frag = npmm[frag_from:frag_to]
                npmm[frag_from:frag_to] = npmm_frag + np.minimum(header.max_val - npmm_frag, vals_frag)

                del npmm_frag
                del npmm
                fhd.flush()
                gc.collect()

            del vals_frag
            gc.collect()

        del unique_w
        del unique_f
        del counts_f
        gc.collect()

def create_fasta_index(
        project_name : str,
        fasta_file   : str,
        index_file   : str,
        kmer_len     : int,
        overwrite    : bool,
        FLUSH_EVERY  : int  =   100_000_000,
        max_frag_size: int  = 1_000_000_000,
        buffer_size  : int  = io.DEFAULT_BUFFER_SIZE,
        debug        : bool = False) -> None:

    header = Header(project_name, fasta_file, kmer_len)

    print(f"project_name {header.project_name} kmer_len {header.kmer_len:15,d} kmer_size {header.kmer_size:15,d} max_size {header.max_size:15,d} bytes {header.max_size//1024:15,d} Kb {header.max_size//1024//1024:15,d} Mb {header.max_size//1024//1024//1024:15,d} Gb")
    # print(header)

    if os.path.exists(index_file) and overwrite:
        os.remove(index_file)

    index_file_tmp = f"{index_file}.tmp"
    if os.path.exists(index_file_tmp):
        os.remove(index_file_tmp)

    print("default buffer size", io.DEFAULT_BUFFER_SIZE)

    # print("opening")
    with open(index_file_tmp, "w+b", buffering=buffer_size) as fhd:
        fhd.seek(header.max_size - 1)
        fhd.write(b'\0')
        # print(f"{f.tell():,d} bytes {f.tell()//1024:,d} Kb {f.tell()//1024//1024:,d} Mb {f.tell()//1024//1024//1024:,d} Gb")

    # https://docs.python.org/3/library/mmap.html

    num_kmers        = 0
    list_pos         = 0
    kmers            = np.zeros(dtype=np.uint64, shape=(FLUSH_EVERY,))
    # kmers            = [None] * FLUSH_EVERY
    last_chrom_num   = None
    frag_size        = header.data_size // 10
    if frag_size > max_frag_size: frag_size = max_frag_size

    # for chrom_num, name, pos, count, mm in gen_kmers(fasta_file, kmer_len, opener):
    for chrom_num, name, fwd, rev in gen_kmers(fasta_file, kmer_len):
        pos               = fwd if fwd < rev else rev
        num_kmers        += 1

        if last_chrom_num != chrom_num or list_pos >= FLUSH_EVERY: #todo: do every N million bp instead of whole chromosomes
            if last_chrom_num != chrom_num: #todo: do every N million bp instead of whole chromosomes
                print(f"  new chrom {name} {list_pos:15,d}")
            else:
                print(f"  {FLUSH_EVERY:15,d} {name} {list_pos:15,d}")

            last_chrom_num = chrom_num
            if list_pos > 0:
                process_kmers(kmers[:list_pos], header, index_file_tmp, buffer_size=buffer_size, frag_size=frag_size, debug=(kmer_len <= 5 and DEBUG) or debug)

                if (kmer_len <= 5 and DEBUG) or debug:
                    print(f"  fwd          {fwd         :3d} rev      {rev     :3d} pos       {pos      :3d}")

                list_pos = 0
                break

        kmers[list_pos] = pos
        list_pos += 1

    if list_pos > 0:
        process_kmers(kmers[:list_pos], header, index_file_tmp, buffer_size=buffer_size, frag_size=frag_size, debug=(kmer_len <= 5 and DEBUG) or debug)
        if (kmer_len <= 5 and DEBUG) or debug:
            print(f"  fwd          {fwd         :3d} rev      {rev     :3d} pos       {pos      :3d} word_len     {header.word_len    :3d} ")

    del kmers
    gc.collect()

    header.num_kmers        = num_kmers

    print(f"project_name {header.project_name} kmer_len {header.kmer_len:15,d}  num_kmers {header.num_kmers:15,d} kmer_size {header.kmer_size:15,d}  max_size {header.max_size:15,d}")

    print("  indexing finished")

    print("  creating header")
    with open(index_file_tmp, "w+b", buffering=buffer_size) as fhd:
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
    header   = Header(project_name, None, None)

    with open(index_file, "r+b") as fhd:
        header.read(fhd)
        print(header)

        print(f"project_name {header.project_name} kmer_len {header.kmer_len:15,d} bytes     {header.max_size//1024:15,d} Kb {header.max_size//1024//1024:15,d} Mb {header.max_size//1024//1024//1024:15,d} Gb")
        print(f"project_name {header.project_name} kmer_len {header.kmer_len:15,d} num_kmers {header.num_kmers:15,d} kmer_size {header.kmer_size:15,d} max_size {header.max_size:15,d}")

        header.check(fhd)
        print("OK")

        # https://docs.python.org/3/library/mmap.html
        npmm = header.as_array(fhd)

        if (header.kmer_len <= 5 and DEBUG) or debug:
            for byte_pos in range(npmm.shape[0]):
                byte_val = npmm[byte_pos].item()
                print(byte_val, end=" ")
            print()

        npmm._mmap.close()
        del npmm

        # https://www.programiz.com/python-programming/methods/built-in/memoryview
        # for byte_pos, byte_val in enumerate(memoryview(mm[header.HEADER_LEN:header.HEADER_LEN+header.data_size])):

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


    for kmer_len in kmer_lens:
        print(f"project_name {project_name} kmer_len {kmer_len:15,d}")
        fasta_file = f"{project_name}-{kmer_len:02d}.fasta.gz"
        # fasta_file = None # None for stdin

        index_file = f"{project_name}-{kmer_len:02d}.fasta.gz.{kmer_len:02d}.bin"

        print(f"project_name {project_name} kmer_len {kmer_len:15,d}")
        create_fasta_index(project_name, fasta_file, index_file, kmer_len, overwrite=overwrite)
        # read_fasta_index(project_name, index_file)

        print()

def main() -> None:
    if len(sys.argv) == 1:
        run_test(overwrite=False)

    else:
        fasta_file =     sys.argv[1]
        kmer_len   = int(sys.argv[2])

        project_name = fasta_file
        index_file   = f"{fasta_file}.{kmer_len:02d}.bin"

        buffer_size = 2**16
        # buffer_size = io.DEFAULT_BUFFER_SIZE

        print(f"project_name {project_name:s} fasta_file {fasta_file:s} index_file {index_file:s} kmer_len {kmer_len:15,d}")

        create_fasta_index(project_name, fasta_file, index_file, kmer_len, buffer_size=buffer_size, overwrite=True, debug=False)

        read_fasta_index(project_name, index_file)

    print()


if __name__ == "__main__":
    main()
