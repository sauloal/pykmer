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

from typing import Callable, Generator, Tuple, Union, List, Dict, NewType, BinaryIO, TextIO
import typing

import gc

import numpy as np


"""
    TODO
    - https://numpy.org/doc/stable/reference/generated/numpy.memmap.html
    - block compress
    - json header

    Profiling
    time pypy -m cProfile -s cumulative ./indexer.py

    Compress
    for BIN in *.bin; do if [[ ! -f "$BIN.gz" ]]; then gzip -v $BIN -c > $BIN.gz; fi; done


    K=13; time (gunzip -c -k example-$K.fasta.gz | pypy ./indexer.py )

    time pypy ./indexer.py - K=13


    4.0K example-03.fasta.gz  2.1K example-03.fasta.gz.03.08.bin
      8K example-05.fasta.gz  3.1K example-05.fasta.gz.05.08.bin
     80K example-07.fasta.gz   19K example-07.fasta.gz.07.08.bin
    1.3M example-09.fasta.gz  259K example-09.fasta.gz.09.08.bin
     19M example-11.fasta.gz  4.1M example-11.fasta.gz.11.08.bin
    291M example-13.fasta.gz   65M example-13.fasta.gz.13.08.bin
    4.7G example-15.fasta.gz  1.1G example-15.fasta.gz.15.08.bin
     74G example-17.fasta.gz   17G example-17.fasta.gz.17.08.bin

    225M S_lycopersicum_chromosomes.4.00.fa.gz



    pypy K=15
    pygzip  C=8 PIPE    C=8 pygz    C=8 STDIN   C=8 Popen
    real    56m11.378s  33m48.443s  35m50.502s  66m 5.744s
    user    41m56.097s  33m36.810s  39m46.295s  44m20.704s
    sys     06m45.105s  00m07.941s  06m21.575s   7m47.551s



    Example
        real        user        sys        speed
    K=11   0m 6.213s   0m 6.172s   0m0.040s   8,554,287 bp/s
    K=13   2m19.482s   2m16.338s   0m1.341s  11,796,623 bp/s
    K=15  27m21.105s  26m23.109s  0m14.690s  11,800,204 bp/s
    K=17
    K=19


    K=15; F=S_lycopersicum_chromosomes.4.00.fa.gz.$K.kin; bgzip -i -I $F.gz.gzi -l 9 -c $F > $F.gz

    time pypy ./indexer.py S_lycopersicum_chromosomes.4.00.fa.gz 15
            real         speed         size  bgzip
    K= 3    16m21.205s   797,621 bp/s    2K   69b
    K= 5    16m 6.006s   809,751 bp/s    3K   94b
    K= 7    16m27.471s   787,715 bp/s   16K  287b
    K= 9
    K=11
    K=13
    K=15  29m 3.686s   481,910 bp/s      1G  156M
    K=17  88m12.597s   167,771 bp/s     17G  574M
    X K=19                             257G
    X K=21                               4T
"""


"""
    IND         0     1            2            3     4            5     6     7     8
    VALD              1            3                  15                             255
    VALB              1            11                 1111                           11111111
    MAX               0-1          0-3                0-128                          0-255
    BIT_MASKS  = [None, 0b0000_0001, 0b0000_0011, None, 0b0000_1111, None, None, None, 0b1111_1111]

    IND          0     1   2  3     4   5     6     7     8
    VALB               1   11       1111                  11111111
    MAX                0-1 0-3      0-128                 0-255
    BIT_MASKS = [None, 1,  3, None, 15, None, None, None, 255]

    HEADER_FMT = '<' + 'B'*2 + 'Q'*4 + 'Q'*(2**8)
    HEADER_LEN = struct.calcsize(HEADER_FMT)
"""


DEBUG      = False
USE_PYTHON = True


ALFA       = 'ACGT'
CONV       = [None] * 255
ALFA_SET   = set(ALFA)
for pos, nuc in enumerate(ALFA):
    CONV[ord(nuc.upper())] = pos
    CONV[ord(nuc.lower())] = pos


Datetime    = NewType('Datetime'   , datetime.datetime)
Hist        = NewType('Hist'       , List[int])
Chromosomes = NewType('Chromosomes', List[Tuple[str, int]])


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
    HEADER_VER      : bytes     = 1
    # HEADER_VER      : bytes     = b'KMER001'
    # HEADER_VAL_FMT  : str       = '<Q'
    # HEADER_VAL_SIZE : int       = struct.calcsize(HEADER_VAL_FMT)
    HEADER_VAL_NAMES: List[str] = ['kmer_len']
    # HEADER_LEN      : int       = len(HEADER_VER) + HEADER_VAL_SIZE
    HEADER_DATA     : List[str] = [
        "project_name"       ,
        "kmer_len"           ,
        "input_file_name"    , "input_file_size"  , "input_file_ctime" , "input_file_cheksum" ,
                               "output_file_size" , "output_file_ctime", "output_file_cheksum",
        "num_kmers"          , "chromosomes"      ,
        "creation_time_start", "creation_time_end", "creation_duration",
        "hostname"           , "checksum_script"  ,
        "hist"               ,
        "hist_sum"           , "hist_count"       , "hist_min"         , "hist_max",
        "vals_sum"           , "vals_count"       , "vals_min"         , "vals_max"
    ]
    EXT = 'kin'
    TMP = 'tmp'

    def __init__(self,
            project_name   : str,
            input_file     : Union[str, None] = None,
            kmer_len       : Union[int, None] = None,
            index_file     : Union[str, None] = None,
            buffer_size    : int = io.DEFAULT_BUFFER_SIZE):

        self.project_name   :str         = project_name
        self.input_file_name:str         = input_file
        self.kmer_len       :int         = kmer_len
        self._buffer_size   :int         = buffer_size

        if index_file is not None:
            self._parse_index_file_name(index_file)

        assert self.kmer_len
        assert self.kmer_len > 0
        assert self.kmer_len % 2 == 1

        self.input_file_size       :int         = None
        self.input_file_ctime      :float       = None
        self.input_file_cheksum    :str         = None

        self.ouput_file_size       :int         = None
        self.ouput_file_ctime      :float       = None
        self.ouput_file_cheksum    :str         = None

        self.num_kmers             :int         = None
        self.chromosomes           :Chromosomes = None

        self._creation_time_start_t:Datetime    = datetime.datetime.now()
        self.creation_time_start   :str         = None
        self.creation_time_end     :str         = None
        self.creation_duration     :str         = None

        self.hostname              :str         = None
        self.checksum_script       :str         = None

        self.hist                  :Hist        = None
        self.hist_sum              :int         = None
        self.hist_count            :int         = None
        self.hist_min              :int         = None
        self.hist_max              :int         = None

        self.vals_sum              :int         = None
        self.vals_count            :int         = None
        self.vals_min              :int         = None
        self.vals_max              :int         = None

    @property
    def index_file(self) -> str: return f"{self.input_file_name}.{self.kmer_len:02d}.{self.EXT}"

    @property
    def index_tmp_file(self) -> str: return f"{self.index_file}.{self.TMP}"

    @property
    def metadata_file(self) -> str: return f"{self.index_file}.json"

    @property
    def kmer_size(self) -> int: return 4 ** self.kmer_len

    @property
    def data_size(self) -> int: return self.kmer_size

    @property
    def max_size(self) -> int: return self.data_size

    @property
    def file_ver(self) -> str: return self.HEADER_VER

    @property
    def max_val(self) -> int: return np.iinfo(np.uint8).max

    # @property
    # def header_len(self) -> int: return self.HEADER_LEN

    def _parse_index_file_name(self, index_file):
        ext = index_file[(2 + 1) + (len(self.EXT) + 1) * -1:]

        if self.input_file_name is None:
            self.input_file_name = ext[len(self.EXT) * -1:]
        
        if self.kmer_len is None:
            self.kmer_len        = int(ext[1:3])

    def _get_mmap(self, fhd: BinaryIO, offset: int = 0, mode: str ="r+") -> np.memmap:
        return np.memmap(fhd, dtype=np.uint8, mode=mode, offset=offset, shape=(self.data_size,))


    def update_stats(self, fhd: BinaryIO) -> None:
        print("updating stats")

        arr             = self.get_array_from_fhd(fhd)
        hist_v, _       = np.histogram(arr, bins=self.max_val, range=(1,self.max_val))

        self.hist       = hist_v.tolist()
        # print(self.hist)

        self.hist_sum   = np.sum(hist_v).item()
        self.hist_count = np.count_nonzero(hist_v)
        self.hist_min   = np.min(hist_v).item()
        self.hist_max   = np.max(hist_v).item()

        self.vals_sum   = np.sum(arr).item()
        self.vals_count = np.count_nonzero(arr)
        self.vals_min   = np.min(arr).item()
        self.vals_max   = np.max(arr).item()

    def update_stats_index_file(self) -> None:
        with self.open_index_file() as fhd:
            self.update_stats(fhd)

    def update_stats_index_tmp_file(self) -> None:
        with self.open_index_tmp_file() as fhd:
            self.update_stats(fhd)


    def update_metadata(self, index_file) -> None:
        print("updating metadata")

        self.input_file_size     = os.path.getsize(self.input_file_name)
        self.input_file_ctime    = os.path.getctime(self.input_file_name)
        self.input_file_cheksum  = gen_checksum(self.input_file_name)

        self.output_file_size    = os.path.getsize(index_file)
        self.output_file_ctime   = os.path.getctime(index_file)
        self.output_file_cheksum = gen_checksum(index_file)

        self.hostname            = socket.gethostname()
        self.checksum_script     = gen_checksum(os.path.basename(__file__))

        time_end                 = datetime.datetime.now()
        self.creation_time_start = str(self._creation_time_start_t)
        self.creation_time_end   = str(time_end)
        self.creation_duration   = str(time_end - self._creation_time_start_t)


    def open_file(self, index_file: str, mode: str = "r+b") -> BinaryIO:
        return open(index_file, mode, buffering=self._buffer_size)

    def open_index_file(self, mode: str = "r+b") -> BinaryIO:
        return self.open_file(self.index_file, mode=mode)

    def open_index_tmp_file(self, mode: str = "r+b") -> BinaryIO:
        return self.open_file(self.index_tmp_file, mode=mode)


    def _init_clean(self, overwrite: bool=False) -> None:
        if os.path.exists(self.index_file):
            if overwrite:
                os.remove(self.index_file)
            else:
                raise ValueError(f"file {self.index_file} already exists and overwritting disabled")

        if os.path.exists(self.metadata_file):
            os.remove(self.metadata_file)

        if os.path.exists(self.index_tmp_file):
            os.remove(self.index_tmp_file)

    def init_file(self, index_file: str, mode: str = "r+b") -> None:
        # print("opening")
        if not os.path.exists(index_file):
            with open(index_file, 'w'):
                pass

        with self.open_file(index_file, mode=mode) as fhd:
            fhd.seek(self.max_size - 1)
            fhd.write(b'\0')
        # print(f"{f.tell():,d} bytes {f.tell()//1024:,d} Kb {f.tell()//1024//1024:,d} Mb {f.tell()//1024//1024//1024:,d} Gb")

    def init_index_file(self, overwrite: bool=False, mode: str = "r+b"):
        self._init_clean(self, overwrite=overwrite)
        return self.init_file(self.index_file, mode=mode)

    def init_index_tmp_file(self, overwrite: bool=False, mode: str = "r+b"):
        self._init_clean(overwrite=overwrite)
        return self.init_file(self.index_tmp_file, mode=mode)


    def get_array_from_fhd(self, fhd: BinaryIO, mode: str = "r+") -> np.memmap:
        npmm = self._get_mmap(fhd, offset=0, mode=mode)
        return npmm

    def get_array_from_index_file(self, fhd_mode: str = "r+b", mm_mode: str = "r+") -> np.memmap:
        with self.open_index_file(mode=fhd_mode) as fhd:
            yield self.get_array_from_fhd(fhd, mode=mm_mode)

    def get_array_from_index_tmp_file(self, fhd_mode: str = "r+b", mm_mode: str = "r+") -> np.memmap:
        with self.open_index_tmp_file(mode=fhd_mode) as fhd:
            yield self.get_array_from_fhd(fhd, mode=mm_mode)


    def write_metadata_file(self, index_file: str) -> None:
        assert self.num_kmers
        assert self.chromosomes

        self.update_metadata(index_file)
        with self.open_file(index_file) as fhd:
            self.update_stats(fhd)

        header_data = {k: getattr(self,k) for k in self.HEADER_DATA}

        with open(self.metadata_file, 'wt') as fhm:
            json.dump(header_data, fhm, indent=1, sort_keys=1)

    def write_metadata_index_file(self) -> None:
        self.write_metadata_file(self.index_file)

    def write_metadata_index_tmp_file(self) -> None:
        self.write_metadata_file(self.index_tmp_file)


    def read_metadata(self):
        with open(self.metadata_file, 'rt') as fhd:
            header_data = json.load(fhd)

        # print(header_data)

        for k in self.HEADER_DATA:
            v = header_data[k]
            # print(f"{k:20s}: {str(v)[:50]}")
            setattr(self, k, v)


    def check_data(self, fhd: BinaryIO) -> None:
        self.read_metadata()

        other = self.__class__(self.project_name, self.input_file_name, self.kmer_len)
        other.read_metadata()
        assert self.project_name     == other.project_name
        assert self.input_file_name  == other.input_file_name
        assert self.kmer_len         == other.kmer_len
        assert self.num_kmers        == other.num_kmers

        other.update_stats(fhd)
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

    def check_data_file(self, filename: str) -> None:
        with self.open_file(filename, mode="r+b") as fhd:
            self.check_data(fhd)

    def check_data_index(self) -> None:
        self.check_data_file(self.index_file)

    def check_data_index_tmp(self) -> None:
        self.check_data_file(self.index_tmp_file)


    def __str__(self) -> str:
        res = []
        for k in ["file_ver", "kmer_len", "kmer_size", "data_size", "max_size"] + self.HEADER_DATA:
            v = getattr(self, k)
            if isinstance(v, int):
                res.append(f"{k:20s}: {v:15,d}")
            else:
                res.append(f"{k:20s}: {str(v)[:50]}")
        return "\n".join(res) + "\n"


def gen_checksum(filename: str, chunk_size: int = 2**16) -> str:
    file_hash = hashlib.sha256()
    with open(filename, "rb") as f:
        chunk = f.read(chunk_size)
        while len(chunk) > 0:
            file_hash.update(chunk)
            chunk = f.read(chunk_size)

    return file_hash.hexdigest()

def test_np(kmer_len: int, seq: str) -> None:
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

def test_np_ord(w, l) -> None:
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

def test_np_example() -> None:
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

def parse_fasta(fhd: TextIO, print_every: int = 25_000_000) -> Generator[Tuple[str, str, int], None, None]:
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

def read_fasta(input_file: Union[str, None]) -> Generator[Tuple[str, str, int], None, None]:
    filehandler = None

    if input_file is None:
        print("READING FASTA FROM STDIN")
        filehandler = lambda: sys.stdin

    else:
        print(f"READING FASTA FROM {input_file}")
        filehandler = lambda: open(input_file, 'rt')

        if input_file.endswith(".gz"):
            if USE_PYTHON:
                print(f"READING FASTA FROM PYGZ {input_file}")
                filehandler = lambda: gzip.open(input_file, 'rt')

            else:
                print(f"READING FASTA FROM PIPE {input_file}")
                filehandler = lambda: os.popen(f"gunzip -k -c {input_file}", mode="r")

                # p1     = subprocess.Popen(f"gunzip -k -c {input_file}", stdout=subprocess.PIPE, shell=True)
                # for row in parse_fasta(p1.stdout):
                #     yield row
                # print("loop done")

        with filehandler() as fhd:
            for row in parse_fasta(fhd):
                yield row

def gen_kmers(input_file: str, kmer_len: int) -> Generator[Tuple[int, str, int, int, int], None, None]:
    pos_val: List[int] = [4**(kmer_len-p-1) for p in range(kmer_len)]

    for chrom_num, (name, seq, seq_len) in enumerate(read_fasta(input_file)):
        # mm = None
        # print(f"{chrom_num+1:11,d} {name}")
        print(f"{chrom_num+1:03d} {name} {seq_len:15,d}")

        ints:List[Union[int,None]] = []
        fwd:int = 0
        rev:int = 0
        for chrom_kmer_num in range(0, seq_len - kmer_len + 1):
            ints:List[Union[int,None]] = seq[chrom_kmer_num:chrom_kmer_num+kmer_len]

            if None in ints: continue

            fwd:int = 0
            rev:int = 0
            for p, j in enumerate(ints):
                fwd += pos_val[         p  ]*   j
                rev += pos_val[kmer_len-p-1]*(3-j)

            # print(f"{num+1:03d} {name} {i} {kmer} {ints} {ints_p} {fwd:03d} {stni} {stni_p} {rev:03d}")
            # print(f"{num+1:03d} {name} {fwd:03d} {rev:03d} {ints}")
            # print(" ints   ", ints)
            # print(" stni   ", stni)
            # print(" pos_val", pos_val)
            # print(" ints_p ", ints_p, fwd)
            # print(" stni_p ", stni_p, rev)

            yield chrom_num, name, seq_len, chrom_kmer_num, fwd, rev

def process_kmers(
        kmers      : np.ndarray,
        header     : Header,
        frag_size  : int  = 1_000_000_000,
        debug      : bool = False) -> np.ndarray:

    print(f"    summarizing       {len(kmers):15,d}")
    unique, counts = np.unique(kmers, return_counts=True)

    if debug:
        print(f"      kmers.shape     {kmers.shape[0]:15,d}")
        print(f"      kmers.min       {np.max(kmers):15,d}")
        print(f"      kmers.sum       {np.sum(kmers):15,d}")
        print(f"      unique.shape    {unique.shape[0]:15,d}")
        print(f"      unique.min      {np.min(unique):15,d}")
        print(f"      unique.max      {np.max(unique):15,d}")
        print(f"      counts.shape    {counts.shape[0]:15,d}")
        print(f"      counts.min      {np.min(counts):15,d}")
        print(f"      counts.max      {np.max(counts):15,d}")
        print(f"      counts.sum      {np.sum(counts):15,d}")
        print( "      unique         ", unique)
        print( "      counts         ", counts)

    # del kmers_a
    gc.collect()
    sys.stdout.flush()

    # https://stackoverflow.com/questions/29611185/avoid-overflow-when-adding-numpy-arrays

    if frag_size > header.data_size: frag_size = header.data_size

    # print(f"    creating zeros")
    # sys.stdout.flush()

    hist = None
    for frag_from in range(0, header.data_size, frag_size):
        frag_to   = frag_from + frag_size
        if frag_to > header.data_size: frag_to = header.data_size

        print(f"        adding [{frag_from:15,d}:{frag_to:15,d}] out of {header.data_size:15,d} (length {frag_to-frag_from:15,d} frag_size {frag_size:15,d}) - {frag_to / header.data_size * 100.0:6.2f} %")

        if frag_from == frag_to:
            print(f"          frag_from ({frag_from}) == ({frag_to}) frag_to - breaking")
            sys.stdout.flush()
            break

        unique_w = np.where((unique >= frag_from) & (unique < frag_to))[0]
        unique_f = unique[unique_w]
        counts_f = counts[unique_w]

        if debug:
            print("        unique_w          ", unique_w)
            print("        unique_w.shape    ", unique_w.shape[0])
            print("        unique_w.min      ", np.min(unique_w))
            print("        unique_w.max      ", np.max(unique_w))

            print("        unique_f          ", unique_f)
            print("        unique_f.shape    ", unique_f.shape[0])
            print("        unique_f.min      ", np.min(unique_f))
            print("        unique_f.max      ", np.max(unique_f))

            print("        counts_f          ", counts_f)
            print("        counts_f.shape    ", counts_f.shape[0])
            print("        counts_f.min      ", np.min(counts_f))
            print("        counts_f.max      ", np.max(counts_f))
            print("        counts_f.sum      ", np.sum(counts_f))

        vals_sum  = np.sum(counts_f)

        if  vals_sum == 0:
            print("          zeroes")
            sys.stdout.flush()

        else:
            # print("          saving")
            # sys.stdout.flush()

            counts_f[counts_f > header.max_val] = header.max_val

            if debug:
                print("        counts_f          ", counts_f)
                print("        counts_f.shape    ", counts_f.shape[0])
                print("        counts_f.min      ", np.min(counts_f))
                print("        counts_f.max      ", np.max(counts_f))
                print("        counts_f.sum      ", np.sum(counts_f))

            vals_frag                       = np.zeros(frag_to-frag_from, dtype=np.uint8)
            vals_frag[unique_f - frag_from] = counts_f.astype(np.uint8)

            if debug:
                print("        vals_frag         ", vals_frag)
                print("        vals_frag.shape   ", vals_frag.shape[0])
                print("        vals_frag.min     ", np.min(vals_frag))
                print("        vals_frag.max     ", np.max(vals_frag))
                print("        vals_frag.sum     ", np.sum(vals_frag))

            with header.open_index_tmp_file() as fhd:
                npmm                    = header.get_array_from_fhd(fhd)
                npmm_frag               = npmm[frag_from:frag_to]
                hist_before             = np.histogram(npmm_frag, bins=np.iinfo(np.uint8).max, range=(1,np.iinfo(np.uint8).max))[0]
                npmm[frag_from:frag_to] = npmm_frag + np.minimum(header.max_val - npmm_frag, vals_frag)
                hist_after              = np.histogram(npmm_frag, bins=np.iinfo(np.uint8).max, range=(1,np.iinfo(np.uint8).max))[0]

                # if hist is None: hist_before[0] = 2 * np.iinfo(np.uint8).max + 2
                hist_diff               = hist_after - hist_before

                if debug:
                    print("        hist_before       ", hist_before)
                    print("        hist_after        ", hist_after)
                    print("        hist_diff         ", hist_diff)
                # print(npmm.tolist())
                
                if hist is None: hist   = hist_diff
                else:            hist  += hist_diff

                del hist_before
                del hist_after
                del hist_diff
                del npmm_frag
                del npmm
                fhd.flush()
                gc.collect()

            # with open(index_file, "r+b", buffering=buffer_size) as fhd:
            #     npmm                    = header.as_array(fhd)
            #     print(npmm.tolist())

            del vals_frag
            gc.collect()

        del unique_w
        del unique_f
        del counts_f
        gc.collect()

    return hist

def create_fasta_index(
        project_name : str,
        input_file   : str,
        kmer_len     : int,
        overwrite    : bool,
        FLUSH_EVERY  : int  =   100_000_000,
        min_frag_size: int  =   500_000_000,
        max_frag_size: int  = 1_000_000_000,
        buffer_size  : int  = io.DEFAULT_BUFFER_SIZE,
        debug        : bool = False) -> None:

    header = Header(project_name, input_file=input_file, kmer_len=kmer_len, buffer_size=buffer_size)

    print(f"project_name {header.project_name} kmer_len {header.kmer_len:15,d} kmer_size {header.kmer_size:15,d} max_size {header.max_size:15,d} bytes {header.max_size//1024:15,d} Kb {header.max_size//1024//1024:15,d} Mb {header.max_size//1024//1024//1024:15,d} Gb")
    # print(header)

    print("default buffer size", buffer_size)

    header.init_index_tmp_file(overwrite=overwrite)

    # https://docs.python.org/3/library/mmap.html

    num_kmers        = 0
    list_pos         = 0
    kmers            = np.zeros(dtype=np.uint64, shape=(FLUSH_EVERY,))
    hist             = None
    last_chrom_num   = None
    frag_size        = header.data_size // 10

    if frag_size > max_frag_size   : frag_size = max_frag_size
    if frag_size < min_frag_size   : frag_size = min_frag_size
    if frag_size > header.data_size: frag_size = header.data_size

    # for chrom_num, name, pos, count, mm in gen_kmers(input_file, kmer_len, opener):
    chromosomes = []
    for chrom_num, name, seq_len, chrom_kmer_num, fwd, rev in gen_kmers(input_file, kmer_len):
        pos               = fwd if fwd < rev else rev
        num_kmers        += 1
        # if chrom_kmer_num >= 1_000_000: continue

        if last_chrom_num != chrom_num or list_pos >= FLUSH_EVERY: #todo: do every N million bp instead of whole chromosomes
            last_list_pos = list_pos

            if last_chrom_num == chrom_num: #todo: do every N million bp instead of whole chromosomes
                print(f"  {FLUSH_EVERY:15,d} {name} {last_list_pos:15,d}")

            if list_pos > 0:
                hist_k = process_kmers(kmers[:list_pos], header, frag_size=frag_size, debug=(kmer_len <= 5 and DEBUG) or debug)

                if hist is None: hist  = hist_k
                else:            hist += hist_k

                # print("hist_k     ", hist_k)
                # print("hist       ", hist)

                if (kmer_len <= 5 and DEBUG) or debug:
                    print(f"  fwd          {fwd         :3d} rev      {rev     :3d} pos       {pos      :3d}")

                list_pos = 0
                # if chrom_num == 2: break

            if last_chrom_num != chrom_num: #todo: do every N million bp instead of whole chromosomes
                print(f"  new chrom {name} {chrom_num:03d} {last_list_pos:15,d}")
                chromosomes.append((name, seq_len)) #todo: ignore if fastq file

            last_chrom_num = chrom_num

        kmers[list_pos] = pos
        list_pos += 1

    if list_pos > 0:
        hist_k = process_kmers(kmers[:list_pos], header, frag_size=frag_size, debug=(kmer_len <= 5 and DEBUG) or debug)

        if hist is None: hist  = hist_k
        else:            hist += hist_k

        # print("hist_k     ", hist_k)
        # print("hist       ", hist)

        if (kmer_len <= 5 and DEBUG) or debug:
            print(f"  fwd          {fwd         :3d} rev      {rev     :3d} pos       {pos      :3d}")

    del kmers
    gc.collect()

    header.num_kmers   = num_kmers
    header.chromosomes = chromosomes
    # header.hist      = hist.tolist()
    # print(hist.tolist())

    print(f"project_name {header.project_name} kmer_len {header.kmer_len:15,d} num_kmers {header.num_kmers:15,d} kmer_size {header.kmer_size:15,d} max_size {header.max_size:15,d}")

    print("  indexing finished. creating header")

    header.write_metadata_index_tmp_file()

    hist_list = hist.tolist()
    assert header.hist == hist_list, f"hist_list\n{hist_list}\nheader.hist\n{header.hist}"

    print("  DONE")

    print("renaming")
    os.rename(header.index_tmp_file, header.index_file)

    print("done")

def read_fasta_index(
        project_name: str,
        input_file  : Union[str, None] = None,
        kmer_len    : Union[int, None] = None,
        index_file  : Union[str, None] = None,
        debug       : bool = False) -> None:

    header   = Header(project_name, input_file=input_file, kmer_len=kmer_len, index_file=index_file)

    header.read_metadata()
    print(header)

    print(f"project_name {header.project_name} kmer_len {header.kmer_len:15,d} bytes     {header.max_size//1024:15,d} Kb        {header.max_size//1024//1024:15,d} Mb       {header.max_size//1024//1024//1024:15,d} Gb")
    print(f"project_name {header.project_name} kmer_len {header.kmer_len:15,d} num_kmers {header.num_kmers     :15,d} kmer_size {header.kmer_size           :15,d} max_size {header.max_size                  :15,d}")

    header.check_data_index()
    print("OK")

    for npmm in header.get_array_from_index_file():
        # https://docs.python.org/3/library/mmap.html

        if (header.kmer_len <= 5 and DEBUG) or debug:
            for byte_pos in range(npmm.shape[0]):
                byte_val = npmm[byte_pos].item()
                print(byte_val, end=" ")
            print()

        npmm._mmap.close()
        del npmm

    # https://www.programiz.com/python-programming/methods/built-in/memoryview
    # for byte_pos, byte_val in enumerate(memoryview(mm[header.HEADER_LEN:header.HEADER_LEN+header.data_size])):

def run_test(overwrite: bool = False):
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
        input_file = f"{project_name}-{kmer_len:02d}.fasta.gz"
        # input_file = None # None for stdin

        index_file = f"{project_name}-{kmer_len:02d}.fasta.gz.{kmer_len:02d}.bin"

        print(f"project_name {project_name} kmer_len {kmer_len:15,d}")
        create_fasta_index(project_name, input_file, index_file, kmer_len, overwrite=overwrite)
        # read_fasta_index(project_name, index_file)

        print()

def main() -> None:
    if len(sys.argv) == 1:
        run_test(overwrite=False)

    else:
        input_file =     sys.argv[1]
        kmer_len        = int(sys.argv[2])

        project_name    = input_file

        buffer_size = 2**16
        # buffer_size = io.DEFAULT_BUFFER_SIZE

        print(f"project_name {project_name:s} input_file {input_file:s} kmer_len {kmer_len:15,d}")

        create_fasta_index(project_name, input_file, kmer_len, buffer_size=buffer_size, overwrite=True, debug=False)

        read_fasta_index(project_name, input_file=input_file, kmer_len=kmer_len)

    print()


if __name__ == "__main__":
    main()
