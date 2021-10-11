import os
import io
import json
import datetime
import socket
import hashlib

from typing import Callable, Generator, Tuple, Union, List, Dict, NewType, BinaryIO, TextIO

import numpy as np


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
