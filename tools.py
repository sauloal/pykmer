import os
import io
import json
import datetime
import socket
import hashlib

from collections import OrderedDict
from typing import Tuple, Union, List, NewType, BinaryIO, Dict, Any, Iterator

import numpy as np
import bgzip

Datetime    = NewType('Datetime'   , datetime.datetime)
Hist        = NewType('Hist'       , List[int])
Chromosomes = NewType('Chromosomes', List[Tuple[str, int]])
Metadata    = NewType('Chromosomes', Dict[str, Any])

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
            f"ela   time {self.time_ela_s  } val {self.val_last :15,d} speed {self.speed_ela  :15,d}\n"
            f"delta time {self.time_delta_s} val {self.val_delta:15,d} speed {self.speed_delta:15,d}"
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
        "input_file_name"    , "input_file_path"  ,
        "input_file_size"    , "input_file_ctime" , "input_file_cheksum" ,
        "output_file_size"   , "output_file_ctime", "output_file_cheksum",
        "num_kmers"          , "chromosomes"      ,
        "creation_time_start", "creation_time_end", "creation_duration"  , "creation_speed",
        "hostname"           , "checksum_script"  ,
        "hist"               ,
        "hist_sum"           , "hist_count"       , "hist_min"         , "hist_max",
        "vals_sum"           , "vals_count"       , "vals_min"         , "vals_max"
    ]
    IND_EXT  = 'kin'
    DESC_EXT = 'json'
    TMP      = 'tmp'

    def __init__(self,
            project_name   : str,
            input_file     : Union[str, None] = None,
            kmer_len       : Union[int, None] = None,
            index_file     : Union[str, None] = None,
            buffer_size    : int = io.DEFAULT_BUFFER_SIZE):

        self.project_name          :str         = project_name
        self.input_file_name       :str         = os.path.basename(input_file) if input_file else input_file
        self.input_file_path       :str         = os.path.abspath( input_file) if input_file else input_file
        self.kmer_len              :int         = kmer_len
        self._buffer_size          :int         = buffer_size

        self.input_file_size       :int         = None
        self.input_file_ctime      :float       = None
        self.input_file_cheksum    :str         = None

        self.output_file_size      :int         = None
        self.output_file_ctime     :float       = None
        self.output_file_cheksum   :str         = None

        self.num_kmers             :int         = None
        self.chromosomes           :Chromosomes = None

        self.timer                 :Timer       = Timer()
        self.creation_time_start   :str         = None
        self.creation_time_end     :str         = None
        self.creation_duration     :str         = None
        self.creation_speed        :int         = None

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

        if index_file is not None:
            self._parse_index_file_name(index_file)
            self.read_metadata()

        assert self.kmer_len
        assert self.kmer_len > 0
        assert self.kmer_len % 2 == 1

    @property
    def index_file(self) -> str:
        if os.path.exists(f"{self.index_file_root}.bgz"):
            return f"{self.index_file_root}.bgz"
        else:
            return self.index_file_root

    @property
    def index_file_root(self) -> str:
        return f"{self.input_file_path}.{self.kmer_len:02d}.{self.IND_EXT}"

    @property
    def index_tmp_file(self) -> str: return f"{self.index_file_root}.{self.TMP}"

    @property
    def metadata_file(self) -> str: return f"{self.index_file_root}.{self.DESC_EXT}"

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


    def _parse_index_file_name(self, index_file: str) -> None:
        index_file = index_file[:-4] if index_file.endswith('.bgz') else index_file
        ext_len    = ((2 + 1) + (len(self.IND_EXT) + 1))
        ext        = index_file[ext_len * -1:]
        # print(f"index_file {index_file}")
        # print(f"ext        {ext}")

        if self.input_file_name is None:
            self.input_file_name = index_file[:ext_len * -1]
            # print(f"input_file_name {self.input_file_name}")
        
        if self.kmer_len is None:
            self.kmer_len        = int(ext[1:3])
            # print(f"kmer_len {self.kmer_len}")

    def _get_mmap(self, fhd: BinaryIO, offset: int=0, mode: str="r+") -> Iterator[np.memmap]:
        nmmp = np.memmap(fhd, dtype=np.uint8, mode=mode, offset=offset, shape=(self.data_size,))
        yield nmmp
        del nmmp


    def update_stats(self, fhd: BinaryIO) -> None:
        print("updating stats")

        for arr in self.get_array_from_fhd(fhd):
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
        for fhd in self.open_index_file():
            self.update_stats(fhd)

    def update_stats_index_tmp_file(self) -> None:
        for fhd in self.open_index_tmp_file():
            self.update_stats(fhd)

    def update_metadata(self, index_file: str) -> None:
        print("updating metadata")

        self.input_file_size     = os.path.getsize(self.input_file_path)
        self.input_file_ctime    = os.path.getctime(self.input_file_path)
        self.input_file_cheksum  = gen_checksum(self.input_file_path)

        self.output_file_size    = os.path.getsize(index_file)
        self.output_file_ctime   = os.path.getctime(index_file)
        self.output_file_cheksum = gen_checksum(index_file)

        self.hostname            = socket.gethostname()
        self.checksum_script     = gen_checksum(os.path.basename(__file__))

        time_end                 = datetime.datetime.now()
        self.creation_time_start = str(self.timer.time_begin)
        self.creation_time_end   = str(time_end)
        self.creation_duration   = str(time_end - self.timer.time_begin)
        self.creation_speed      = self.timer.speed_ela


    def open_file(self, index_file: str, mode: str="r+b") -> Iterator[BinaryIO]:
        if index_file.endswith('.bgz'):
            with open(index_file, mode, buffering=self._buffer_size) as raw:
                with bgzip.BGZipReader(raw) as fhd:
                    yield fhd
        else:
            with open(index_file, mode, buffering=self._buffer_size) as fhd:
                yield fhd

    def open_index_file(self, mode: str="r+b") -> Iterator[BinaryIO]:
        return self.open_file(self.index_file, mode=mode)

    def open_index_tmp_file(self, mode: str="r+b") -> Iterator[BinaryIO]:
        return self.open_file(self.index_tmp_file, mode=mode)


    def _init_clean(self, overwrite: bool=False) -> None:
        if os.path.exists(self.index_file):
            if overwrite:
                os.remove(self.index_file)
            else:
                raise ValueError(f"file {self.index_file} already exists and overwritting disabled")

        if os.path.exists(self.index_file_root):
            if overwrite:
                os.remove(self.index_file_root)
            else:
                raise ValueError(f"file {self.index_file_root} already exists and overwritting disabled")

        if os.path.exists(self.metadata_file):
            os.remove(self.metadata_file)

        if os.path.exists(self.index_tmp_file):
            os.remove(self.index_tmp_file)

    def init_file(self, index_file: str, mode: str="r+b") -> None:
        # print("opening")
        if not os.path.exists(index_file):
            with open(index_file, 'w'):
                pass

        for fhd in self.open_file(index_file, mode=mode):
            fhd.seek(self.max_size - 1)
            fhd.write(b'\0')
        # print(f"{f.tell():,d} bytes {f.tell()//1024:,d} Kb {f.tell()//1024//1024:,d} Mb {f.tell()//1024//1024//1024:,d} Gb")

    def init_index_file(self, overwrite: bool=False, mode: str="r+b") -> None:
        self._init_clean(self, overwrite=overwrite)
        return self.init_file(self.index_file, mode=mode)

    def init_index_tmp_file(self, overwrite: bool=False, mode: str="r+b") -> None:
        self._init_clean(overwrite=overwrite)
        return self.init_file(self.index_tmp_file, mode=mode)


    def get_array_from_fhd(self, fhd: BinaryIO, mode: str="r+") -> Iterator[np.memmap]:
        for npmm in self._get_mmap(fhd, offset=0, mode=mode):
            yield npmm

    def get_array_from_index_file(self, fhd_mode: str="r+b", mm_mode: str="r+") -> Iterator[np.memmap]:
        for fhd in self.open_index_file(mode=fhd_mode):
            return self.get_array_from_fhd(fhd, mode=mm_mode)

    def get_array_from_index_tmp_file(self, fhd_mode: str="r+b", mm_mode: str="r+") -> Iterator[np.memmap]:
        for fhd in self.open_index_tmp_file(mode=fhd_mode):
            return self.get_array_from_fhd(fhd, mode=mm_mode)


    def write_metadata_file(self, index_file: str) -> None:
        assert self.num_kmers
        assert self.chromosomes

        self.update_metadata(index_file)
        for fhd in self.open_file(index_file):
            self.update_stats(fhd)

        header_data = {k: getattr(self,k) for k in self.HEADER_DATA}

        with open(self.metadata_file, 'wt') as fhm:
            json.dump(header_data, fhm, indent=1, sort_keys=1)

    def write_metadata_index_file(self) -> None:
        self.write_metadata_file(self.index_file)

    def write_metadata_index_tmp_file(self) -> None:
        self.write_metadata_file(self.index_tmp_file)


    def read_metadata(self) -> None:
        with open(self.metadata_file, 'rt') as fhd:
            header_data = json.load(fhd)

        # print(header_data)

        for k in self.HEADER_DATA:
            v = header_data[k]
            # print(f"{k:20s}: {str(v)[:50]}")
            setattr(self, k, v)


    def check_data(self, fhd: BinaryIO) -> None:
        self.read_metadata()

        other = self.__class__(self.project_name, input_file=self.input_file_path, kmer_len=self.kmer_len)
        other.read_metadata()
        assert self.project_name     == other.project_name
        assert self.input_file_name  == other.input_file_name
        assert self.input_file_path  == other.input_file_path
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

    def __iter__(self) -> Iterator[int]:
        for fhd in self.open_index_file():
            c = fhd.read(1)
            while c:
                yield c[0]
                c = fhd.read(1)

    def calculate_distance(self, other: "Header", min_count: int=1):
        s_count = 0
        o_count = 0
        c_count = 0
        for pos, (s_char, o_char) in enumerate(zip(self, other)):
            s_valid = s_char >= min_count
            o_valid = o_char >= min_count
            s_count += 1 if s_valid             else 0
            o_count += 1 if o_valid             else 0
            c_count += 1 if s_valid and o_valid else 0
            # print(f"{pos:15,d} {s_char:3d} {o_char:3d} {s_count:15,d} {o_count:15,d} {c_count:15,d}")
        return s_count, o_count, c_count


    def to_dict(self) -> Dict[str, Any]:
        data = OrderedDict()
        for k in ["file_ver", "kmer_len", "kmer_size", "data_size", "max_size"] + self.HEADER_DATA:
            v = getattr(self, k)
            data[k] = v
        return data

    def to_json(self, indent: int=1, sort_keys: bool=True) -> str:
        return json.dumps(self.to_dict(), indent=indent, sort_keys=sort_keys)


    def __str__(self) -> str:
        res = []
        for k, v in self.to_dict().items():
            if isinstance(v, int):
                res.append(f"{k:20s}: {v:15,d}")
            else:
                res.append(f"{k:20s}: {str(v)[:50]}")
        return "\n".join(res) + "\n"

    def __repr__(self) -> str:
        return str(self)

def gen_checksum(filename: str, chunk_size: int=2**16) -> str:
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
