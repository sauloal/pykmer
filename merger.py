#!/usr/bin/env python3

import os
import sys
import time
import argparse

import pathlib
from pathlib import Path
from typing import Dict, List, Tuple, NewType

#https://stackoverflow.com/questions/18478287/making-object-json-serializable-with-regular-encoder/18561055#18561055

import numpy as np

import json
from   json import JSONEncoder

from multiprocessing.pool import AsyncResult
from multiprocessing import Pool, TimeoutError


def _default(self, obj):
    if hasattr(obj.__class__, "to_dict"):
        return getattr(obj.__class__, "to_dict", _default.default)(obj)
    elif isinstance(obj, (pathlib.Path, pathlib.PosixPath, pathlib.PurePath, pathlib.PurePosixPath, pathlib.PureWindowsPath, pathlib.WindowsPath)):
        return str(obj)

_default.default    = JSONEncoder().default
JSONEncoder.default = _default


from tools import Header, Metadata

ResultType = NewType('ResultType', Dict[Tuple[int,int],Tuple[int,int,int]])
DataType   = NewType('DataType'  , List[Metadata])

EXTS = (
    '.'+Header.IND_EXT,
    '.'+Header.IND_EXT+'.'+Header.COMP_EXT,
    '.kma',
    '.kma.'+Header.COMP_EXT
)

DEFAULT_MIN_COUNT   = Header.DEFAULT_MIN_COUNT
DEFAULT_MAX_COUNT   = Header.DEFAULT_MAX_COUNT
DEFAULT_BUFFER_SIZE = Header.DEFAULT_BUFFER_SIZE
DEFAULT_BLOCK_SIZE  = Header.DEFAULT_BLOCK_SIZE
DEFAULT_THREADS     = 4

parser = argparse.ArgumentParser(description='Merge kmer databases.')
parser.add_argument('Project_Name' , metavar='P', type=str ,                                         help='Project name')
parser.add_argument('Kmer_1'       , metavar='K', type=Path,                              nargs=1  , help='List of kin files')
parser.add_argument('Kmer_N'       , metavar='K', type=Path,                              nargs='+', help='List of kin files')
parser.add_argument('--min-count'  ,              type=int , default=DEFAULT_MIN_COUNT  , nargs='?', help=f'Minimum Kmer Count [{DEFAULT_MIN_COUNT}]')
parser.add_argument('--max-count'  ,              type=int , default=DEFAULT_MAX_COUNT  , nargs='?', help=f'Maximum Kmer Count [{DEFAULT_MAX_COUNT}]')
parser.add_argument('--buffer-size',              type=int , default=DEFAULT_BUFFER_SIZE, nargs='?', help=f'Buffer size [{DEFAULT_BUFFER_SIZE}]')
parser.add_argument('--block-size' ,              type=int , default=DEFAULT_BLOCK_SIZE , nargs='?', help=f'Block size [{DEFAULT_BLOCK_SIZE}]')
parser.add_argument('--threads'    ,              type=int , default=DEFAULT_THREADS    , nargs='?', help=f'Threads [{DEFAULT_THREADS}]')


def calculate_distance(
        k_index_file: str,
        l_index_file: str,
        min_count   : int=DEFAULT_MIN_COUNT,
        max_count   : int=DEFAULT_MAX_COUNT,
        buffer_size : int=DEFAULT_BUFFER_SIZE,
        block_size  : int=DEFAULT_BLOCK_SIZE
    ) -> Tuple[int,int,int]:
    k_index_file = str(k_index_file)
    l_index_file = str(l_index_file)
    
    k_header     = Header(k_index_file, index_file=k_index_file, buffer_size=buffer_size)
    l_header     = Header(l_index_file, index_file=l_index_file, buffer_size=buffer_size)

    k_count, l_count, s_count = k_header.calculate_distance(l_header, min_count=min_count, max_count=max_count, block_size=block_size, threading=True)

    return k_count, l_count, s_count

def merge(
        project_name: str,
        indexes     : List[Path],
        min_count   : int=DEFAULT_MIN_COUNT,
        max_count   : int=DEFAULT_MAX_COUNT,
        buffer_size : int=DEFAULT_BUFFER_SIZE,
        block_size  : int=DEFAULT_BLOCK_SIZE,
        threads     : int=DEFAULT_THREADS
    ) -> Tuple[DataType, ResultType]:

    assert min_count    >=   1 
    assert max_count    <= 255
    assert buffer_size  >    0
    assert block_size   >    0
    assert len(indexes) >    0

    outfile = Path(f"{project_name}.{min_count:03d}-{max_count:03d}.kma")

    assert not Path(project_name).exists(), f"project name ({project_name}) is a file. maybe forgot to pass project name as first argument?"
    assert not outfile.exists()           , f"project output file ({outfile}) already exists. not overwriting."

    indexes = [Path(p) for p in indexes]
    assert all([i.exists() for i in indexes])

    data: DataType = [None] * len(indexes)
    kmer_len       = None
    for k, kin in enumerate(indexes):
        print(f"verifying {kin}")
        kins = str(kin)
        assert kins.endswith(EXTS), f"all files must be .{Header.IND_EXT}[.bgz]: {kin}"
        assert os.path.exists(kin), f"all files must exist: {kin}"

        desc = kins[:-1*(len(Header.COMP_EXT)+1)] if kins.endswith('.'+Header.COMP_EXT) else kin
        desc = Path(f"{desc}.{Header.DESC_EXT}")

        assert desc.exists(), f"all .{Header.IND_EXT}[.{Header.COMP_EXT}] files must have a associated .{Header.IND_EXT}.{Header.DESC_EXT}: {desc}"

        header = Header(kins, index_file=kins, buffer_size=buffer_size)

        if kmer_len is None:
            kmer_len = header.kmer_len

        assert header.kmer_len == kmer_len, f"kmer_length differs. expected {kmer_len}, got {header.kmer_len}"

        data[k] = {
            "pos"             : k,
            "index_file"      : kin,
            "description_file": desc,
            "header"          : header
        }

    print()
    # print(json.dumps(data, indent=1))
    # exit(0)

    matrix:ResultType = [None] * len(data)
    matrix            = np.ndarray(shape=(len(data),len(data),3), dtype=np.uint64)
    with Pool(processes=threads) as pool:
        tasks: Dict[Tuple[int,int], Tuple[bool,AsyncResult]] = {}
        for k in range(len(data)-1):
            k_data:Metadata = data[k]
            # k_header:Header = k_data["header"]
            
            print("comparing", k_data["index_file"])
            for l in range(k+1, len(data)):
                l_data:Metadata = data[l]

                print(" versus", l_data["index_file"])
                sys.stdout.flush()

                # l_header:Header = l_data["header"]

                task = pool.apply_async(calculate_distance, (k_data["index_file"], l_data["index_file"]), {"buffer_size": buffer_size, "min_count": min_count, "max_count": max_count, "block_size": block_size})
                tasks[(k,l)] = [True, task]

                # k_count, l_count, s_count = k_header.calculate_distance(l_header, min_count=min_count, max_count=max_count, block_size=block_size)

                # print(f"   matrix Total #1 {k_count:15,d} Total #2 {l_count:15,d} Shared {s_count:15,d}")
                # matrix[k,l,:] = (k_count, l_count, s_count)
                # matrix[l,k,:] = (l_count, k_count, s_count)

        pool.close()

        while len(tasks) > 0:
            print(f"waiting for {len(tasks):7,d} tasks")
            sys.stdout.flush()
            for key, (running, task) in tasks.items():
                # print("  checking", key)
                try:
                    k_count, l_count, s_count = task.get(timeout=0)
                except TimeoutError:
                    continue
                (k,l) = key
                print(f"   matrix Total #{k:3d} {k_count:15,d} Total #{l:3d} {l_count:15,d} Shared {s_count:15,d}")
                sys.stdout.flush()
                matrix[k,l,:] = (k_count, l_count, s_count)
                matrix[l,k,:] = (l_count, k_count, s_count)
                tasks[key][0] = False
            tasks = {k:v for k,v in tasks.items() if v[0]}
            # print("  sleeping", key)
            # sys.stdout.flush()
            time.sleep(10)

        print(f"\nfinished waiting. joining\n")
        pool.join()
        print(f"\nfinished waiting. joined\n")

    for v in data:
        v["header"] = v["header"].to_dict(lean=True)

    output = {
        "project_name": project_name,
        "min_count"   : min_count,
        "max_count"   : max_count,
        "data"        : data,
    }

    outfile_json     = Path(f"{outfile}.json")
    outfile_json_tmp = Path(f"{outfile_json}.tmp")
    print(f"saving {outfile_json}")
    with outfile_json_tmp.open(mode="wt") as fhd:
        json.dump(output, fhd, sort_keys=True, indent=1)
    outfile_json_tmp.rename(outfile_json)

    print(f"saving {outfile}")
    outfile_tmp = Path(f"{outfile}.tmp")
    with outfile_tmp.open(mode="wb") as fhd:
        np.savez_compressed(fhd, matrix=matrix)
    outfile_tmp.rename(outfile)

    return data, matrix


def main() -> None:
    args                     = parser.parse_args()

    project_name :str        = args.Project_Name
    indexes      :List[Path] = args.Kmer_1 + args.Kmer_N
    min_count    :int        = args.min_count
    max_count    :int        = args.max_count
    buffer_size  :int        = args.buffer_size
    block_size   :int        = args.block_size
    threads      :int        = args.threads

    if len(indexes) <= 1:
        print("needs at least 2 files")
        sys.exit(1)

    indexes.sort()
    # print(indexes)

    merge(
        project_name,
        indexes,
        min_count   = min_count,
        max_count   = max_count,
        buffer_size = buffer_size,
        block_size  = block_size,
        threads     = threads
    )

if __name__ == "__main__":
    main()