#!/usr/bin/env python3

import os
import sys
import gzip

from collections import OrderedDict
from typing import Dict, List, Tuple, NewType

#https://stackoverflow.com/questions/18478287/making-object-json-serializable-with-regular-encoder/18561055#18561055

import json
from json import JSONEncoder

def _default(self, obj):
    # return getattr(obj.__class__, "to_json", _default.default)(obj)
    return getattr(obj.__class__, "to_dict", _default.default)(obj)

_default.default = JSONEncoder().default
JSONEncoder.default = _default


from tools import Timer, Header, Metadata

ResultType = NewType('ResultType', Dict[Tuple[int,int],Tuple[int,int,int]])
DataType   = NewType('DataType'  , List[Metadata])

EXTS = (
    '.'+Header.IND_EXT,
    '.'+Header.IND_EXT+'.'+Header.COMP_EXT,
    '.kma',
    '.kma.'+Header.COMP_EXT
)


def merge(project_name: str, indexes: List[str], min_count: int=1, max_count: int=255, buffer_size: int=Header.DEFAULT_BUFFER_SIZE, blk_size: int=100_000_000) -> Tuple[DataType, ResultType]:
    outfile = f"{project_name}.kma"

    assert not os.path.exists(project_name), f"project name ({project_name}) is a file. maybe forgot to pass project name as first argument?"
    assert not os.path.exists(outfile     ), f"project output file ({outfile}) already exists. not overwriting."

    data: DataType = [None] * len(indexes)
    kmer_len       = None
    for k, kin in enumerate(indexes):
        print(f"verifying {kin}")
        assert kin.endswith(EXTS) , f"all files must be .{Header.IND_EXT}[.bgz]: {kin}"
        assert os.path.exists(kin), f"all files must exist: {kin}"

        desc = kin[:-1*(len(Header.COMP_EXT)+1)] if kin.endswith('.'+Header.COMP_EXT) else kin
        desc = f"{desc}.{Header.DESC_EXT}"

        assert os.path.exists(desc), f"all .{Header.IND_EXT}[.{Header.COMP_EXT}] files must have a associated .{Header.IND_EXT}.{Header.DESC_EXT}: {desc}"

        header = Header(kin, index_file=kin, buffer_size=buffer_size)

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
    # exit(0)
    # print(json.dumps(data, indent=1))

    matrix:ResultType = [None] * len(data)
    for k in range(len(data)-1):
        k_data:Metadata = data[k]
        k_header:Header = k_data["header"]
        print("comparing", k_data["index_file"])
        matrix[k] = [None] * len(data)
        for l in range(k+1, len(data)):
            l_data:Metadata = data[l]
            print(" versus", l_data["index_file"])
            sys.stdout.flush()
            l_header:Header = l_data["header"]
            # k_header.calculate_distance(l_header)
            k_count, l_count, s_count = k_header.calculate_distance(l_header, min_count=min_count, max_count=max_count, blk_size=blk_size)
            print(f"   matrix Total #1 {k_count:15,d} Total #2 {l_count:15,d} Shared {s_count:15,d}")
            matrix[k][l] = (k_count, l_count, s_count)


    output = {
        "project_name": project_name,
        "min_count"   : min_count,
        "max_count"   : max_count,
        "data"        : data,
        "matrix"      : matrix
    }

    print(f"saving {project_name}.kma")
    with gzip.open(outfile+'.tmp', "wt") as fhd:
        json.dump(output, fhd, sort_keys=True, indent=1)
    os.rename(outfile+'.tmp', outfile)

    return data, matrix


def main() -> None:
    project_name = sys.argv[1]
    indexes      = sys.argv[2:]

    if len(indexes) <= 1:
        print("needs at least 2 files")
        sys.exit(1)

    indexes.sort()
    # print(indexes)

    merge(project_name, indexes, buffer_size=Header.DEFAULT_BUFFER_SIZE)

if __name__ == "__main__":
    main()