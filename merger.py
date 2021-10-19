#!/usr/bin/env python3

import os
import sys

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


def merge(indexes: List[str], min_count: int=1) -> Tuple[DataType, ResultType]:
    data: DataType = [None] * len(indexes)
    kmer_len = None
    for k, kin in enumerate(indexes):
        print(f"verifying {kin}")
        assert kin.endswith(EXTS) , f"all files must be .{Header.IND_EXT}[.bgz]: {kin}"
        assert os.path.exists(kin), f"all files must exist: {kin}"

        desc = kin[:-1*(len(Header.COMP_EXT)+1)] if kin.endswith('.'+Header.COMP_EXT) else kin
        desc = f"{desc}.{Header.DESC_EXT}"

        assert os.path.exists(desc), f"all .{Header.IND_EXT}[.{Header.COMP_EXT}] files must have a associated .{Header.IND_EXT}.{Header.DESC_EXT}: {desc}"

        header = Header(kin, index_file=kin)

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

    res:ResultType = OrderedDict()
    for k in range(len(data)-1):
        k_data:Metadata = data[k]
        k_header:Header = k_data["header"]
        print("comparing", k_data["index_file"])
        for l in range(k+1, len(data)):
            l_data:Metadata = data[l]
            print(" versus", l_data["index_file"])
            sys.stdout.flush()
            l_header:Header = l_data["header"]
            k_header.calculate_distance(l_header)
            k_count, l_count, s_count = k_header.calculate_distance(l_header, min_count=min_count)
            print("   res", l_data["index_file"], k_count, l_count, s_count)
            res[(k,l)] = (k_count, l_count, s_count)

    return data, res


def main() -> None:
    indexes = sys.argv[1:]

    if len(indexes) <= 1:
        print("needs at least 2 files")
        sys.exit(1)

    indexes.sort()
    # print(indexes)

    merge(indexes)

if __name__ == "__main__":
    main()