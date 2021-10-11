#!/usr/bin/env python3

import os
import sys
import io
import gzip

from typing import Iterator, Tuple, Union, List, TextIO

import gc

import numpy as np

from tools import Timer, Header

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

def parse_fasta(fhd: TextIO, print_every: int = 25_000_000) -> Iterator[Tuple[str, str, int]]:
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

def read_fasta(input_file: Union[str, None]) -> Iterator[Tuple[str, str, int]]:
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

def gen_kmers(input_file: str, kmer_len: int) -> Iterator[Tuple[int, str, int, int, int]]:
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

            for fhd in header.open_index_tmp_file():
                npmm                    = header.get_array_from_fhd(fhd).next()
                npmm_frag               = npmm[frag_from:frag_to]
                hist_before             = np.histogram(npmm_frag, bins=header.max_val, range=(1,header.max_val))[0]
                npmm[frag_from:frag_to] = npmm_frag + np.minimum(header.max_val - npmm_frag, vals_frag)
                hist_after              = np.histogram(npmm_frag, bins=header.max_val, range=(1,header.max_val))[0]

                # if hist is None: hist_before[0] = 2 * header.max_val + 2
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
