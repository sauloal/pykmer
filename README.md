# kmer-index

## Install

`pip install -r requirements.txt`

`apt-get install xvfb`

## Run

`225M S_lycopersicum_chromosomes.4.00.fa.gz`

```bash
K=15

time pypy ./indexer.py S_lycopersicum_chromosomes.4.00.fa.bgz $K

F=S_lycopersicum_chromosomes.4.00.fa.bgz.$K.kin

bgzip -i -I $F.bgz.gzi -l 9 -c $F > $F.bgz

./gzireader.py $F.bgz.gzi

./merger.py $F.bgz $F.bgz

pypy ./merger.py data/*.15.kin.bgz

./calculate_distance.sh merged.001-255.kma

./calculate_distance.sh merged.001-050.kma
```

## Benchmark

|     K  |  Real Time  |  Speed         |  Kin Size  |  bgzip  |    xz  |   7z  |  lz4  |
|--------|-------------|----------------|------------|---------|--------|-------|-------|
|     3  |     16m21s  |  797,621 bp/s  |        2K  |    69b  |        |       |       |
|        |     16m 6s  |  809,751 bp/s  |        3K  |    94b  |        |       |       |
|     7  |     16m27s  |  787,715 bp/s  |       16K  |   287b  |        |       |       |
|     9  |     18m28s  |  706,750 bp/s  |      256K  |    11K  |        |       |       |
|    11  |     18m45s  |  702,199 bp/s  |        4M  |     2M  |        |       |       |
|    13  |     19m54s  |  677,203 bp/s  |       64M  |    26M  |        |       |       |
|    15  |     27m50s  |  503,287 bp/s  |        1G  |   156M  |  124M  |  130M |  180M |
|    17  |     111m5s  |  128,452 bp/s  |       17G  |   574M  |        |       |       |
|  X 19  |             |                |      257G  |         |        |       |       |
|  X 21  |             |                |        4T  |         |        |       |       |

## Merge

```bash
time pypy ./merger.py merged data/*.15.kin.bgz
# merged.001-255.kma

time pypy ./merger.py merged data/*.15.kin.bgz --max-count=50
# merged.001-050.kma
```

```text
comparing data/Solanum_tuberosum_PGSC_DM_v4.03_pseudomolecules.fa.bgz.15.kin.bgz
 versus data/Vitis_vinifera_Genoscope_12X_2010_02_12_chr.fa.bgz.15.kin.bgz
               1               1/  1,073,741,824 (  0.00%)               2     100,000,001/  1,073,741,824 (  9.31%)               3     200,000,001/  1,073,741,824 ( 18.63%)
               4     300,000,001/  1,073,741,824 ( 27.94%)               5     400,000,001/  1,073,741,824 ( 37.25%)               6     500,000,001/  1,073,741,824 ( 46.57%)
               7     600,000,001/  1,073,741,824 ( 55.88%)               8     700,000,001/  1,073,741,824 ( 65.19%)               9     800,000,001/  1,073,741,824 ( 74.51%)
              10     900,000,001/  1,073,741,824 ( 83.82%)              11   1,000,000,001/  1,073,741,824 ( 93.13%)
   matrix Total #1     172,022,482 Total #2     145,297,478 Shared      84,710,204
saving merged.kma.json
saving merged.kma

real    333m56.907s 6h
user    303m45.443s
sys      29m45.602s

real    351m20.043s
user    1256m13.652s
sys     113m55.620s
```

## Output

`.kin`

The file contents are:

```text
uint8_t coverage
```

Each position correspond to a `k-mer` encoded in binary.

```text
AAAA = 0000 = 0
              00 00 00 00

ACGT = 0123 = 0*256 + 1*64 + 2*16 + 3*4 =   0 +  64 + 32 + 12 = 108
              00      01     10     11                          00 01 10 11

TGCA = 3210 = 3*256 + 2*64 + 1*16 + 0*4 = 768 + 128 + 16 +  0 = 912
              11      10     01     00                          11 10 01 00

CCAA = 2200 = 2*256 + 2*64 + 0*16 + 0*4 = 512 + 128 +  0 +  0 = 640 = 40
AACC = 0022 = 0*256 + 0*64 + 2*16 + 2*4 =   0 +   0 + 32 +  8 =  40 = 40
              00      00     10     10                          00 00 10 10
```

`.kin.json`

```json
{
    "checksum_script": "93d6365c05ea8ca7fd302b95a73657d47ddc467f2eab264df45652c4ae28d344",
    "chromosomes": [
        [ "SL4.0ch00", 9643250 ],
        [ "SL4.0ch12", 66688036]
    ],
    "creation_duration": "0:27:14.185148",
    "creation_time_end": "2021-10-10 21:39:27.735455",
    "creation_time_start": "2021-10-10 21:12:13.550307",
    "hist": [ 88567750, 35330720, 753, 120840 ],
    "hist_count": 255,
    "hist_max": 88567750,
    "hist_min": 731,
    "hist_sum": 186366572,
    "hostname": "SAULOACER",
    "input_file_cheksum": "365255b77847ebbbbbabe0403b214bbf5e0e7a49ba1801f268c77876174e5184",
    "input_file_ctime": 1632606653.5075698,
    "input_file_name": "S_lycopersicum_chromosomes.4.00.fa.gz",
    "input_file_size": 235901169,
    "kmer_len": 15,
    "num_kmers": 782469030,
    "output_file_cheksum": "ac535097ca9d1cca82a05b1822e6348982810748056c776ece6fe18ec78be2d5",
    "output_file_ctime": 1633894762.362578,
    "output_file_size": 1073741824,
    "project_name": "S_lycopersicum_chromosomes.4.00.fa.gz",
    "vals_count": 186366572,
    "vals_max": 255,
    "vals_min": 0,
    "vals_sum": 730970717
}
```

Where:

- `chromosomes` are the chromosomes and their lengths.

- `hist` are the number of kmers with the coverate equal to the position in
the `array + 1`. `array[0] = 100` means 100 kmers have coverage of `1`.

- Maximum coverage is `255` with 1 byte per k-mer

```text
-rw-r--r-- 1 saulo saulo 9.4K Oct 21 05:38 merged.kma
-rw-r--r-- 1 saulo saulo 305M Oct 21 05:38 merged.kma.json
```

`.kma`

Numpy [`Savez`](https://numpy.org/doc/stable/reference/generated/numpy.savez_compressed.html#numpy.savez_compressed)
with matrix as `matrix` key.

`.kma.json`

```json
{
 "max_count": 255,
 "min_count": 1,
 "project_name": "merged",
 "data": [
  {
   "description_file": "data/Arabidopsis_thaliana_10.fa.bgz.15.kin.json",
   "header": {
    "checksum_script": "17a685dbeeef334e1dcb6e5f38f84715a2c7cec5028454af071451e9a52bcadc",
    "chromosomes": [
     [
      "Chr1 CHROMOSOME dumped from ADB: Jun/20/09 14:53; last updated: 2009-02-02",
      30427671
     ],
    ],
    "creation_duration": "0:03:07.736857",
    "creation_speed": 829989,
    "creation_time_end": "2021-10-16 22:36:27.578499",
    "creation_time_start": "2021-10-16 22:33:19.841642",
    "data_size": 1073741824,
    "file_ver": "KMER001",
    "flush_every": 500000000,
    "frag_size": 357914000,
    "hist": [
     55881827,
     954
    ],
    "hist_count": 255,
    "hist_max": 55881827,
    "hist_min": 2,
    "hist_sum": 75864811,
    "hostname": "SAULOACER",
    "input_file_cheksum": "1f9a12cce0fa57be18ffbec8d46dfe224dd56336638daabc6458417e3b754482",
    "input_file_ctime": 1634346116.4379272,
    "input_file_name": "Arabidopsis_thaliana_10.fa.bgz",
    "input_file_path": "/home/saulo/Programs/kmer-index/data/Arabidopsis_thaliana_10.fa.bgz",
    "input_file_size": 36390576,
    "kmer_len": 15,
    "kmer_size": 1073741824,
    "max_size": 1073741824,
    "num_kmers": 119478452,
    "output_file_cheksum": "4208a6f7e1d6240ea4e2de81d070d2a719d1f74a3c20767da350f8def4f7d497",
    "output_file_ctime": 1634416582.7479272,
    "output_file_size": 1073741824,
    "project_name": "data/Arabidopsis_thaliana_10.fa.bgz",
    "vals_count": 75864811,
    "vals_max": 255,
    "vals_min": 0,
    "vals_sum": 119184917
   },
   "index_file": "data/Arabidopsis_thaliana_10.fa.bgz.15.kin.bgz",
   "pos": 0
  }
 ]
}
```

## Analysis

### Convert

`./calculate_distance.sh merged.001-255.kma`

`./calculate_distance.sh merged.001-050.kma`

or

`xvfb-run python3 $PWD/calculate_distance.py merged.kma`

### Cluster Output

| Filename                                    | Description                   |
|---------------------------------------------|-------------------------------|
| merged.kma                                  | Input Matrix file             |
| merged.kma.json                             | Input Matrix Header           |
| merged.kma.names.tsv                        | [Optional] Input Sample names |
| merged.kma.dist.jaccard.mat.condensed.np    | Condensed matrix np           |
| merged.kma.dist.jaccard.mat.condensed.txt   | Condensed matrix tsv          |
| merged.kma.dist.jaccard.mat.redundant.lsmat | Redundant matrix tsv          |
| merged.kma.dist.jaccard.mat.redundant.np    | Redundant matrix np           |
| merged.kma.dist.jaccard.npz                 | Distance matrix npz           |
| merged.kma.dist.jaccard.newick              | Tree - Newick                 |
| merged.kma.dist.jaccard.png                 | Tree - PNG                    |
| merged.kma.dist.jaccard.tree                | Tree - ASCII Art              |

## Development notes

### Profiling

```bash
time pypy -m cProfile -s cumulative ./indexer.py
```

### Compress

```bash
for BIN in *.kin; do
    if [[ ! -f "$BIN.bgz" ]]; then
        bgzip -i -I $BIN.bgz.gzi -l 9 -c $BIN > $BIN.bgz
    fi
done
```

### Benchmark Calculation

```bash
time pypy ./indexer.py S_lycopersicum_chromosomes.4.00.fa.gz 15

K=15
F=S_lycopersicum_chromosomes.4.00.fa.gz.$K.kin
bgzip -i -I $F.bgz.gzi -l 9 -c $F > $F.bgz

```

## TODO

- ~~<https://numpy.org/doc/stable/reference/generated/numpy.memmap.html>~~
- ~~block compress~~
- ~~json header~~
- ~~merge databases~~
- ~~create matrix~~
- multithread per chromosome
- multithread merge

### ~~Pipe~~

```bash
K=15
time (gunzip -c -k example-$K.fasta.gz | pypy ./indexer.py )
```

```text
pypy K=15
pygzip  PIPE        pygz        STDIN       Popen
real    56m11.378s  33m48.443s  35m50.502s  66m 5.744s
```
