# kmer-index

## Install

`pip install -r requirements.txt`

## Run

`225M S_lycopersicum_chromosomes.4.00.fa.gz`

```bash
K=15

time pypy ./indexer.py S_lycopersicum_chromosomes.4.00.fa.bgz $K

F=S_lycopersicum_chromosomes.4.00.fa.bgz.$K.kin

bgzip -i -I $F.bgz.gzi -l 9 -c $F > $F.bgz

./gzireader.py $F.bgz.gzi

./merger.py $F.bgz $F.bgz
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

- `hist` are the number of kmers with the coverate equal to the position in the `array + 1`. `array[0] = 100` means 100 kmers have coverage of `1`.

- Maximum coverage is `255` with 1 byte per k-mer

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

## TODO

- ~~<https://numpy.org/doc/stable/reference/generated/numpy.memmap.html>~~
- ~~block compress~~
- ~~json header~~
- multithread per chromosome
- merge databases
- create matrix
