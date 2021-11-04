# KWIP

<https://kwip.readthedocs.io/en/latest/kwip.html#overview>

<https://github.com/kdm9/kWIP>

<https://github.com/dib-lab/khmer>

## Help

```text
./kwip 
kwip version 0.2.0

USAGE: ./kwip [options] hashes

OPTIONS:
-t, --threads       Number of threads to utilise. [default N_CPUS]
-k, --kernel        Output file for the kernel matrix. [default None]
-d, --distance      Output file for the distance matrix. [default stdout]
-U, --unweighted    Use the unweighted inner proudct kernel. [default off]
-w, --weights       Bin weight vector file (input, or output w/ -C).
-C, --calc-weights  Calculate only the bin weight vector, not kernel matrix.
-h, --help          Print this help message.
-V, --version       Print the version string.
-v, --verbose       Increase verbosity. May or may not acutally do anything.
-q, --quiet         Execute silently but for errors.

Each sample's oxli Countgraph should be specified after arguments:
./kwip [options] sample1.ct sample2.ct ... sampleN.ct
```

## Install

### Install Khmer

```bash
pip install khmer
```

#### load into counting

```bash
load-into-counting.py -h

|| This is the script load-into-counting.py in khmer.
|| You are running khmer version 2.1.1
|| You are also using screed version 1.0.5
||
|| If you use this script in a publication, please cite EACH of the following:
||
||   * MR Crusoe et al., 2015. http://dx.doi.org/10.12688/f1000research.6924.1
||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
||
|| Please see http://khmer.readthedocs.io/en/latest/citations.html for details.

usage: load-into-counting.py [--version] [--info] [-h] [-k KSIZE] [-U UNIQUE_KMERS] [--fp-rate FP_RATE] [-M MAX_MEMORY_USAGE]
                             [--small-count] [-T THREADS] [-b] [-s FORMAT] [-f] [-q]
                             output_countgraph_filename input_sequence_filename [input_sequence_filename ...]

Build a k-mer countgraph from the given sequences.

positional arguments:
  output_countgraph_filename
                        The name of the file to write the k-mer countgraph to.
  input_sequence_filename
                        The names of one or more FAST[AQ] input sequence files.

optional arguments:
  --version             show program's version number and exit
  --info                print citation information
  -h, --help            show this help message and exit
  -k KSIZE, --ksize KSIZE
                        k-mer size to use (default: 32)
  -U UNIQUE_KMERS, --unique-kmers UNIQUE_KMERS
                        approximate number of unique kmers in the input set (default: 0)
  --fp-rate FP_RATE     Override the automatic FP rate setting for the current script (default: None)
  -M MAX_MEMORY_USAGE, --max-memory-usage MAX_MEMORY_USAGE
                        maximum amount of memory to use for data structure (default: None)
  --small-count         Reduce memory usage by using a smaller counter for individual kmers. (default: False)
  -T THREADS, --threads THREADS
                        Number of simultaneous threads to execute (default: 1)
  -b, --no-bigcount     The default behaviour is to count past 255 using bigcount. This flag turns bigcount off, limiting counts to 255.
                        (default: True)
  -s FORMAT, --summary-info FORMAT
                        What format should the machine readable run summary be in? (`json` or `tsv`, disabled by default) (default: None)
  -f, --force           Overwrite output file if it exists (default: False)
  -q, --quiet

Note: with `-b`/`--no-bigcount` the output will be the exact size of the k-mer
countgraph and this script will use a constant amount of memory. In exchange
k-mer counts will stop at 255. The memory usage of this script with `-b` will
be about 1.15x the product of the `-x` and `-N` numbers.

Example:

    load-into-counting.py -k 20 -x 5e7 out data/100k-filtered.fa

Multiple threads can be used to accelerate the process, if you have extra cores
to spare.

Example:

    load-into-counting.py -k 20 -x 5e7 -T 4 out data/100k-filtered.fa
```

Example

```bash
load-into-counting.py \
    --ksize 21 \
    -x 5e7 \
    --threads 4 \
    --no-bigcount \
    --summary-info json
    out
    in.fa
```

#### Load Graph

```text
load-graph.py -h

|| This is the script load-graph.py in khmer.
|| You are running khmer version 2.1.1
|| You are also using screed version 1.0.5
||
|| If you use this script in a publication, please cite EACH of the following:
||
||   * MR Crusoe et al., 2015. http://dx.doi.org/10.12688/f1000research.6924.1
||   * J Pell et al., http://dx.doi.org/10.1073/pnas.1121464109
||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
||
|| Please see http://khmer.readthedocs.io/en/latest/citations.html for details.

usage: load-graph.py [--version] [--info] [-h] [-k KSIZE] [-U UNIQUE_KMERS] [--fp-rate FP_RATE] [-M MAX_MEMORY_USAGE]
                     [-T THREADS] [--no-build-tagset] [-f]
                     output_nodegraph_filename input_sequence_filename [input_sequence_filename ...]

Load sequences into the compressible graph format plus optional tagset.

positional arguments:
  output_nodegraph_filename
                        output k-mer nodegraph filename.
  input_sequence_filename
                        input FAST[AQ] sequence filename

optional arguments:
  --version             show program's version number and exit
  --info                print citation information
  -h, --help            show this help message and exit
  -k KSIZE, --ksize KSIZE
                        k-mer size to use (default: 32)
  -U UNIQUE_KMERS, --unique-kmers UNIQUE_KMERS
                        approximate number of unique kmers in the input set (default: 0)
  --fp-rate FP_RATE     Override the automatic FP rate setting for the current script (default: None)
  -M MAX_MEMORY_USAGE, --max-memory-usage MAX_MEMORY_USAGE
                        maximum amount of memory to use for data structure (default: None)
  -T THREADS, --threads THREADS
                        Number of simultaneous threads to execute (default: 1)
  --no-build-tagset, -n
                        Do NOT construct tagset while loading sequences (default: False)
  -f, --force           Overwrite output file if it exists (default: False)
```

Example

```bash

# the next command will create a '50m.ct' and a '50m.tagset',
# representing the de Bruijn graph
load-graph.py -k 32 -N 4 -x 16e9 50m iowa-corn-50m.fa.gz
```

### Install KWIP

```bash
# Install kwip
wget https://github.com/kdm9/kWIP/releases/download/0.2.0/kwip-binaries_0.2.0.tar.gz
tar axvf kwip-binaries_0.2.0.tar.gz
mv bin/kwip .
rmdir bin
```

Example

```bash
kwip \
    --threads 4 \                        # Use 4 threads
    --kernel rice.kern \                # Output kernel matrix to ./rice.kern
    --distance rice.dist \                # Output distance matrix to ./rice.dist
    ./hashes/rice_sample_*.ct.gz  # Path to sample hashes, with wildcard
```

## Run

### Run Khmer

```bash
for file in ../data/*.fa.bgz; do
    echo "file ${file}"
    if [[ -f "${file}.khmer" ]]; then
        echo "  exists"
    else
        echo "  counting"
        time load-into-counting.py \
            --ksize 21 \
            -x 5e7 \
            --threads 4 \
            --no-bigcount \
            --summary-info json \
            ${file}.khmer \
            ${file}
        # real    0m37.343s
        # user    2m23.223s
        # sys     0m1.160s
        # 191M ../data/Solanum_tuberosum_PGSC_DM_v3_superscaffolds.fa.bgz.khmer
        # 159  ../data/Solanum_tuberosum_PGSC_DM_v3_superscaffolds.fa.bgz.khmer.info
        # 245  ../data/Solanum_tuberosum_PGSC_DM_v3_superscaffolds.fa.bgz.khmer.info.json
    fi
done
```

### Run Kwip

```bash
./kwip \
    --threads  4 \
    --kernel   test.kern \
    --distance test.dist \
    --verbose \
    ../data/*.khmer

# real     7m 5.453s
# user    18m48.806s
# sys      0m55.190s
```
