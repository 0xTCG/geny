<p align="center" style="font-size:30px;">🩸☠️🦠</p>

# Geny

Geny is a tool for genotyping and calling alleles of the KIR region in the human genome
from whole-genome sequencing data.

## Instalation

Geny needs Python 3.7+ (3.11+ recommended) and [minimap2](https://github.com/lh3/minimap2) to run.

To install Geny, just run:
```bash
pip install -r requirements.txt
```

## Usage

To genotype a sample, just run:
```bash
python geny.py <file.bam> -l <log.txt>
```

The result will be shown on the standard output and the execution log will be stored in `log.txt`. Example:

```
$ python geny.py HG01928.final.cram -l log.txt
Mapping reads: 0...5,000,000...
...
Solution:
- KIR2DL1  0030205     (support=  332)
- KIR2DL1  0030205     (support=  326)
- KIR2DL3  0010101     (support=  404)
- KIR2DL3  0010101     (support=  436)
- KIR2DL4  0010201     (support=  276)
- KIR2DL4  0010201     (support=  325)
- KIR2DP1  0020106     (support=  414)
- KIR2DP1  0020106     (support=  356)
- KIR2DS4  0010102     (support=  464)
- KIR2DS4  0010102     (support=  472)
- KIR3DL1  0150203     (support=  306)
- KIR3DL1  0250101     (support=  307)
- KIR3DL2  0020102     (support=   63)
- KIR3DL2  0020102     (support=   72)
- KIR3DL2  071         (support=   86)
- KIR3DL3  00902       (support=  187)
- KIR3DL3  00902       (support=  127)
- KIR3DP1  0030212     (support=  106)
```

Geny also supports CRAM files. Note that the files should be indexed and aligned to GRCh38/hg38.

If you only have FASTQ data, run Geny as:
```bash
python geny.py <file.fq> -c <coverage>
```

`-c/--coverage` is the expected coverage of the sample **per chromosome**.
For diploid genomes, that is the total sequencing coverage divided by 2.

## Parameters

```
usage: geny.py [-h] [-V] [-l LOG] [-c COVERAGE] [-t THREADS] file

Geny: a tool for genotyping KIR genes

positional arguments:
  file                  input BAM/CRAM/FASTA/FASTQ file

options:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit
  -l LOG, --log LOG     log file location
  -c COVERAGE, --coverage COVERAGE
                        sample coverage
  -t THREADS, --threads THREADS
                        number of threads to use
```

## Paper & Simulations

Geny preprint is available on [bioRxiv]().

Please cite:
> Zhou, Q., Ghezelji, M., et al. Geny: A Genotyping Tool for Allelic Decomposition of Killer Cell Immunoglobulin-Like Receptor Genes. bioRxiv (2024). https://doi.org/10.1186/s13015-022-00210-2

<!-- BibTeX entry:
```
@article{zhou2024geny,
  title={Geny: A Genotyping Tool for Allelic Decomposition of Killer Cell Immunoglobulin-Like Receptor Genes},
  author={I{\v{s}}eri{\'c}, Hamza and Alkan, Can and Hach, Faraz and Numanagi{\'c}, Ibrahim},
  year={2024},
}
``` -->

Experimental data is available in [paper](paper/) directory.

## Changelog

- **Geny v1.0** (Feb 2024)

## Contact

Please reach out to
[Qinghui Zhou](mailto:qinghuiz_at_uvic_dot_ca)
or
[Ibrahim Numanagić](mailto:inumanag_at_uvic_dot_ca).
