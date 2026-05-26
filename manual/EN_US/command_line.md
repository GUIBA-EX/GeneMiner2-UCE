# Command Line Interface

This repository contains only the command-line GeneMiner2 workflow. GUI project files, GUI documentation, screenshots, and bundled demo datasets have been removed.

## Building from Source

Building the standalone CLI bundle requires Haxe, a C++ compiler, Python, and PyInstaller. Real analyses should use the bundled `cli/geneminer2` entry point after `make`. Please refer to [Usage](#usage) for command-line parameters.

First, install the default C++ compiler and [zlib](https://zlib.net/) for the platform. On Ubuntu 20.04 and later, this can be done as follows:

```
sudo apt install build-essential zlib1g zlib1g-dev
```

Advanced users can install [zlib-ng](https://github.com/zlib-ng/zlib-ng) or [cloudflare-zlib](https://github.com/cloudflare/zlib) in place of zlib to improve performance.

Second, install [Haxe 4](https://haxe.org/) and [configure the C++ target](https://haxe.org/manual/target-cpp-getting-started.html). On Ubuntu 24.04 and later, the following lines are sufficient:

```
sudo apt install haxe
haxelib setup ~/haxelib
haxelib install hxcpp
```

Third, install Python dependencies. If conda is available, run the following commands. Otherwise, it would be necessary to install them manually using the system package manager.

```
conda create -c conda-forge -n geneminer python=3.11 numpy=2.1.3 biopython cython matplotlib pyinstaller scipy setuptools wheel
conda activate geneminer
```

Finally, download the source code of this fork and build the binaries:

```
git clone --depth 1 https://github.com/GUIBA-EX/GeneMiner2.git
cd GeneMiner2
make
```

The binaries are placed under the `cli` directory. Run `cli/geneminer2 -h` to view the help messages. Several important tools must be available at run time. With conda, these tools can be installed as follows:

```
conda install -c bioconda aster blast clustalo fasttree iqtree mafft magicblast minimap2 muscle raxml-ng trimal veryfasttree
```

AliFilter can be used as an optional alignment filtering program. It is not installed by the command above; download the AliFilter executable separately and make sure `AliFilter` is available in `PATH` before using `--alignment-filter alifilter`.

## Usage

To run an analysis, GeneMiner2 requires a tab-delimited sample list and reference sequences in FASTA format. The sample list has the format `<Species Name><Tab><Read File 1>` (single read) or `<Species Name><Tab><Read File 1><Tab><Read File 2>` (paired-end reads), each line denoting a sample.

For paired-end reads, a sample list looks as follows:

```
Sample_A	/data/reads/Sample_A_R1.fq.gz	/data/reads/Sample_A_R2.fq.gz
Sample_B	/data/reads/Sample_B_R1.fq.gz	/data/reads/Sample_B_R2.fq.gz
Sample_C	/data/reads/Sample_C_R1.fq.gz	/data/reads/Sample_C_R2.fq.gz
```

The reference sequences have to be under a separate directory. For each gene, place all of its reference sequences in `<Gene Name>.fasta`. For example, to extract matK and psbA genes, create `matK.fasta` and `psbA.fasta` under an empty directory, and write reference sequences into the respective file.

Next, assuming the sample list is saved to `/home/user/project/samples.tsv`, reference loci are saved under `/home/user/project/references`, and the desired output location is `/home/user/project/output`, run GeneMiner2 with default settings:

```
cli/geneminer2 -f /home/user/project/samples.tsv -r /home/user/project/references -o /home/user/project/output
```

GeneMiner2 will build a coalescent tree at `/home/user/project/output/Coalescent.tree`.

When `--assembly-mode uce` is used without explicit subcommands, the default workflow skips the reference-based `trim` step and runs `filter refilter assemble combine tree`. This prevents newly recovered UCE flanking regions from being trimmed again during reference-based trimming. Add the `trim` subcommand explicitly if reference-based trimming is still desired.

Command line parameters:

- `-f`: Sample list in tsv format
- `-r`: Reference directory
- `-o`: Output directory
- `-p`: Number of parallel processes
- `--max-reads`: Max million reads per file, disabled by default
- `-kf`: Filter k-mer size
- `-s`: Filter step size
- `-ka`: Assembly k-mer size, automatic by default
- `--min-ka`: Minimum auto-estimated assembly k-mer size
- `--max-ka`: Maximum auto-estimated assembly k-mer size
- `-e`: Error threshold
- `-sb`: Soft boundary (`0`, `auto` or `unlimited`)
- `-i`: Search depth
- `--assembly-mode`: Assembly mode (`reference` or `uce`). `reference` is the default and preserves the existing reference-guided assembly and boundary control; `uce` relaxes reference-boundary trimming and prioritizes longer read-supported UCE flanking sequences.
- In UCE assembly mode, re-filtering keeps paired-end mates together when either mate passes the locus filter. This helps retain flanking reads paired with short UCE probes.
- `--uce-side-candidates`: Number of one-sided branch candidates combined during UCE assembly. Larger values can expose longer low-support flanks but increase runtime and the number of candidate paths.
- `--uce-rescue-reads`: UCE mode only. Run one additional raw-read recruitment round using preliminary contigs plus the original references, then re-filter and re-assemble.
- `--uce-rescue-min-contig-length`: Minimum preliminary contig length used for UCE raw-read rescue.
- `-c`: Consensus threshold (as a decimal, between 0.0 and 1.0)
- `-ts`: Trim source (`assembly` or `consensus`)
- `-tm`: Trim method (`all`, `longest`, `terminal` or `isoform`)
- `-tr`: Retention length threshold (as a decimal, between 0.0 and 1.0)
- `-cs`: Source of combine results (`assembly`, `consensus` or `trimmed`)
- `--msa-program`: `mafft`, `muscle` or `clustalo`
- `--msa-threads`: Threads used by each multiple-sequence-alignment job (default = 1). GeneMiner2 limits concurrent alignment jobs so the total requested alignment threads do not exceed `-p`.
- `--alignment-filter`: Alignment column filtering program before tree reconstruction (`trimal`, `alifilter` or `none`; default = `trimal`)
- `--filter-processes`: Maximum number of concurrent trimAl or AliFilter jobs (default = `-p`)
- `--alifilter-model`: AliFilter model specification or `model.json` path when `--alignment-filter alifilter` is used
- `--strict-combine-errors`: Stop `combine` if any locus fails during multiple-sequence alignment, alignment cleanup, or alignment filtering. By default, failed loci are skipped with a warning.
- `--no-trimal`: Deprecated alias for `--alignment-filter none`
- `-cd`: Maximum difference after clean-up
- `-cn`: Minimum number of sequences after clean-up
- `-m`: `coalescent` or `concatenation`
- `-b`: Number of bootstrap replicates
- `--phylo-program`: `fasttree`, `veryfasttree`, `iqtree` or `raxmlng`

For example, after running the command above, you can ask GeneMiner2 to build a concatenation tree based on previous results using `tree` subcommand and `-m concatenation` argument:

```
cli/geneminer2 tree -f /home/user/project/samples.tsv -r /home/user/project/references -o /home/user/project/output -m concatenation
```

Similarly, it is also possible to ask GeneMiner2 to run specific steps. For example, given these parameters:

| Parameter                  | Value             |
| -------------------------- | ----------------- |
| Source Sequence            | Consensus Results |
| Retention Length Threshold | 50%               |
| Trim Method                | All Fragments     |
| Alignment Program          | muscle            |
| Maximum difference >       | 0.2               |
| Number of sequences        | 5                 |

The following line runs **Trim With Reference** and **Combine Results**:

```
cli/geneminer2 trim combine -f /home/user/project/samples.tsv -r /home/user/project/references -o /home/user/project/output -ts consensus -tm all -tr 0.5 -cd 0.2 -cn 5 --msa-program muscle
```

The command-line interface accepts decimal values between 0.0 and 1.0 for percentage-like thresholds. Several internal options, such as `--min-coverage`, are also exposed for advanced users.
