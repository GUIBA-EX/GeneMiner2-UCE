# Command Line Interface

For most users of upstream GeneMiner2, downloading prebuilt binaries from Sourceforge is sufficient. The UCE, AliFilter, and combine-stage performance options documented here are fork-specific command-line features; use this repository's scripts or build the CLI from this repository to access them. Refer to [Usage](#usage) to see the documentation.

## Building from Source

Building GeneMiner2 from source is not trivial. Most users of upstream GeneMiner2 can try the prebuilt `.tar.gz` package from Sourceforge. The fork-specific UCE, AliFilter, and combine-stage performance options require running this repository's scripts or building the CLI from this repository. The rest of this section is for advanced users. Please refer to [Usage](#usage) for command-line parameters.

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
conda create -c conda-forge -n geneminer python=3.11 numpy=2.1.3 biopython cython matplotlib pyinstaller scipy
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

Assuming the repository is cloned to `/home/user/GeneMiner2`, the sample list for [DEMO 3](../../DEMO/DEMO3/DEMO3.md) looks as follows:

```
Bupleurum_chinense	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_chinense_1.fq.gz	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_chinense_2.fq.gz
Bupleurum_fruticosum	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_fruticosum_1.fq.gz	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_fruticosum_2.fq.gz
Bupleurum_krylovianum	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_krylovianum_1.fq.gz	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_krylovianum_2.fq.gz
Bupleurum_malconense	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_malconense_1.fq.gz	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_malconense_2.fq.gz
Bupleurum_wenchuanense	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_wenchuanense_1.fq.gz	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_wenchuanense_2.fq.gz
Bupleurum_yunnanense	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_yunnanense_1.fq.gz	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_yunnanense_2.fq.gz
```

The reference sequences have to be under a separate directory. For each gene, place all of its reference sequences in `<Gene Name>.fasta`. For example, to extract matK and psbA genes, create `matK.fasta` and `psbA.fasta` under an empty directory, and write reference sequences into the respective file.

Next, assuming the sample list is saved to `/home/user/GeneMiner2/DEMO/DEMO3/samples.tsv`, Angiosperms353 genes saved under `/home/user/Angiosperms353`, and the desired output location `/home/user/GeneMiner2/DEMO/DEMO3/output`, run GeneMiner2 with default settings:

```
cli/geneminer2 -f /home/user/GeneMiner2/DEMO/DEMO3/samples.tsv -r /home/user/Angiosperms353 -o /home/user/GeneMiner2/DEMO/DEMO3/output
```

GeneMiner2 will build a coalescent tree at `/home/user/GeneMiner2/DEMO/DEMO3/output/Coalescent.tree`.

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
- `--no-trimal`: Deprecated alias for `--alignment-filter none`
- `-cd`: Maximum difference after clean-up
- `-cn`: Minimum number of sequences after clean-up
- `-m`: `coalescent` or `concatenation`
- `-b`: Number of bootstrap replicates
- `--phylo-program`: `fasttree`, `veryfasttree`, `iqtree` or `raxmlng`

For example, after running the command above, you can ask GeneMiner2 to build a concatenation tree based on previous results using `tree` subcommand and `-m concatenation` argument:

```
cli/geneminer2 tree -f /home/user/GeneMiner2/DEMO/DEMO3/samples.tsv -r /home/user/Angiosperm353 -o /home/user/GeneMiner2/DEMO/DEMO3/output -m concatenation
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
cli/geneminer2 trim combine -f /home/user/GeneMiner2/DEMO/DEMO3/samples.tsv -r /mnt/data/Angiosperm353 -o /home/user/GeneMiner2/DEMO/DEMO3/output -ts consensus -tm all -tr 0.5 -cd 0.2 -cn 5 --msa-program muscle
```

All parameters and output are analogus to their counterparts in the graphical version, except that the command line interface only accepts decimal values between 0.0 and 1.0 for percentages. Additionally, several internal options (such as `--min-coverage`) are also exposed, providing extra flexibility for advanced users.
