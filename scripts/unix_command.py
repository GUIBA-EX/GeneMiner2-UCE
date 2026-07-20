from Bio.SeqIO.FastaIO import SimpleFastaParser
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, as_completed, wait
import argparse
import csv
import hashlib
import math
import os
import shutil
import shlex
import subprocess
import sys
import threading


COMMAND_HELP = '''
filter    Reference-based filtering of raw reads
refilter  Refinement of filtered reads
assemble  Gene assembly using wDBG
mito      Mitochondrial GB mining, UCE assembly, overlap merging and circularity QC
profiling Marker read-level profiling (one recruitment, no assembly)
population Build a cohort UCE reference and generate complete, one-per-UCE and LD-pruned SNP panels
gene-annotate Run protein-guided miniprot annotation on recovered gene candidates
gene-resolve Resolve annotated gene candidates into strict SC-OG and unresolved family sets
gene-tree    Infer a species tree from resolved strict or multicopy gene trees
te        Reference-free short-read repeatome discovery, curation, annotation, and quantification
gene      Recover candidate nuclear gene-family copies with original-rust and write conservative cohort summaries
consensus Consensus generation on heterozygous sites
trim      Flank sequence removal
combine   Gene alignment, concatenation and cleanup
tree      Phylogenetic tree reconstruction
stats     UCE recovery statistics and heatmaps
'''

DEPTH_DEPRECATION_EXPLAINER = '''
  Gene assembly involves two different types of depth measurements: k-mer frequencies and read depths. In early versions of GeneMiner2, the two measurements were confounded, causing slight deviation in contig depth calculation. Several new options have since been added to avoid ambiguity.

  --depth-limit INT
    This option limits the highest depth of a gene to prevent pathological assemblies. It should not significantly affect quality unless set to a low value.

  --depth-low-water-mark INT
    This option corresponds to the depth at which GeneMiner2 begins to relax its filtering criteria, so highly diverged sequences can still be recovered. Its value should be sufficiently large, since most irrelevant reads will be removed at a later stage. This option is not for strict quality control!

  --error-threshold INT
    This option exists in previous versions and corresponds to the minimum frequency of a k-mer. Increasing its value will improve quality at the cost of contiguity.

  --min-coverage INT
    This option specifically targets the minimum read depth of contigs. Any contig with a lower read depth will be removed. Note that GeneMiner2 can recover certain types of sequence without a full supporting read. These reads are not counted towards the read depth, making this option very stringent. Use this option if the goal is to ensure read coverage of assembled contigs.
'''

HELP_EPILOG = 'quality control of assembled genes:' + DEPTH_DEPRECATION_EXPLAINER

SCRIPT_ROOT = os.path.join(sys._MEIPASS, os.pardir) if hasattr(sys, '_MEIPASS') else os.path.dirname(__file__)
REFERENCE_EXTENSIONS = ('.fa', '.fas', '.fasta')

def is_reference_file_name(name):
    """瞅瞅这文件名儿是不是咱认的参考序列。"""
    return os.path.splitext(name)[1].lower() in REFERENCE_EXTENSIONS

def materialize_profile_reference(args):
    """Allow profiling to receive one .fa/.fasta file while MainFilter receives a directory."""
    if "profiling" not in getattr(args, "command", ()) or not os.path.isfile(args.r):
        return
    extension = os.path.splitext(args.r)[1].lower()
    if extension not in (".fa", ".fasta"):
        raise RuntimeError("profiling reference must use the .fa or .fasta extension")
    reference_dir = os.path.join(args.o.strip(), ".marker_profile_reference")
    os.makedirs(reference_dir, exist_ok=True)
    link_path = os.path.join(reference_dir, os.path.basename(args.r))
    if os.path.lexists(link_path):
        os.remove(link_path)
    os.symlink(os.path.realpath(args.r), link_path)
    args.r = reference_dir

def prepare_mito_reference(args):
    """Build a per-locus mitochondrial reference directory with the Rust helper."""
    if not getattr(args, "mito_genbank", ""):
        return
    if not os.path.isfile(args.mito_genbank):
        raise RuntimeError(f"Unable to read mitochondrial GenBank file: {args.mito_genbank}")
    if args.mito_flank < 0:
        raise RuntimeError("--mito-flank must be non-negative")
    reference_dir = os.path.join(args.o.strip(), ".gm2_mito_reference")
    mito_bin = find_executable("mito_workflow", internal=True)
    if os.path.isdir(reference_dir):
        shutil.rmtree(reference_dir)
    subprocess.run([mito_bin, "prepare-reference", "--input", args.mito_genbank,
                    "--out-dir", reference_dir, "--flank", str(args.mito_flank),
                    "--tile-length", str(args.mito_tile_length),
                    "--tile-step", str(args.mito_tile_step)], check=True)
    args.r = reference_dir

def run_mito_finalize(args, samples):
    """Finalize the single GM2 UCE assembly with overlaps and paired-read links."""
    if not getattr(args, "mito_genbank", ""):
        return
    mito_bin = find_executable("mito_workflow", internal=True)
    reference_dir = args.r
    reference_genome = os.path.join(reference_dir, "metadata", "mitochondrial_reference.fasta")
    def finalize_sample(sample):
        sample_dir = os.path.join(args.o.strip(), sample)
        command = [
            mito_bin, "finalize",
            "--reference-genome", reference_genome,
            "--contigs", os.path.join(sample_dir, "contigs_all", "mitochondrion.fasta"),
            "--paired-reads", os.path.join(sample_dir, "filtered", "mitochondrion.fq"),
            "--out-dir", os.path.join(sample_dir, "mito"),
            "--minimum-overlap", str(args.mito_min_overlap),
            "--minimum-identity", str(args.mito_min_overlap_identity),
            "--terminal-window", str(args.mito_terminal_window),
            "--link-kmer", str(args.mito_link_kmer),
            "--minimum-link-hits", str(args.mito_min_link_hits),
            "--minimum-pair-support", str(args.mito_min_pair_support),
            "--bridge-kmer", str(args.mito_bridge_kmer),
            "--bridge-minimum-depth", str(args.mito_bridge_min_depth),
            "--maximum-bridge", str(args.mito_max_bridge),
            "--minimum-junction-support", str(getattr(args, "mito_min_junction_support", 3)),
            "--require-circular", "true",
        ]
        subprocess.run(command, check=True)

    workers = min(max(1, getattr(args, 'p', 1)), len(samples))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(finalize_sample, sample): sample for sample in samples}
        for future in as_completed(futures):
            sample = futures[future]
            try:
                future.result()
            except Exception as error:
                raise RuntimeError(
                    f"Mitochondrial finalization failed for {sample}: {error}"
                ) from error
def profile_cache_key(paths, kmer_size):
    digest = hashlib.sha256()
    digest.update(str(kmer_size).encode())
    for path in paths:
        if not path:
            continue
        resolved = os.path.realpath(path)
        digest.update(b'\0')
        digest.update(resolved.encode())
        with open(resolved, 'rb') as handle:
            while chunk := handle.read(1024 * 1024):
                digest.update(chunk)
    return digest.hexdigest()[:16]

def prepare_profile_cache_key(args):
    """Hash immutable profiling inputs once before per-sample work begins."""
    if "profiling" not in getattr(args, "command", ()):
        return
    ref_files = list(iter_reference_files(args.r))
    if len(ref_files) != 1:
        raise RuntimeError("profiling requires exactly one marker reference FASTA")
    bundled_candidates = (
        os.path.normpath(os.path.join(SCRIPT_ROOT, os.pardir, "tools", "themisto-v3.2.2", "themisto_linux-v3.2.2", "themisto")),
        os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(sys.executable)), os.pardir, os.pardir, "tools", "themisto-v3.2.2", "themisto_linux-v3.2.2", "themisto")),
    )
    bundled_themisto = next((path for path in bundled_candidates if os.path.isfile(path)), "")
    args.profile_themisto_bin = args.profile_themisto or bundled_themisto or find_executable("themisto")
    args.profile_reference_path = ref_files[0].path
    args.profile_cache_key = profile_cache_key(
        (args.profile_reference_path, args.profile_group_map, args.profile_decoy, args.profile_themisto_bin),
        args.profile_kmer_size,
    )

def find_executable(prog, internal=False):
    """把要用的程序划拉出来，找不着就麻溜儿报错。"""
    bundled_paths = (
        os.path.join(SCRIPT_ROOT, prog),
        os.path.normpath(os.path.join(SCRIPT_ROOT, os.pardir, "cli", "bin", prog)),
    )
    bin_path = next((path for path in bundled_paths if os.path.isfile(path) and os.access(path, os.X_OK)), None)

    if not bin_path:
        if internal:
            raise RuntimeError(f"A GeneMiner component is missing from {bundled_paths!r}")
        bin_path = shutil.which(prog)

    if not bin_path:
        raise RuntimeError(f"Unable to find {prog} executable")

    return bin_path

def get_ref_genes(ref_dir):
    """把参考目录里的基因名和后缀都归拢出来。"""
    genes = set()

    for entry in iter_reference_files(ref_dir):
        genes.add(os.path.splitext(entry.name))

    return genes

def get_sample_ext(data_path):
    """瞅一眼样本文件，整明白该用 FASTQ 还是 FASTA 后缀。"""
    data_name, data_ext = os.path.splitext(data_path)

    if data_ext == '.gz':
        data_name, data_ext = os.path.splitext(data_name)

    if data_ext == '.fq' or data_ext == '.fastq':
        return '.fq'
    else:
        return '.fasta'

def iter_reference_files(ref_dir):
    """按名儿顺溜地把参考序列文件一个个递出来。"""
    with os.scandir(ref_dir) as entries:
        for entry in sorted(entries, key=lambda x: x.name):
            if not entry.is_file():
                continue

            if is_reference_file_name(entry.name):
                yield entry

def reference_cache_key(ref_dir, kmer_size, step_size):
    """给参考数据和参数摁个指纹，省得缓存整串了。"""
    digest = hashlib.sha256()
    digest.update(os.path.abspath(ref_dir).encode())
    digest.update(b'\0')
    digest.update(str(kmer_size).encode())
    digest.update(b'\0')
    digest.update(str(step_size).encode())

    for entry in iter_reference_files(ref_dir):
        stat = entry.stat()
        digest.update(b'\0')
        digest.update(entry.name.encode())
        digest.update(b'\0')
        digest.update(str(stat.st_size).encode())
        digest.update(b'\0')
        digest.update(str(stat.st_mtime_ns).encode())

    return digest.hexdigest()[:16]

def get_reference_kmer_dict_path(args, out_loc):
    """把参考 k-mer 字典该搁哪儿算明白。"""
    if not args.reuse_reference_cache:
        return os.path.join(out_loc, f'kmer_dict_k{args.kf}.dict')

    cache_dir = args.reference_cache_dir or os.path.join(out_loc, '.gm2_reference_cache')
    cache_name = f'reference_kmer_k{args.kf}_s{args.step_size}_{reference_cache_key(args.r, args.kf, args.step_size)}.dict'
    return os.path.join(cache_dir, cache_name)

def get_assembler_reference_cache_dir(args, out_loc):
    """瞅瞅要不要复用组装参考缓存，再把地方定下来。"""
    if not args.reuse_reference_cache:
        return None

    return os.path.join(args.reference_cache_dir or os.path.join(out_loc, '.gm2_reference_cache'), 'assembler')

def prepare_workdir(args):
    """读样本表、收拾样本名，再把干活目录都支棱起来。"""
    samples = {}
    tsv_loc = args.f

    try:
        sp_id = 0
        tsv_loc = os.path.realpath(tsv_loc, strict=True)

        with open(tsv_loc, 'r') as f:
            for row in csv.reader(f, delimiter="\t", quotechar='"'):
                if not row:
                    continue

                sp_id += 1
                sp_name = "".join('_' if c in ' -' else c for c in row[0].strip() if c.isalnum() or c in " -_.").capitalize()

                if not sp_name:
                    print(f"Invalid sample name '{row[0]}'")
                    return {}

                if len(row) == 1:
                    print(f"Sample '{row[0]}' has no data files")
                    return {}

                samples[f'{sp_id}_{sp_name}'] = (row[1], row[1] if len(row) == 2 else row[2])

    except OSError as e:
        print(f"Unable to read sample list '{tsv_loc}': {e}")
        return {}

    out_loc = args.o.strip()

    print(f"Preparing working directory '{out_loc}'")

    try:
        os.makedirs(out_loc, exist_ok=True)
    except OSError as e:
        print(f"Unable to create working directory '{out_loc}': {e}")
        return {}

    for name in samples.keys():
        sp_path = os.path.join(out_loc, name)

        try:
            os.makedirs(sp_path, exist_ok=True)
        except OSError as e:
            print(f"Unable to create directory '{sp_path}': {e}")
            return {}

    return samples

def write_fasta_record(out, header, sequence, line_width=80):
    """把一条序列规规矩矩写成 FASTA，空的咱可不写。"""
    sequence = ''.join(sequence.split()).upper()

    if not sequence:
        return False

    out.write(f'>{header}\n')

    for i in range(0, len(sequence), line_width):
        out.write(sequence[i:i + line_width] + '\n')

    return True

def build_uce_rescue_refs(ref_dir, sample_dir, rescue_ref_dir, min_contig_len):
    """拿原参考和靠谱 contig 拼一套 UCE 救援参考。"""
    results_dir = os.path.join(sample_dir, 'results')
    summary_rows = read_uce_summary(os.path.join(sample_dir, 'uce_assembly_summary.csv'))
    added_contigs = 0

    if not os.path.isdir(results_dir):
        return 0

    if os.path.isdir(rescue_ref_dir):
        shutil.rmtree(rescue_ref_dir, ignore_errors=True)

    os.makedirs(rescue_ref_dir, exist_ok=True)

    with os.scandir(ref_dir) as entries:
        for entry in entries:
            if not entry.is_file():
                continue

            if not is_reference_file_name(entry.name):
                continue

            gene = os.path.splitext(entry.name)[0]
            contig_path = os.path.join(results_dir, f'{gene}.fasta')
            rescue_path = os.path.join(rescue_ref_dir, f'{gene}.fasta')
            contig_index = 0
            wrote_any = False

            with open(rescue_path, 'w') as out:
                with open(entry.path, 'r') as ref_in:
                    for title, seq in SimpleFastaParser(ref_in):
                        wrote_any |= write_fasta_record(out, title, seq)

                if (os.path.isfile(contig_path)
                        and uce_summary_row_is_accepted(summary_rows.get(gene))):
                    with open(contig_path, 'r') as contig_in:
                        for _, seq in SimpleFastaParser(contig_in):
                            if len(seq) < min_contig_len:
                                continue

                            contig_index += 1
                            added_contigs += 1
                            wrote_any |= write_fasta_record(out, f'{gene}_gm2_rescue_contig_{contig_index}', seq)

            if not wrote_any:
                os.remove(rescue_path)

    return added_contigs

def read_uce_summary(summary_path):
    """把 UCE 汇总表按 locus 收拢成一摞，后头好查。"""
    rows = {}

    if not os.path.isfile(summary_path):
        return rows

    with open(summary_path, newline='') as f:
        for row in csv.DictReader(f):
            locus = row.get('locus')

            if locus:
                rows[locus] = row

    return rows


def uce_summary_row_is_accepted(row):
    """瞅瞅这条 locus 收没收，老版汇总表也照样认。"""
    if not row:
        return False

    accepted = str(row.get('accepted', '')).strip().lower()
    if accepted:
        return accepted in {'1', 'true', 'yes'}

    low_quality = str(row.get('low_quality', '')).strip().lower()
    return row.get('status') == 'success' and low_quality not in {'1', 'true', 'yes'}

def int_or_blank(value):
    """能整成整数就整，整不了就撂空白。"""
    try:
        return int(value)
    except (TypeError, ValueError):
        return ''


def float_or_blank(value):
    """能整成小数就整，整不了就撂空白。"""
    try:
        return float(value)
    except (TypeError, ValueError):
        return ''

def delta_or_blank(after, before):
    """前后数值能对上就算差，对不上咱就留空。"""
    after_value = int_or_blank(after)
    before_value = int_or_blank(before)

    if after_value == '' or before_value == '':
        return ''

    return after_value - before_value

def read_density_or_blank(row):
    """从汇总行里把读段密度掰扯明白，没数就留空。"""
    length = int_or_blank(row.get('selected_contig_length'))
    read_count = int_or_blank(row.get('unique_read_count'))
    if length != '' and length > 0 and read_count != '':
        return read_count / length

    unique_density = float_or_blank(row.get('unique_read_density'))
    if unique_density != '':
        return unique_density

    read_count = int_or_blank(row.get('read_count'))
    if length == '' or read_count == '' or length <= 0:
        return ''

    return read_count / length

def density_ratio_or_blank(before, after):
    """算救援前后的密度倍数，底数不靠谱就甭硬算。"""
    before_density = read_density_or_blank(before)
    after_density = read_density_or_blank(after)

    if before_density == '' or after_density == '':
        return ''

    if before_density <= 0:
        return ''

    return after_density / before_density

def rescue_density_below_ratio(before, after, min_density_ratio):
    """瞅瞅救援后的密度掉没掉过警戒线。"""
    density_ratio = density_ratio_or_blank(before, after)

    if density_ratio == '':
        return False

    return density_ratio < min_density_ratio

UCE_ASSEMBLY_SUMMARY_FIELDS = [
    'locus',
    'status',
    'accepted',
    'rejection_reason',
    'selected_contig_length',
    'read_supported_span',
    'slice_supported_bases',
    'slice_support_breadth',
    'max_slice_support_gap',
    'read_count',
    'unique_read_count',
    'multi_mapping_read_count',
    'read_density',
    'unique_read_density',
    'support_fraction',
    'flank_balance',
    'kmer_median_depth',
    'kmer_depth_cv',
    'kmer_max_depth_ratio',
    'candidate_count',
    'low_quality',
]

UCE_RESCUE_SUMMARY_FIELDS = [
    'sample',
    'locus',
    'rescue_status',
    'before_status',
    'after_status',
    'before_length',
    'after_length',
    'length_delta',
    'before_read_count',
    'after_read_count',
    'read_count_delta',
    'before_read_density',
    'after_read_density',
    'density_ratio',
    'before_read_supported_span',
    'after_read_supported_span',
    'error',
]

SAMPLE_STATE_BACKUP_ITEMS = [
    'results',
    'result_dict.txt',
    'uce_assembly_summary.csv',
    'contigs_all',
    'contigs_all_low',
    'filtered',
    'filtered_pe',
    'ref_reads_count_dict.txt',
]

def write_uce_assembly_summary_rows(summary_path, rows):
    """把 UCE 组装汇总按固定列、固定顺序写利索。"""
    with open(summary_path, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=UCE_ASSEMBLY_SUMMARY_FIELDS)
        writer.writeheader()

        for locus in sorted(rows):
            row = rows[locus]
            writer.writerow({field: row.get(field, '') for field in UCE_ASSEMBLY_SUMMARY_FIELDS})

def write_result_dict_from_uce_summary(sample_dir, rows):
    """照 UCE 汇总重整 result_dict，跳过的 locus 不往里塞。"""
    result_path = os.path.join(sample_dir, 'result_dict.txt')

    with open(result_path, 'w') as out:
        for locus in sorted(rows):
            row = rows[locus]

            if row.get('status') == 'skipped':
                continue

            out.write(f"{locus},{row.get('status', '')},{row.get('read_count', '')},\n")

def restore_locus_file(sample_dir, backup_dir, subdir, locus):
    """把单个 locus 的结果文件从备份里原样倒腾回来。"""
    rel_path = os.path.join(subdir, f'{locus}.fasta')
    src = os.path.join(backup_dir, rel_path)
    dest = os.path.join(sample_dir, rel_path)

    if os.path.isfile(src):
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copy2(src, dest)
    elif os.path.isfile(dest):
        os.remove(dest)

def locus_file_name_matches(name, locus, paired=False):
    """瞅瞅文件名跟这个 locus 对不对号，双端名儿也管。"""
    stem = os.path.splitext(name)[0]
    if paired:
        return stem in (f'{locus}_1', f'{locus}_2')
    return stem == locus

def restore_locus_directory_files(sample_dir, backup_dir, subdir, locus):
    """只还原这个 locus 的读段文件，救成的那些咱不碰。"""
    source_dir = os.path.join(backup_dir, subdir)
    destination_dir = os.path.join(sample_dir, subdir)
    names = set()

    for directory in (source_dir, destination_dir):
        if os.path.isdir(directory):
            names.update(entry.name for entry in os.scandir(directory) if entry.is_file())

    for name in names:
        if not locus_file_name_matches(name, locus, paired=subdir == 'filtered_pe'):
            continue

        source = os.path.join(source_dir, name)
        destination = os.path.join(destination_dir, name)
        if os.path.isfile(source):
            os.makedirs(destination_dir, exist_ok=True)
            shutil.copy2(source, destination)
        elif os.path.isfile(destination):
            os.remove(destination)

def restore_locus_read_count(sample_dir, backup_dir, locus):
    """把这个 locus 原先的读段计数塞回当前计数表。"""
    filename = 'ref_reads_count_dict.txt'
    source = os.path.join(backup_dir, filename)
    destination = os.path.join(sample_dir, filename)

    def read_rows(path):
        """把非空计数行都划拉出来，文件没有就算了。"""
        if not os.path.isfile(path):
            return []
        with open(path) as handle:
            return [line for line in handle if line.strip()]

    backup_rows = read_rows(source)
    current_rows = read_rows(destination)
    backup_locus_rows = [line for line in backup_rows if line.split(',', 1)[0] == locus]
    merged_rows = [line for line in current_rows if line.split(',', 1)[0] != locus]
    merged_rows.extend(backup_locus_rows)

    if merged_rows:
        with open(destination, 'w') as handle:
            handle.writelines(merged_rows)
    elif os.path.isfile(destination):
        os.remove(destination)

def format_float_or_blank(value, digits=6):
    """小数收拾利索再输出，空值还让它空着。"""
    if value == '':
        return ''

    return f'{value:.{digits}f}'.rstrip('0').rstrip('.')

def revert_invalid_rescue_loci(sample_dir, backup_dir, before_rows, rescue_rows, min_density_ratio):
    """救援后变孬的 locus 给它退回原样，别硬留着。"""
    reverted = {}
    final_rows = {locus: row.copy() for locus, row in rescue_rows.items()}

    for locus, before in before_rows.items():
        after = rescue_rows.get(locus)
        if not uce_summary_row_is_accepted(before):
            continue

        if not uce_summary_row_is_accepted(after):
            status = 'reverted_failed_rescue'
        elif rescue_density_below_ratio(before, after, min_density_ratio):
            status = 'reverted_density_drop'
        else:
            continue

        for subdir in ('results', 'contigs_all', 'contigs_all_low'):
            restore_locus_file(sample_dir, backup_dir, subdir, locus)
        for subdir in ('filtered', 'filtered_pe'):
            restore_locus_directory_files(sample_dir, backup_dir, subdir, locus)
        restore_locus_read_count(sample_dir, backup_dir, locus)

        final_rows[locus] = before.copy()
        reverted[locus] = status

    if reverted:
        write_uce_assembly_summary_rows(os.path.join(sample_dir, 'uce_assembly_summary.csv'), final_rows)
        write_result_dict_from_uce_summary(sample_dir, final_rows)

    return reverted

def write_sample_uce_rescue_summary(sample_dir, sample, before_rows, after_rows, rescue_status, error='', status_by_locus=None, error_by_locus=None):
    """把一个样本救援前后的变化明明白白写进表里。"""
    out_path = os.path.join(sample_dir, 'uce_rescue_summary.csv')
    loci = sorted(set(before_rows) | set(after_rows))
    status_by_locus = {} if status_by_locus is None else status_by_locus
    error_by_locus = {} if error_by_locus is None else error_by_locus

    with open(out_path, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=UCE_RESCUE_SUMMARY_FIELDS)
        writer.writeheader()

        for locus in loci:
            before = before_rows.get(locus, {})
            after = after_rows.get(locus, {})
            before_length = before.get('selected_contig_length', '')
            after_length = after.get('selected_contig_length', '')
            before_count = before.get('read_count', '')
            after_count = after.get('read_count', '')
            before_density = read_density_or_blank(before)
            after_density = read_density_or_blank(after)

            writer.writerow({
                'sample': sample,
                'locus': locus,
                'rescue_status': status_by_locus.get(locus, rescue_status),
                'before_status': before.get('status', ''),
                'after_status': after.get('status', ''),
                'before_length': before_length,
                'after_length': after_length,
                'length_delta': delta_or_blank(after_length, before_length),
                'before_read_count': before_count,
                'after_read_count': after_count,
                'read_count_delta': delta_or_blank(after_count, before_count),
                'before_read_density': format_float_or_blank(before_density),
                'after_read_density': format_float_or_blank(after_density),
                'density_ratio': format_float_or_blank(density_ratio_or_blank(before, after)),
                'before_read_supported_span': before.get('read_supported_span', ''),
                'after_read_supported_span': after.get('read_supported_span', ''),
                'error': error_by_locus.get(locus, error),
            })

def backup_sample_state(sample_dir):
    """救援前先把样本现场挪走备份，留条后路。"""
    backup_dir = os.path.join(sample_dir, '.uce_rescue_backup')

    if os.path.isdir(backup_dir):
        shutil.rmtree(backup_dir, ignore_errors=True)

    os.makedirs(backup_dir, exist_ok=True)

    for item in SAMPLE_STATE_BACKUP_ITEMS:
        src = os.path.join(sample_dir, item)

        if os.path.exists(src):
            shutil.move(src, os.path.join(backup_dir, item))

    return backup_dir

def restore_sample_state(sample_dir, backup_dir):
    """救援整岔劈了，就把样本现场从备份还原回来。"""
    if not os.path.isdir(backup_dir):
        return

    for item in SAMPLE_STATE_BACKUP_ITEMS:
        dest = os.path.join(sample_dir, item)

        if os.path.isdir(dest):
            shutil.rmtree(dest, ignore_errors=True)
        elif os.path.isfile(dest):
            os.remove(dest)

    with os.scandir(backup_dir) as entries:
        for entry in entries:
            dest = os.path.join(sample_dir, entry.name)
            shutil.move(entry.path, dest)

    shutil.rmtree(backup_dir, ignore_errors=True)

def discard_sample_state_backup(backup_dir):
    """救援稳当了就把临时备份收拾掉，别占地方。"""
    if os.path.isdir(backup_dir):
        shutil.rmtree(backup_dir, ignore_errors=True)

def write_failed_samples(out_loc, failures):
    """把失败样本单列出来；一个没有就把旧表撤了。"""
    out_path = os.path.join(out_loc, 'failed_samples.tsv')

    if not failures:
        if os.path.isfile(out_path):
            os.remove(out_path)

        return

    with open(out_path, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['sample', 'stage', 'error'])
        writer.writerows(failures)

def get_rescue_sample_names(samples, failures):
    """前头没整成的样本挑出去，剩下的再进救援。"""
    # 前头都失败了就别硬救了，省得越整越乱。
    failed = {sample for sample, _, _ in failures}
    return [name for name in samples if name not in failed]

def get_uce_rescue_parallelism(total_threads, sample_count):
    """按总线程和样本数掂量救援咋分工最合适。"""
    rescue_threads = max(1, min(4, total_threads))
    rescue_workers = max(1, min(4, sample_count, total_threads // rescue_threads))
    return rescue_workers, rescue_threads

def build_uce_rescue_filter_commands(filter_bin, rescue_ref_dir, sample_dir, q1, q2, args, rescue_kmer_dict_path):
    """把 UCE 救援过滤要跑的两趟命令整齐备好。"""
    dict_cmd = [filter_bin, '-r', rescue_ref_dir, '-o', sample_dir, '-kf', str(args.kf),
                '-s', str(args.step_size), '-gr', '-lkd', rescue_kmer_dict_path, '-m', '2']
    reads_cmd = [filter_bin, '-r', rescue_ref_dir, '-q1', q1, '-q2', q2, '-o', sample_dir,
                 '-kf', str(args.kf), '-s', str(args.step_size), '-gr', '-subdir', 'filtered_pe',
                 '-m', '5', '-lb', '-lkd', rescue_kmer_dict_path]

    if args.max_reads > 0:
        reads_cmd.extend(['-m_reads', str(args.max_reads)])

    return dict_cmd, reads_cmd

def build_assembler_command(assembler_bin, args, sample_dir, ref_dir, soft_boundary, thr, backend='uce-rust'):
    """照实现类型拼组装命令，老原版不掺 UCE 新参数。"""
    command = [
        assembler_bin, '-r', ref_dir, '-o', sample_dir, '-ka', str(args.ka),
        '-k_min', str(args.min_ka), '-k_max', str(args.max_ka),
        '-limit_count', str(args.error_threshold), '-iteration', str(args.search_depth),
        '-sb', soft_boundary, '-cov_min', str(args.min_coverage), '-p', str(thr),
    ]

    # 原版算法后端只认老参数；Rust 复刻版额外认 reference cache，UCE 参数一概不塞。
    if backend in ('original', 'original-rust'):
        if backend == 'original-rust':
            assembler_cache_dir = getattr(args, 'assembler_reference_cache_dir', None)
            original_ref_dir = getattr(args, 'r', None)
            if assembler_cache_dir and original_ref_dir and os.path.abspath(ref_dir) != os.path.abspath(original_ref_dir):
                assembler_cache_dir = None
            if assembler_cache_dir:
                command.extend(['--assembler-reference-cache-dir', assembler_cache_dir])
        return command

    command.extend([
        '--assembly-mode', args.assembly_mode,
        '--uce-side-candidates', str(args.uce_side_candidates),
        '--uce-max-contig-length', str(args.uce_max_contig_length),
        '--uce-min-read-density', str(args.uce_min_read_density),
        '--uce-density-check-min-length', str(args.uce_density_check_min_length),
        '--uce-max-depth-cv', str(args.uce_max_depth_cv),
        '--uce-max-depth-ratio', str(args.uce_max_depth_ratio),
    ])

    if backend == 'uce-rust':
        command.extend([
            '--uce-path-strategy', getattr(args, 'uce_path_strategy', 'backbone'),
            '--uce-backbone-lookahead', str(getattr(args, 'uce_backbone_lookahead', 24)),
            '--assembler-read-chunk-size', str(getattr(args, 'assembler_read_chunk_size', 8192)),
            '--assembler-kmer-count-threads', str(getattr(args, 'assembler_kmer_count_threads', 0)),
            '--assembler-graph-format', getattr(args, 'assembler_graph_format', 'none'),
        ])

    assembler_cache_dir = getattr(args, 'assembler_reference_cache_dir', None)
    original_ref_dir = getattr(args, 'r', None)
    if assembler_cache_dir and original_ref_dir and os.path.abspath(ref_dir) != os.path.abspath(original_ref_dir):
        assembler_cache_dir = None

    if assembler_cache_dir:
        command.extend(['--assembler-reference-cache-dir', assembler_cache_dir])

    return command

def do_filter_assemble(args, samples, do_filter, do_refilter, do_assemble, ignore_hook=lambda *_, **__: None):
    """把过滤、再过滤、组装和救援这一大趟活儿串起来。"""
    out_loc = args.o.strip()
    is_profiling = "profiling" in getattr(args, "command", ())
    is_mito = bool(getattr(args, "is_mito_workflow", False))
    kmer_dict_path = get_reference_kmer_dict_path(args, out_loc)
    args.assembler_reference_cache_dir = get_assembler_reference_cache_dir(args, out_loc)
    rescue_enabled = args.uce_rescue_reads
    failed_samples = []
    rescue_workers, rescue_threads = get_uce_rescue_parallelism(args.p, len(samples))

    if rescue_enabled and (args.assembly_mode != 'uce' or not do_filter or not do_refilter or not do_assemble):
        raise RuntimeError('--uce-rescue-reads requires --assembly-mode uce and the filter, refilter and assemble steps')

    if args.soft_boundary == 'auto':
        soft_boundary = -1
    elif args.soft_boundary == 'unlimited':
        soft_boundary = 10000
    else:
        try:
            soft_boundary = int(args.soft_boundary)
        except ValueError:
            raise RuntimeError(f"Invalid soft boundary {args.soft_boundary} (must be an integer)")

        if soft_boundary < 0:
            raise RuntimeError(f"Invalid soft boundary {args.soft_boundary} (must be positive or zero)")

    soft_boundary = str(soft_boundary)

    if do_filter:
        filter_bin = find_executable('MainFilterNew', internal=True)

        os.makedirs(os.path.dirname(kmer_dict_path), exist_ok=True)

        if os.path.isfile(kmer_dict_path) and args.reuse_reference_cache:
            print(f'Reusing reference k-mer cache: {kmer_dict_path}')
        elif os.path.isfile(kmer_dict_path):
            os.remove(kmer_dict_path)

        if not os.path.isfile(kmer_dict_path):
            try:
                subprocess.run([filter_bin, '-r', args.r, '-o', out_loc, '-kf', str(args.kf), '-s', str(args.step_size),
                                '-gr', '-lkd', kmer_dict_path, '-m', '2'], check=True)
            except subprocess.SubprocessError as e:
                raise RuntimeError(f"Unable to build k-mer dictionary: {e}")

        def run_filter(name):
            """给这个样本捞参考相关读段，再把输出归拢好。"""
            q1, q2 = samples[name]
            read_count_path = os.path.join(out_loc, name, 'ref_reads_count_dict.txt')
            out_dir = os.path.join(out_loc, name, 'filtered_pe')

            if os.path.isfile(read_count_path):
                os.remove(read_count_path)

            if os.path.isdir(out_dir):
                shutil.rmtree(out_dir, ignore_errors=True)

            params = [filter_bin, '-r', args.r, '-q1', q1, '-q2', q2, '-o', os.path.join(out_loc, name),
                      '-kf', str(args.kf), '-s', str(args.step_size), '-gr', '-subdir', 'filtered_pe',
                      '-m', '4' if is_profiling or is_mito else '5', '-lb', '-lkd', kmer_dict_path]

            if args.max_reads > 0:
                params.extend(['-m_reads', str(args.max_reads)])

            subprocess.run(params, check=True)

            if not os.path.isfile(read_count_path):
                raise RuntimeError('Filter failed')

            if is_mito:
                collapsed_dir = out_dir + '_collapsed'
                if os.path.isdir(collapsed_dir):
                    shutil.rmtree(collapsed_dir, ignore_errors=True)
                subprocess.run([find_executable('mito_workflow', internal=True), 'collapse-baits',
                                '--input-dir', out_dir, '--out-dir', collapsed_dir,
                                '--output-name', 'mitochondrion'], check=True)
                shutil.rmtree(out_dir, ignore_errors=True)
                os.replace(collapsed_dir, out_dir)

            if not do_refilter and os.path.isdir(out_dir):
                merge_dir = os.path.join(out_loc, name, 'filtered')
                sample_ext = get_sample_ext(q1)

                if os.path.isdir(merge_dir):
                    shutil.rmtree(merge_dir, ignore_errors=True)

                os.makedirs(merge_dir, exist_ok=True)

                genes = set()

                with open(read_count_path, 'r') as f:
                    for line in f:
                        line = line.strip()

                        if not line:
                            continue

                        genes.add(line.split(',')[0])

                for gene in genes:
                    read_1 = os.path.join(out_dir, f'{gene}_1{sample_ext}')
                    read_2 = os.path.join(out_dir, f'{gene}_2{sample_ext}')

                    if not os.path.isfile(read_1):
                        continue

                    with open(os.path.join(merge_dir, gene + sample_ext), 'wb') as f:
                        with open(read_1, 'rb') as r:
                            shutil.copyfileobj(r, f)

                        if not os.path.isfile(read_2):
                            continue

                        with open(read_2, 'rb') as r:
                            shutil.copyfileobj(r, f)

    else:
        run_filter = ignore_hook

    if do_refilter:
        refilter_bin = find_executable('main_refilter_new', internal=True)

        def run_refilter(name, thr=1, ref_dir=None):
            """把样本读段再筛一遍，杂的赖的往外挑。"""
            in_dir  = os.path.join(out_loc, name, 'filtered_pe')
            out_dir = os.path.join(out_loc, name, 'filtered')
            ref_dir = args.r if ref_dir is None else ref_dir

            if not os.path.isdir(in_dir):
                raise RuntimeError('No successful filter run, cannot re-filter')

            if os.path.isdir(out_dir):
                shutil.rmtree(out_dir, ignore_errors=True)

            params = [refilter_bin, '-r', ref_dir, '-qd', in_dir, '-o', out_dir, '-kf', str(args.kf),
                      '-p', str(thr), '--log-file', os.path.join(out_loc, name, 'log.txt'),
                      '--min-depth', str(args.depth_low_water_mark), '--max-depth', str(args.depth_limit),
                      '--max-size', str(args.file_size_limit)]
            if not is_mito:
                params.append('--use-gm2-format')

            if args.assembly_mode == 'uce' or is_profiling:
                params.append('--keep-linked-mates')

            subprocess.run(params, check=True)

            if do_filter and os.path.isdir(in_dir) and os.path.isdir(out_dir):
                shutil.rmtree(in_dir, ignore_errors=True)

    else:
        run_refilter = ignore_hook

    if do_assemble:
        assembler_implementation = getattr(args, 'assembler_implementation', 'auto')
        original_assembler_bin = None
        original_rust_assembler_bin = None
        rust_assembler_bin = None

        # original 模式默认用 original-rust；UCE 默认用 uce-rust，original 仍保留给上游 Python 对照。
        if args.assembly_mode == 'original':
            if assembler_implementation in ('auto', 'original-rust'):
                original_rust_assembler_bin = find_executable('main_assembler-original-rust', internal=True)
            elif assembler_implementation != 'uce-rust':
                original_assembler_bin = find_executable('main_assembler-original', internal=True)
            else:
                rust_assembler_bin = find_executable('main_assembler-rust', internal=True)
        elif args.assembly_mode == 'uce':
            if assembler_implementation in ('original', 'original-rust'):
                raise RuntimeError(
                    f'{args.assembly_mode.upper()} mode requires the Rust UCE assembler'
                )
            rust_assembler_bin = find_executable('main_assembler-rust', internal=True)

        def run_assembler(name, thr=1, ref_dir=None):
            """组装这个样本；original 默认 original-rust，UCE 默认 uce-rust。"""
            sample_dir = os.path.join(out_loc, name)
            in_dir = os.path.join(sample_dir, 'filtered')
            out_dir = os.path.join(sample_dir, 'results')
            result_path = os.path.join(sample_dir, 'result_dict.txt')
            uce_summary_path = os.path.join(sample_dir, 'uce_assembly_summary.csv')
            ref_dir = args.r if ref_dir is None else ref_dir

            if not os.path.isdir(in_dir):
                raise RuntimeError('No successful filter run, cannot assemble')

            if is_profiling:
                reads = [
                    entry.path for entry in sorted(os.scandir(in_dir), key=lambda entry: entry.name)
                    if entry.is_file() and get_sample_ext(entry.name) in ('.fq', '.fasta')
                ]
                if not reads:
                    raise RuntimeError('No recruited marker reads found for profiling')
                profile_dir = os.path.join(sample_dir, 'marker_profile')
                if os.path.isdir(profile_dir):
                    shutil.rmtree(profile_dir, ignore_errors=True)
                if len(reads) != 1:
                    raise RuntimeError('profiling requires exactly one merged recruited-read file')
                quant_bin = find_executable('marker_profile', internal=True)
                cache_root = args.profile_index_dir or args.reference_cache_dir or os.path.join(out_loc, '.gm2_reference_cache')
                cache_dir = os.path.join(cache_root, f'profile_themisto_k{args.profile_kmer_size}_{args.profile_cache_key}')
                command = [
                    quant_bin, '--reference', args.profile_reference_path, '--reads', reads[0],
                    '--output', profile_dir, '--cache', cache_dir,
                    '--themisto', args.profile_themisto_bin,
                    '--threads', str(thr), '--kmer-size', str(args.profile_kmer_size),
                    '--threshold', str(args.profile_pseudoalign_threshold),
                    '--relevant-kmer-fraction', str(args.profile_relevant_kmer_fraction),
                    '--index-memory-gb', str(args.profile_index_memory_gb),
                ]
                if args.profile_group_map:
                    command.extend(['--groups', args.profile_group_map])
                if args.profile_decoy:
                    command.extend(['--decoy', args.profile_decoy])
                if args.profile_force_rebuild:
                    command.append('--force-rebuild')
                subprocess.run(command, check=True)
                if not os.path.isfile(os.path.join(profile_dir, 'marker_reference_support.tsv')):
                    raise RuntimeError('profiling failed to produce marker_reference_support.tsv')
                return

            def clear_assembly_outputs():
                """开整前把旧组装产物清出去，省得串锅。"""
                if os.path.isdir(out_dir):
                    shutil.rmtree(out_dir, ignore_errors=True)
                for candidate_dir in ('contigs_all', 'contigs_all_low'):
                    path = os.path.join(sample_dir, candidate_dir)
                    if os.path.isdir(path):
                        shutil.rmtree(path, ignore_errors=True)
                graph_dir = os.path.join(sample_dir, 'assembly_graphs')
                if os.path.isdir(graph_dir):
                    shutil.rmtree(graph_dir, ignore_errors=True)
                for path in (result_path, uce_summary_path, os.path.join(sample_dir, 'marker_profile_summary.csv')):
                    if os.path.isfile(path):
                        os.remove(path)

            def execute_assembler(executable, backend='uce-rust'):
                """真把组装器跑起来，再瞅瞅结果落地没。"""
                clear_assembly_outputs()
                command = build_assembler_command(
                    executable, args, sample_dir, ref_dir, soft_boundary, thr, backend=backend
                )
                subprocess.run(command, check=True)
                if not os.path.isfile(result_path):
                    raise RuntimeError('Assembly failed to produce result_dict.txt')

            if original_assembler_bin is not None:
                execute_assembler(original_assembler_bin, backend='original')
            elif original_rust_assembler_bin is not None:
                execute_assembler(original_rust_assembler_bin, backend='original-rust')
            else:
                execute_assembler(rust_assembler_bin, backend='uce-rust')

    else:
        run_assembler = ignore_hook

    if rescue_enabled:
        def run_uce_rescue(name, thr=1):
            """拿初组装结果再捞一轮读段，救救这个 UCE 样本。"""
            sample_dir = os.path.join(out_loc, name)
            rescue_ref_dir = os.path.join(sample_dir, 'uce_rescue_refs')
            rescue_kmer_dict_path = os.path.join(sample_dir, f'uce_rescue_kmer_dict_k{args.kf}.dict')
            summary_path = os.path.join(sample_dir, 'uce_assembly_summary.csv')
            read_count_path = os.path.join(sample_dir, 'ref_reads_count_dict.txt')
            filtered_pe_dir = os.path.join(sample_dir, 'filtered_pe')
            filtered_dir = os.path.join(sample_dir, 'filtered')

            before_rows = read_uce_summary(summary_path)
            added_contigs = build_uce_rescue_refs(args.r, sample_dir, rescue_ref_dir, args.uce_rescue_min_contig_length)

            if added_contigs == 0:
                print(f'No preliminary UCE contigs for {name}; skipping raw-read rescue.')
                write_sample_uce_rescue_summary(sample_dir, name, before_rows, before_rows, 'skipped')
                return

            print(f'Running one-round UCE raw-read rescue for {name} using {added_contigs} preliminary contig(s).')
            backup_dir = backup_sample_state(sample_dir)

            try:
                if os.path.isfile(rescue_kmer_dict_path):
                    os.remove(rescue_kmer_dict_path)

                q1, q2 = samples[name]
                dict_cmd, reads_cmd = build_uce_rescue_filter_commands(
                    filter_bin,
                    rescue_ref_dir,
                    sample_dir,
                    q1,
                    q2,
                    args,
                    rescue_kmer_dict_path,
                )

                subprocess.run(dict_cmd, check=True)

                if os.path.isfile(read_count_path):
                    os.remove(read_count_path)

                if os.path.isdir(filtered_pe_dir):
                    shutil.rmtree(filtered_pe_dir, ignore_errors=True)

                if os.path.isdir(filtered_dir):
                    shutil.rmtree(filtered_dir, ignore_errors=True)

                subprocess.run(reads_cmd, check=True)

                if not os.path.isfile(read_count_path):
                    raise RuntimeError('UCE rescue filter failed')

                run_refilter(name, thr=thr, ref_dir=rescue_ref_dir)
                run_assembler(name, thr=thr, ref_dir=rescue_ref_dir)

            except Exception as e:
                restore_sample_state(sample_dir, backup_dir)
                write_sample_uce_rescue_summary(sample_dir, name, before_rows, before_rows, 'failed_rolled_back', str(e))
                print(f'Warning: UCE raw-read rescue failed for {name}; first-round assembly was restored: {e}')
                return

            else:
                after_rows = read_uce_summary(summary_path)
                reverted_loci = revert_invalid_rescue_loci(
                    sample_dir,
                    backup_dir,
                    before_rows,
                    after_rows,
                    args.uce_rescue_min_density_ratio,
                )
                status_by_locus = reverted_loci
                error_by_locus = {}
                for locus, status in reverted_loci.items():
                    if status == 'reverted_density_drop':
                        error_by_locus[locus] = (
                            f'rescue unique read density ratio below '
                            f'{args.uce_rescue_min_density_ratio:g}; first-round contig restored')
                    else:
                        error_by_locus[locus] = (
                            'rescue result missing or rejected; first-round contig restored')
                final_rows = read_uce_summary(summary_path)
                write_sample_uce_rescue_summary(sample_dir, name, before_rows, final_rows, 'success', status_by_locus=status_by_locus, error_by_locus=error_by_locus)
                discard_sample_state_backup(backup_dir)
    else:
        run_uce_rescue = ignore_hook

    if args.p > 1:
        avail_cpu = args.p
        asm_thr   = max(min(args.p // 2, 6), 2)
        filt_thr  = 1 if args.p < 4 else 2

        def calc_task_thr():
            """瞅着手头空闲 CPU，给下个任务匀点线程。"""
            min_thr = min(asm_thr, filt_thr) if filter_list else asm_thr
            return avail_cpu if avail_cpu - asm_thr < min_thr else asm_thr

        filter_list   = []
        refilter_list = []
        assemble_list = []

        executor      = ThreadPoolExecutor(max_workers=math.ceil(avail_cpu / filt_thr))
        running_tasks = {}
        task_metadata = {} # (stage, threads)

        if do_filter:
            filter_list.extend(reversed(samples.keys()))
        elif do_refilter:
            refilter_list.extend(reversed(samples.keys()))
        elif do_assemble:
            assemble_list.extend(reversed(samples.keys()))

        while True:
            while refilter_list and avail_cpu >= filt_thr:
                sample = refilter_list.pop()
                task_thr = calc_task_thr()
                avail_cpu -= task_thr
                running_tasks[sample] = executor.submit(run_refilter, sample, thr=task_thr)
                task_metadata[sample] = (2, task_thr)

            while assemble_list and avail_cpu >= asm_thr:
                sample = assemble_list.pop()
                task_thr = calc_task_thr()
                avail_cpu -= task_thr
                running_tasks[sample] = executor.submit(run_assembler, sample, thr=task_thr)
                task_metadata[sample] = (3, task_thr)

            while filter_list and avail_cpu >= filt_thr:
                sample = filter_list.pop()
                avail_cpu -= filt_thr
                running_tasks[sample] = executor.submit(run_filter, sample)
                task_metadata[sample] = (1, filt_thr)

            if not running_tasks:
                break

            wait(running_tasks.values(), return_when=FIRST_COMPLETED)

            processed_samples = set()

            for sample, task in running_tasks.items():
                if not task.done():
                    continue

                processed_samples.add(sample)

                stage, thr_cnt = task_metadata[sample]
                avail_cpu += thr_cnt

                try:
                    task.result()
                except Exception as e:
                    print(f'An error occurred while processing {sample}: {e}')
                    failed_samples.append((sample, {1: 'filter', 2: 'refilter', 3: 'assemble'}.get(stage, 'unknown'), str(e)))
                    continue

                if stage == 1:
                    if do_refilter:
                        refilter_list.append(sample)
                    elif do_assemble:
                        assemble_list.append(sample)
                elif stage == 2 and do_assemble:
                    assemble_list.append(sample)

            for sample in processed_samples:
                del running_tasks[sample]
                del task_metadata[sample]

        executor.shutdown()

    else:
        for name in samples.keys():
            try:
                if do_filter:
                    run_filter(name)
                if do_refilter:
                    run_refilter(name)
                if do_assemble:
                    run_assembler(name)
            except Exception as e:
                print(f'An error occurred while processing {name}: {e}')
                failed_samples.append((name, 'filter/refilter/assemble', str(e)))
                continue

    if rescue_enabled:
        rescue_samples = get_rescue_sample_names(samples, failed_samples)
        print(f'Running UCE raw-read rescue with up to {rescue_workers} sample(s) in parallel and {rescue_threads} thread(s) per sample.')

        if rescue_workers > 1:
            with ThreadPoolExecutor(max_workers=rescue_workers) as executor:
                running_rescues = {
                    executor.submit(run_uce_rescue, name, thr=rescue_threads): name
                    for name in rescue_samples
                }

                for task in as_completed(running_rescues):
                    name = running_rescues[task]

                    try:
                        task.result()
                    except Exception as e:
                        print(f'An error occurred during UCE raw-read rescue for {name}: {e}')
                        failed_samples.append((name, 'uce_rescue', str(e)))
        else:
            for name in rescue_samples:
                try:
                    run_uce_rescue(name, thr=rescue_threads)
                except Exception as e:
                    print(f'An error occurred during UCE raw-read rescue for {name}: {e}')
                    failed_samples.append((name, 'uce_rescue', str(e)))
                    continue

    write_failed_samples(out_loc, failed_samples)

    if failed_samples:
        raise RuntimeError(f'{len(failed_samples)} sample task(s) failed; see {os.path.join(out_loc, "failed_samples.tsv")}')

def make_phyluce_sample_name(sample):
    """把样本名收拾成 PHYLUCE 能认、还不犯膈应的样儿。"""
    name = ''.join(c if ord(c) < 128 and (c.isalnum() or c == '_') else '_' for c in sample).strip('_')

    if not name:
        name = 'sample'

    if not name[0].isalpha():
        name = 'sample_' + name

    return name

def get_contig_read_count(header):
    """从 contig 标题里抠出读段数，抠不着就按零算。"""
    parts = header.split('_')

    if len(parts) >= 6 and parts[0] == 'contig' and parts[5].isdigit():
        return parts[5]

    if parts and parts[-1].isdigit():
        return parts[-1]

    return '0'

def write_uce_contigs_for_phyluce(args, samples):
    """把收下的 UCE contig 改好名儿，整成 PHYLUCE 能接的文件。"""
    out_loc = args.o.strip()
    uce_dir = os.path.join(out_loc, 'uce_contigs')

    if os.path.isdir(uce_dir):
        shutil.rmtree(uce_dir, ignore_errors=True)

    os.makedirs(uce_dir, exist_ok=True)

    name_map_path = os.path.join(uce_dir, 'sample_name_map.tsv')

    with open(name_map_path, 'w', newline='') as map_file:
        writer = csv.writer(map_file, delimiter='\t')
        writer.writerow(['geneminer_sample', 'phyluce_sample', 'contigs_file', 'contig_count'])
        used_names = set()

        for sample in samples.keys():
            phyluce_sample = make_phyluce_sample_name(sample)

            if phyluce_sample in used_names:
                suffix = 2
                base_name = phyluce_sample

                while f'{base_name}_{suffix}' in used_names:
                    suffix += 1

                phyluce_sample = f'{base_name}_{suffix}'

            used_names.add(phyluce_sample)
            results_dir = os.path.join(out_loc, sample, 'results')
            summary_rows = read_uce_summary(os.path.join(out_loc, sample, 'uce_assembly_summary.csv'))
            accepted_loci = {locus for locus, row in summary_rows.items()
                             if uce_summary_row_is_accepted(row)}
            out_path = os.path.join(uce_dir, phyluce_sample + '.contigs.fasta')
            contig_count = 0

            with open(out_path, 'w') as out:
                if os.path.isdir(results_dir):
                    for entry in sorted(os.scandir(results_dir), key=lambda e: e.name):
                        if not entry.is_file() or not is_reference_file_name(entry.name):
                            continue

                        locus = os.path.splitext(entry.name)[0]
                        if locus not in accepted_loci:
                            continue

                        with open(entry.path) as f:
                            for header, seq in SimpleFastaParser(f):
                                contig_count += 1
                                read_count = get_contig_read_count(header)
                                out.write(f'>NODE_{contig_count}_length_{len(seq)}_cov_{read_count}.0_{locus}\n')
                                out.write(seq + '\n')

            writer.writerow([sample, phyluce_sample, out_path, contig_count])

def write_uce_assembly_summary(args, samples):
    """把各样本 UCE 组装表拢成一张总表。"""
    out_loc = args.o.strip()
    out_path = os.path.join(out_loc, 'uce_assembly_summary.csv')
    fieldnames = ['sample', *UCE_ASSEMBLY_SUMMARY_FIELDS]

    with open(out_path, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()

        for sample in samples.keys():
            summary_path = os.path.join(out_loc, sample, 'uce_assembly_summary.csv')

            if not os.path.isfile(summary_path):
                continue

            with open(summary_path, newline='') as f:
                for row in csv.DictReader(f):
                    row['sample'] = sample
                    writer.writerow({name: row.get(name, '') for name in fieldnames})

def write_uce_rescue_summary(args, samples):
    """把各样本救援记录并成总表，没记录就不留空壳。"""
    out_loc = args.o.strip()
    out_path = os.path.join(out_loc, 'uce_rescue_summary.csv')
    wrote_any = False

    with open(out_path, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=UCE_RESCUE_SUMMARY_FIELDS)
        writer.writeheader()

        for sample in samples.keys():
            summary_path = os.path.join(out_loc, sample, 'uce_rescue_summary.csv')

            if not os.path.isfile(summary_path):
                continue

            with open(summary_path, newline='') as f:
                for row in csv.DictReader(f):
                    writer.writerow({name: row.get(name, '') for name in UCE_RESCUE_SUMMARY_FIELDS})
                    wrote_any = True

    if not wrote_any:
        os.remove(out_path)

def write_uce_outputs(args, samples):
    """把 UCE 后续要用的 contig 和汇总产物一气儿写全。"""
    write_uce_contigs_for_phyluce(args, samples)
    write_uce_assembly_summary(args, samples)

    if args.uce_rescue_reads:
        write_uce_rescue_summary(args, samples)

def run_population(args):
    """把 population 模式参数攒齐，交给主程序开整。"""
    population_bin = find_executable('main_population', internal=True)
    engine = getattr(args, 'engine', 'pseudoref')
    command = [
        population_bin,
        '--output', args.o.strip(),
        '--samples', args.f,
        '--engine', engine,
        '--reference-strategy', args.population_reference_strategy,
        '--start-at', args.population_start_at,
        '--threads', str(args.p),
        '--min-mapq', str(args.population_min_mapq),
        '--min-baseq', str(args.population_min_baseq),
        '--min-dp', str(args.population_min_dp),
        '--min-gq', str(args.population_min_gq),
        '--min-qual', str(args.population_min_qual),
        '--min-call-rate', str(args.population_min_call_rate),
        '--min-mac', str(args.population_min_mac),
        '--ld-window', str(args.population_ld_window),
        '--ld-step', str(args.population_ld_step),
        '--ld-r2', str(args.population_ld_r2),
        '--admixture-k-min', str(args.population_admixture_k_min),
        '--admixture-k-max', str(args.population_admixture_k_max),
        '--admixture-cv', str(args.population_admixture_cv),
        '--stop-after', args.population_stop_after,
        '--minibwa', args.population_minibwa,
        '--samtools', args.population_samtools,
        '--bcftools', args.population_bcftools,
        '--plink', args.population_plink,
        '--admixture', args.population_admixture,
    ]

    if engine in ('panref', 'panrefv2'):
        command.extend(['--panref-baits', args.r])

    if getattr(args, 'population_panrefv2_include_low_confidence', False):
        command.append('--panrefv2-include-low-confidence')

    if args.population_reference_fasta:
        command.extend(['--reference-fasta', args.population_reference_fasta])

    if args.population_skip_mark_duplicates:
        command.append('--skip-mark-duplicates')

    if args.population_skip_plink:
        command.append('--skip-plink')

    if args.population_skip_admixture:
        command.append('--skip-admixture')

    subprocess.run(command, check=True)

def run_te(args):
    """Run the reference-free short-read repeatome workflow."""
    repeat_bin = find_executable('main_repeat', internal=True)
    mainfilter_bin = find_executable('MainFilterNew', internal=True)
    command = [
        repeat_bin, '--samples', args.f, '--output', args.o.strip(),
        '--stage', args.te_stage, '--threads', str(args.p),
        '--kmer', str(args.te_kmer), '--min-kmer-count', str(args.te_min_kmer_count),
        '--catalog-pairs', str(args.te_catalog_pairs), '--mainfilter', mainfilter_bin,
        '--annotation-min-fragment', str(args.te_annotate_min_fragment),
        '--annotation-max-fragment', str(args.te_annotate_max_fragment),
        '--annotation-min-support', str(args.te_annotate_min_support),
        '--annotation-min-identity', str(args.te_annotate_min_identity),
        '--annotation-min-coverage', str(args.te_annotate_min_coverage),
        '--annotation-min-delta', str(args.te_annotate_min_delta),
    ]
    if args.te_read_ledger:
        command.extend(['--read-ledger', args.te_read_ledger])
    if args.te_library:
        command.extend(['--te-library', args.te_library])
    subprocess.run(command, check=True)

def run_stats(args, samples):
    """归拢样本读段和参考信息，跑 UCE 恢复统计。"""
    stats_bin = find_executable('gm2_stats', internal=True)
    command = [stats_bin, '--output', args.o.strip(), '--reference', args.r]
    for name, reads in samples.items():
        # 单端别算两遍，整重了可不行。
        second_read = '' if reads[1] == reads[0] else reads[1]
        command.extend(['--sample', name, reads[0], second_read])
    if args.stats_count_input_reads:
        command.append('--count-input-reads')
    if args.stats_no_heatmap:
        command.append('--no-heatmap')
    subprocess.run(command, check=True)

def generate_consensus(args, samples):
    """把读段贴回组装结果，给每个基因整出一致序列。"""
    out_loc = args.o.strip()

    consensus_bin = find_executable('build_consensus', internal=True)
    minimap2_bin = find_executable('minimap2')

    if args.consensus_threshold <= 0 or args.consensus_threshold > 1:
        raise RuntimeError(f"Invalid consensus threshold {args.consensus_threshold} (must be between 0.0 and 1.0)")

    genes = get_ref_genes(args.r)

    def iterate_gene(sample):
        """把这个样本里能接着干的基因任务挨个递出来。"""
        in_dir = os.path.join(out_loc, sample, 'results')

        if not os.path.isdir(in_dir):
            print(f'Error: Sample {sample} has no assembled genes, cannot generate consensus')
            return

        cns_dir  = os.path.join(out_loc, sample, 'consensus')
        filt_dir = os.path.join(out_loc, sample, 'filtered')

        if os.path.isdir(cns_dir):
            shutil.rmtree(cns_dir, ignore_errors=True)

        os.makedirs(cns_dir, exist_ok=True)

        for name, ext in genes:
            asm_path = os.path.join(in_dir, name + ext)
            read_path = os.path.join(filt_dir, name + get_sample_ext(samples[sample][0]))

            if os.path.isfile(asm_path) and os.path.isfile(read_path):
                sam_path = os.path.join(cns_dir, name + '.sam')
                yield (name, asm_path, read_path, sam_path)

    def process_gene(task):
        """接过一个基因任务，把眼前这道工序跑利索。"""
        gene, asm_path, read_path, sam_path = task

        subprocess.run([minimap2_bin, '-ax', 'sr', '-t', '1', '--sam-hit-only', '--secondary=no',
                        '-o', sam_path, asm_path, read_path],
                       check=True)

        if os.path.isfile(sam_path):
            subprocess.run([consensus_bin, '-i', sam_path, '-c', str(args.consensus_threshold),
                            '-o', os.path.dirname(sam_path), '-s', '0'],
                           check=True)

            os.remove(sam_path)

    if args.p > 1:
        with ThreadPoolExecutor(max_workers=args.p) as executor:
            for _ in executor.map(process_gene, (task
                                                for sample in samples.keys()
                                                for task in iterate_gene(sample))):
                pass

    else:
        for sample in samples.keys():
            for task in iterate_gene(sample):
                process_gene(task)

def blast_trim(args, samples):
    """拿 BLAST 对照参考，把组装序列两边收拾利索。"""
    out_loc = args.o.strip()
    trim_bin = find_executable('build_trimed', internal=True)

    makeblastdb_bin = find_executable('makeblastdb')

    if args.trim_mode == 'isoform':
        blast_bin = find_executable('magicblast')
    else:
        blast_bin = find_executable('blastn')

    if args.trim_retention < 0 or args.trim_retention > 1:
        raise RuntimeError(f"Invalid trim retention threshold {args.trim_retention} (must be between 0.0 and 1.0)")

    trim_modes = {'all': '0', 'longest': '1', 'terminal': '2', 'isoform': '3'}
    trim_mode = trim_modes[args.trim_mode]

    genes = get_ref_genes(args.r)

    os.makedirs(os.path.join(out_loc, 'blast_db'), exist_ok=True)

    def build_blast_db(name_tup):
        """给一个参考基因建 BLAST 库，后头查着快。"""
        name, ext = name_tup
        subprocess.run([makeblastdb_bin, "-in", os.path.realpath(os.path.join(args.r, name + ext)),
                        "-dbtype", "nucl", "-out", name],
                       cwd=os.path.join(out_loc, 'blast_db'), check=True)

    def iterate_gene(sample):
        """把这个样本里能接着干的基因任务挨个递出来。"""
        if args.trim_source == 'consensus':
            in_dir = os.path.join(out_loc, sample, 'consensus')
        else:
            in_dir = os.path.join(out_loc, sample, 'results')

        if not os.path.isdir(in_dir):
            print(f'Error: Sample {sample} has no {args.trim_source} sequences, cannot trim')
            return

        blast_dir = os.path.join(out_loc, sample, 'blast')

        if os.path.isdir(blast_dir):
            shutil.rmtree(blast_dir, ignore_errors=True)

        os.makedirs(blast_dir, exist_ok=True)

        for name, ext in genes:
            asm_path = os.path.join(in_dir, name + '.fasta')
            ref_path = os.path.join(args.r, name + ext)

            if os.path.isfile(asm_path):
                yield (name, asm_path, ref_path, os.path.join(blast_dir, name + '.fasta'))

    def process_gene(task):
        """接过一个基因任务，把眼前这道工序跑利索。"""
        name, asm_path, ref_path, out_path = task
        subprocess.run([
            trim_bin,
            '-i', asm_path,
            '-r', ref_path,
            '-o', out_path,
            '-b', os.path.join(out_loc, 'blast_db', name),
            '-m', trim_mode,
            '-p', str(args.trim_retention * 100),
            '--executable', blast_bin,
        ], check=True)

    gene_count = len(genes) * len(samples)
    trimmed_count = 0

    if args.p > 1:
        with ThreadPoolExecutor(max_workers=args.p) as executor:
            for _ in executor.map(build_blast_db, genes):
                pass

            for _ in executor.map(process_gene, (task
                                                for sample in samples.keys()
                                                for task in iterate_gene(sample))):
                trimmed_count += 1

                if trimmed_count >= 2:
                    print(f'{trimmed_count}/{gene_count} genes trimmed\r', end='')

    else:
        for gene_tup in genes:
            build_blast_db(gene_tup)

        for sample in samples.keys():
            for task in iterate_gene(sample):
                process_gene(task)

                trimmed_count += 1

                if trimmed_count >= 2:
                    print(f'{trimmed_count}/{gene_count} genes trimmed\r', end='')

    print('\n')

def combine_genes(args, samples):
    """按 locus 合样本、做比对再清理，整出能拼接的结果。"""
    out_loc = args.o.strip()
    alignment_filter = get_alignment_filter(args)

    if not args.no_alignment:
        msa_threads = get_msa_threads(args)
        filter_processes = get_filter_processes(args)

        alifilter_model = get_alifilter_model(args)

        if alifilter_model and alignment_filter != 'alifilter':
            raise RuntimeError("--alifilter-model requires --alignment-filter alifilter")

        if args.msa_program == 'clustalo':
            msa_bin = find_executable('clustalo')
        else:
            msa_bin = find_executable('mafft')

        if alignment_filter == 'trimal':
            trimal_bin = find_executable('trimal')
            alifilter_bin = None
        elif alignment_filter == 'alifilter':
            trimal_bin = None
            alifilter_bin = find_executable('AliFilter')
        else:
            trimal_bin = None
            alifilter_bin = None

        merge_seq_bin = find_executable('merge_seq', internal=True)
        alignment_cleaner_bin = find_executable('fix_alignment', internal=True)

        if args.clean_difference < 0 or args.clean_difference > 1:
            raise RuntimeError(f"Invalid maximum difference {args.clean_difference} (must be between 0.0 and 1.0)")

        if args.clean_sequences < 0 or args.clean_sequences > len(samples):
            raise RuntimeError(f"Invalid required number of sequences {args.clean_sequences} (must be between 0 and {len(samples)})")

    combine_dir = os.path.join(out_loc, 'combined_results')

    if os.path.isdir(combine_dir):
        shutil.rmtree(combine_dir, ignore_errors=True)

    os.makedirs(combine_dir, exist_ok=True)

    if not args.no_alignment:
        alignment_dir = os.path.join(out_loc, 'combined_results', 'aligned')
        trim_dir = os.path.join(out_loc, 'combined_trimed')
        combined_fasta = os.path.join(out_loc, 'combined_results.fasta')
        trim_fasta = os.path.join(out_loc, 'combined_trimed.fasta')

        if os.path.isdir(trim_dir):
            shutil.rmtree(trim_dir, ignore_errors=True)

        if os.path.isfile(trim_fasta):
            os.remove(trim_fasta)

        if os.path.isfile(combined_fasta):
            os.remove(combined_fasta)

        os.makedirs(alignment_dir, exist_ok=True)

        if alignment_filter != 'none':
            os.makedirs(trim_dir, exist_ok=True)

    if args.combine_source == 'trimmed':
        in_name = 'blast'
    elif args.combine_source == 'consensus':
        in_name = 'consensus'
    else:
        in_name = 'results'

    genes = {t[0] for t in get_ref_genes(args.r)}
    accepted_loci_by_sample = {}
    if getattr(args, 'assembly_mode', 'original') == 'uce':
        for name in samples.keys():
            summary_rows = read_uce_summary(
                os.path.join(out_loc, name, 'uce_assembly_summary.csv'))
            accepted_loci_by_sample[name] = {
                locus for locus, row in summary_rows.items()
                if uce_summary_row_is_accepted(row)
            }

    def merge_gene(gene):
        """把各样本这个基因的序列归到一个 FASTA 里。"""
        out_path = os.path.join(combine_dir, gene + '.fasta')
        written  = False

        with open(out_path, 'w+') as f:
            for name in samples.keys():
                if (accepted_loci_by_sample
                        and gene not in accepted_loci_by_sample.get(name, set())):
                    continue

                in_path = os.path.join(out_loc, name, in_name, gene + '.fasta')

                if not os.path.isfile(in_path):
                    continue

                with open(in_path, 'r') as r:
                    try:
                        _, seq = next(SimpleFastaParser(r))
                    except StopIteration:
                        continue

                    f.write(f'>{name}\n{seq}\n')

                written = True

        if not written:
            os.remove(out_path)

    def handle_locus_failure(gene, stage, exc=None, output_path=None):
        """单个 locus 整失败了就收拾残局，再按严格模式处理。"""
        if output_path and os.path.isfile(output_path):
            os.remove(output_path)

        message = f"{stage} failed on {gene}"

        if args.strict_combine_errors:
            raise RuntimeError(message) from exc

        if exc:
            message = f"{message}: {exc}"

        print(f"Warning: {message}", file=sys.stderr)
        return False

    def align_gene(gene):
        """给这个基因跑多序列比对，整不成就稳当报错。"""
        in_path = os.path.join(combine_dir, gene + '.fasta')
        out_path = os.path.join(alignment_dir, gene + '.fasta')

        if not os.path.isfile(in_path):
            return False

        try:
            if args.msa_program == 'clustalo':
                subprocess.run([msa_bin, '-i', in_path, '-o', out_path, '--auto', '--force',
                                '--seqtype=DNA', f'--threads={msa_threads}'],
                               stderr=subprocess.DEVNULL, check=True)
            else:
                with open(out_path, 'w') as out:
                    subprocess.run([msa_bin, '--auto', '--quiet', '--nuc',
                                    '--thread', str(msa_threads), in_path],
                                   stdout=out, stderr=subprocess.DEVNULL, check=True)
        except (OSError, subprocess.CalledProcessError, RuntimeError, ValueError) as e:
            return handle_locus_failure(gene, args.msa_program, e, out_path)

        return os.path.isfile(out_path)

    def clean_gene(gene):
        """把比对里差得太离谱的序列挑出去。"""
        gene_path = os.path.join(alignment_dir, gene + '.fasta')

        if os.path.isfile(gene_path):
            try:
                subprocess.run([
                    alignment_cleaner_bin,
                    '-f', gene_path,
                    '-n', str(args.clean_sequences),
                    '-p', str(args.clean_difference),
                ], check=True)
            except (OSError, subprocess.CalledProcessError) as e:
                return handle_locus_failure(gene, 'alignment cleanup', e, gene_path)

        return os.path.isfile(gene_path)

    def filter_gene(gene):
        """拿选定工具把不靠谱的比对位点修剪掉。"""
        in_path = os.path.join(alignment_dir, gene + '.fasta')
        out_path = os.path.join(trim_dir, gene + '.fasta')

        if not os.path.isfile(in_path):
            return False

        if alignment_filter == 'trimal':
            cmd = [trimal_bin, '-in', in_path, '-out', out_path, '-automated1']
        elif alignment_filter == 'alifilter':
            cmd = [alifilter_bin, '-i', in_path, '-o', out_path]

            if alifilter_model:
                cmd.extend(['-m', alifilter_model])
        else:
            return False

        try:
            subprocess.run(cmd, check=True)
        except (OSError, subprocess.CalledProcessError) as e:
            return handle_locus_failure(gene, alignment_filter, e, out_path)

        return os.path.isfile(out_path)

    alignment_count = 0
    gene_count = len(genes)

    if args.no_alignment:
        process_gene = merge_gene
    else:
        msa_slots = max(1, args.p // msa_threads)
        msa_semaphore = threading.Semaphore(msa_slots)
        filter_semaphore = threading.Semaphore(filter_processes)

        def process_gene(gene):
            """这个基因从合并到比对、清理、修剪一趟整完。"""
            merge_gene(gene)

            with msa_semaphore:
                aligned = align_gene(gene)

            if not aligned:
                return False

            cleaned = clean_gene(gene)

            if not cleaned:
                return False

            if alignment_filter != 'none':
                with filter_semaphore:
                    filtered = filter_gene(gene)

                if not filtered:
                    return False

            return True

    if args.p > 1:
        with ThreadPoolExecutor(max_workers=args.p) as executor:
            for aligned in executor.map(process_gene, genes):
                if aligned:
                    alignment_count += 1

                    if alignment_count >= 2:
                        print(f'{alignment_count}/{gene_count} genes aligned\r', end='')

    else:
        for gene in genes:
            aligned = process_gene(gene)
            if aligned:
                alignment_count += 1

                if alignment_count >= 2:
                    print(f'{alignment_count}/{gene_count} genes aligned\r', end='')

    print('\n')

    if not args.no_alignment:
        subprocess.run([merge_seq_bin, '-input', alignment_dir, '-exts', '.fasta', '-missing', '-',
                        '-output', os.path.join(out_loc, 'combined_results.fasta')],
                       check=True)

        if alignment_filter != 'none':
            subprocess.run([merge_seq_bin, '-input', trim_dir, '-exts', '.fasta', '-missing', '-',
                            '-output', os.path.join(out_loc, 'combined_trimed.fasta')],
                           check=True)

def get_alignment_filter(args):
    """瞅参数定下比对过滤工具，老参数也给它兜着。"""
    if getattr(args, 'no_trimal', False):
        return 'none'

    return getattr(args, 'alignment_filter', None) or 'trimal'

def get_alifilter_model(args):
    """把 AliFilter 模型名儿收拾明白，默认模型就不额外传。"""
    model = getattr(args, 'alifilter_model', None)

    if not model:
        return None

    model = model.strip()

    if not model or model.lower() == 'default':
        return None

    return model

def get_msa_threads(args):
    """把每个多序列比对该用几条线程定下来。"""
    return max(1, getattr(args, 'msa_threads', 1))

def get_filter_processes(args):
    """算明白比对过滤最多能同时跑几个。"""
    filter_processes = getattr(args, 'filter_processes', None)

    if filter_processes is None:
        return max(1, args.p)

    return max(1, filter_processes)

def get_locus_alignment_dir(args):
    """按过滤设置找准单 locus 比对结果搁的地方。"""
    out_loc = args.o.strip()

    if get_alignment_filter(args) == 'none':
        return os.path.join(out_loc, 'combined_results', 'aligned')

    return os.path.join(out_loc, 'combined_trimed')

def get_concatenated_alignment_path(args):
    """按过滤设置找准拼接比对文件搁的地方。"""
    out_loc = args.o.strip()

    if get_alignment_filter(args) == 'none':
        return os.path.join(out_loc, 'combined_results.fasta')

    return os.path.join(out_loc, 'combined_trimed.fasta')

def build_single_tree(prog_name, prog_bin, in_path, bootstrap=0, quiet=False, threads=1):
    """照选定的建树程序给一份比对整出一棵树。"""
    if prog_name == 'raxmlng':
        out_path = in_path + ".raxml.bestTree"
        params = [prog_bin, '--msa', in_path, '--msa-format', 'FASTA',
                  '--model', 'GTR+G', '--redo']

        if bootstrap:
            params.extend(['--all', '--bs-trees', str(bootstrap)])
        else:
            params.append('--search')

        if threads > 1:
            params.extend(['--threads', f'auto{{{threads}}}', '--workers', 'auto'])
        else:
            params.extend(['--threads', '1'])

        if os.path.isfile(out_path):
            os.remove(out_path)

        subprocess.run(params, stdout=subprocess.DEVNULL if quiet else None, check=True)

        return out_path

    elif prog_name == 'iqtree':
        out_path = in_path + ".treefile"
        params = [prog_bin, '-s', in_path, '-redo']

        if bootstrap:
            params.extend(['-B', str(bootstrap)])

        if quiet:
            params.append('-quiet')

        if threads > 1:
            params.extend(['-T', 'AUTO', '-ntmax', str(threads)])
        else:
            params.extend(['-T', '1'])

        if os.path.isfile(out_path):
            os.remove(out_path)

        subprocess.run(params, stdout=subprocess.DEVNULL if quiet else None, check=True)

        return out_path

    elif prog_name == 'veryfasttree':
        out_path = in_path + ".veryfasttree.tre"
        params = [prog_bin, '-out', out_path, '-gtr']

        if bootstrap:
            params.extend(['-boot', str(bootstrap)])
        else:
            params.append('-nosupport')

        if quiet:
            params.extend(['-quiet'])

        if threads > 1:
            params.extend(['-threads', str(threads)])

        params.extend(['-nt', in_path])

        if os.path.isfile(out_path):
            os.remove(out_path)

        subprocess.run(params, stderr=subprocess.DEVNULL if quiet else None, check=True)

        return out_path

    else:
        out_path = in_path + ".fasttree.tre"
        params = [prog_bin, '-out', out_path, '-gtr']

        if bootstrap:
            params.extend(['-boot', str(bootstrap)])
        else:
            params.append('-nosupport')

        if quiet:
            params.append('-quiet')

        params.extend(['-nt', in_path])

        if os.path.isfile(out_path):
            os.remove(out_path)

        subprocess.run(params, stderr=subprocess.DEVNULL if quiet else None, check=True)

        return out_path

def run_gene_tree_job(gene, alignment_dir, make_gene_tree):
    """给单个基因跑树任务，成败都规整地交回来。"""
    in_path = os.path.join(alignment_dir, f'{gene}.fasta')

    try:
        return gene, in_path, make_gene_tree(gene), ''
    except Exception as e:
        return gene, in_path, None, str(e)

def write_failed_gene_trees(out_loc, failures):
    """把没整成的基因树记下来，没失败就撤掉旧表。"""
    out_path = os.path.join(out_loc, 'failed_gene_trees.tsv')

    if not failures:
        if os.path.isfile(out_path):
            os.remove(out_path)

        return

    with open(out_path, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['locus', 'alignment', 'error'])
        writer.writerows(failures)

def build_coalescent_tree(args):
    """先逐基因建树，再用 ASTRAL 拢成共祖树。"""
    out_loc = args.o.strip()

    if args.phylo_program == 'raxmlng':
        phylo_bin = find_executable('raxml-ng')
    elif args.phylo_program == 'iqtree':
        phylo_bin = find_executable('iqtree')
    elif args.phylo_program == 'veryfasttree':
        phylo_bin = find_executable('VeryFastTree')
    else:
        phylo_bin = find_executable('FastTree')

    astral_bin = find_executable('astral')

    def find_genes(path):
        """到指定目录瞅瞅现成比对里都有啥基因。"""
        try:
            with os.scandir(path) as it:
                return {os.path.splitext(entry.name)[0] for entry in it if entry.is_file() and is_reference_file_name(entry.name)}
        except OSError:
            return set()

    alignment_dir = get_locus_alignment_dir(args)

    genes = {t[0] for t in get_ref_genes(args.r)} & find_genes(alignment_dir)
    gene_count = len(genes)

    if not genes:
        raise RuntimeError(f"No gene alignments found under '{alignment_dir}'")

    def make_gene_tree(gene):
        """给这个基因悄默声地跑出一棵树。"""
        return build_single_tree(args.phylo_program, phylo_bin, os.path.join(alignment_dir, f'{gene}.fasta'), quiet=True)

    tree_files = set()

    failed_gene_trees = []

    def handle_gene_tree_result(result):
        """把基因树结果分成成功和失败两摞，别整混了。"""
        gene, alignment_path, tree_path, error = result

        if error:
            print(f'Warning: gene tree failed on {gene}: {error}')
            failed_gene_trees.append((gene, alignment_path, error))
            return

        if tree_path and os.path.isfile(tree_path):
            tree_files.add(tree_path)
            tree_count = len(tree_files)

            if tree_count >= 2:
                print(f'{tree_count}/{gene_count} trees built\r', end='')

    if args.p > 1:
        with ThreadPoolExecutor(max_workers=args.p) as executor:
            futures = [executor.submit(run_gene_tree_job, gene, alignment_dir, make_gene_tree) for gene in sorted(genes)]

            for task in as_completed(futures):
                handle_gene_tree_result(task.result())

    else:
        for gene in sorted(genes):
            handle_gene_tree_result(run_gene_tree_job(gene, alignment_dir, make_gene_tree))

    print('\n')
    write_failed_gene_trees(out_loc, failed_gene_trees)

    coal_trees_path = os.path.join(out_loc, 'combined_genes.trees')
    coal_out_path = os.path.join(out_loc, 'Coalescent.tree')
    written = False

    with open(coal_trees_path, 'w') as f:
        # 顺序排明白，回回生成都一个样。
        for path in sorted(tree_files):
            if os.path.getsize(path) <= 2: # Empty tree
                continue

            with open(path, 'r') as r:
                tree = next((line.strip() for line in r if line.strip()), '')

            if not tree:
                continue

            f.write(tree + '\n')

            written = True

    if not written:
        raise RuntimeError(f"Unable to reconstruct coalescent trees because no gene tree is available")

    if os.path.isfile(coal_out_path):
        os.remove(coal_out_path)

    subprocess.run([astral_bin, '-i', coal_trees_path, '-o', coal_out_path, '-t', str(args.p)], check=True)

def build_concatenation_tree(args):
    """拿拼接好的全套比对直接整一棵物种树。"""
    out_loc = args.o.strip()

    if args.phylo_program == 'raxmlng':
        phylo_bin = find_executable('raxml-ng')
    elif args.phylo_program == 'iqtree':
        phylo_bin = find_executable('iqtree')
    elif args.phylo_program == 'veryfasttree':
        phylo_bin = find_executable('VeryFastTree')
    else:
        phylo_bin = find_executable('FastTree')

    in_path = get_concatenated_alignment_path(args)

    if not os.path.isfile(in_path):
        raise RuntimeError(f"Unable to find the concatenated alignment at '{in_path}'")

    final_tree_path = os.path.join(out_loc, 'Concatenation.tree')

    if os.path.isfile(final_tree_path):
        os.remove(final_tree_path)

    out_path = build_single_tree(args.phylo_program, phylo_bin, in_path,
                                 bootstrap=args.bootstrap, threads=args.p)

    if not os.path.isfile(out_path):
        raise RuntimeError(f"Phylogenetic tree reconstruction failed")

    shutil.copyfile(out_path, final_tree_path)


def read_single_fasta_sequence(path):
    with open(path) as handle:
        return ''.join(line.strip() for line in handle if not line.startswith('>'))

def circularly_consistent(left, right, minimum_identity=0.9999):
    """Return true when two circular sequences have the same cut-independent sequence.

    Exact equality remains the normal fast path.  The seeded comparison also accepts
    the very small number of consensus differences expected between consecutive
    read-depth stages, without treating a materially different assembly as stable.
    """
    if len(left) != len(right) or not left:
        return False
    doubled = right + right
    if left in doubled:
        return True
    seed_length = min(31, len(left))
    anchors = sorted({0, len(left) // 4, len(left) // 2, len(left) * 3 // 4})
    starts = set()
    for anchor in anchors:
        seed = left[anchor:anchor + seed_length]
        start = doubled.find(seed)
        while 0 <= start < len(right):
            starts.add((start - anchor) % len(right))
            start = doubled.find(seed, start + 1)
    for start in starts:
        candidate = doubled[start:start + len(left)]
        matches = sum(base == candidate_base for base, candidate_base in zip(left, candidate))
        if matches / len(left) >= minimum_identity:
            return True
    return False

def mito_stage_is_circular(path):
    try:
        with open(path) as handle:
            return any(line.strip() == 'status\tcircular' for line in handle)
    except OSError:
        return False

def run_mito_adaptive(args, samples):
    base_output = args.o.strip()
    stage_root = os.path.join(base_output, '.mito_adaptive')
    reference_dir = args.r
    original_max_reads = args.max_reads
    previous = None
    limit = args.mito_initial_reads
    maximum = args.mito_max_reads
    while True:
        stage_output = os.path.join(stage_root, f'{limit}m')
        args.o = stage_output
        args.max_reads = limit
        do_filter_assemble(args, samples, True, True, True)
        write_uce_outputs(args, samples)
        run_mito_finalize(args, samples)
        current = {}
        for sample in samples:
            sample_dir = os.path.join(stage_output, sample, 'mito')
            summary = os.path.join(sample_dir, 'mitochondrial_assembly_summary.tsv')
            fasta = os.path.join(sample_dir, 'mitochondrial_assembly.fasta')
            if not mito_stage_is_circular(summary):
                current[sample] = None
            else:
                current[sample] = read_single_fasta_sequence(fasta)
        if previous and all(current.get(sample) and circularly_consistent(previous[sample], current[sample]) for sample in samples):
            for sample in samples:
                source = os.path.join(stage_output, sample)
                destination = os.path.join(base_output, sample)
                if os.path.isdir(destination):
                    shutil.rmtree(destination)
                shutil.copytree(source, destination)
            args.o = base_output
            args.max_reads = original_max_reads
            args.r = reference_dir
            return
        previous = current
        if maximum and limit >= maximum:
            raise RuntimeError(f'mito adaptive stop did not confirm a circular assembly by {limit}M reads')
        next_limit = limit * 2
        limit = min(next_limit, maximum) if maximum else next_limit

def run_gene_finalize(args, samples):
    """Summarize original-rust candidate contigs without altering assembly behavior."""
    gene_bin = find_executable('gene_workflow', internal=True)
    gene_root = os.path.join(args.o.strip(), 'gene')
    if os.path.isdir(gene_root):
        shutil.rmtree(gene_root, ignore_errors=True)
    os.makedirs(gene_root, exist_ok=True)

    def classify_sample(sample):
        sample_dir = os.path.join(args.o.strip(), sample)
        command = [
            gene_bin, 'classify', '--reference', args.r,
            '--contigs', os.path.join(sample_dir, 'contigs_all'),
            '--sample', sample, '--out', gene_root,
        ]
        subprocess.run(command, check=True)

    ordered_samples = sorted(samples)
    workers = min(max(1, args.p), len(ordered_samples))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(classify_sample, sample): sample for sample in ordered_samples}
        for future in as_completed(futures):
            sample = futures[future]
            try:
                future.result()
            except Exception as error:
                raise RuntimeError(f'gene candidate classification failed for {sample}: {error}') from error

    command = [gene_bin, 'cohort', '--reference', args.r, '--out', gene_root]
    for sample in ordered_samples:
        command.extend(['--sample', sample])
    subprocess.run(command, check=True)
    if args.gene_protein_reference:
        if not os.path.isdir(args.gene_protein_reference):
            raise RuntimeError('--gene-protein-reference must be a directory')
        annotation_dir = os.path.join(args.o.strip(), 'gene_annotation')
        if os.path.isdir(annotation_dir):
            shutil.rmtree(annotation_dir, ignore_errors=True)
        subprocess.run([gene_bin, 'annotate', '--input', gene_root, '--protein-reference', args.gene_protein_reference, '--out', annotation_dir, '--miniprot', args.gene_miniprot, '--threads', str(max(1, args.p))], check=True)


def run_gene_annotate(args):
    gene_bin = find_executable('gene_workflow', internal=True)
    if not os.path.isdir(args.gene_input):
        raise RuntimeError('--gene-input must be a gene output directory')
    if not os.path.isdir(args.gene_protein_reference):
        raise RuntimeError('gene-annotate requires --gene-protein-reference')
    subprocess.run([gene_bin, 'annotate', '--input', args.gene_input, '--protein-reference', args.gene_protein_reference, '--out', args.o.strip(), '--miniprot', args.gene_miniprot, '--threads', str(max(1, args.p))], check=True)


def run_gene_resolve(args):
    gene_bin = find_executable('gene_workflow', internal=True)
    if not os.path.isdir(args.gene_input):
        raise RuntimeError('--gene-input must be an annotation directory')
    command = [gene_bin, 'resolve', '--input', args.gene_input, '--out', args.o.strip(), '--mafft', args.gene_mafft, '--iqtree', args.gene_iqtree, '--threads', str(max(1, args.p)), '--min-taxa', str(args.gene_min_taxa)]
    if args.gene_outgroup:
        if not os.path.isfile(args.gene_outgroup):
            raise RuntimeError('--gene-outgroup must be a readable file')
        command.extend(['--outgroup', args.gene_outgroup])
    if args.gene_ufboot:
        command.extend(['--ufboot', str(args.gene_ufboot)])
    if args.gene_taper:
        if not os.path.isfile(args.gene_taper):
            raise RuntimeError('--gene-taper must be a readable correction_multi.jl script')
        command.extend(['--taper-script', args.gene_taper, '--julia', args.gene_julia])
    subprocess.run(command, check=True)


def gene_file_sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def run_gene_tree(args):
    resolved = args.gene_input
    if not os.path.isdir(resolved):
        raise RuntimeError('--gene-input must be a gene-resolve output directory')
    if args.gene_species_mode == 'strict':
        trees = os.path.join(resolved, 'astral_input', 'resolved_1to1.trees')
        output_name = 'gene_strict_aster.tree'
    else:
        trees = os.path.join(resolved, 'astralpro_input', 'multicopy.trees')
        mapping = os.path.join(resolved, 'astralpro_input', 'leaf_to_species.tsv')
        output_name = 'gene_multicopy_aster.tree'
        if not os.path.isfile(mapping):
            raise RuntimeError(f'Missing multicopy leaf map: {mapping}')
    if not os.path.isfile(trees) or os.path.getsize(trees) == 0:
        raise RuntimeError(f'No usable {args.gene_species_mode} gene trees found: {trees}')
    os.makedirs(args.o.strip(), exist_ok=True)
    executable = shutil.which(args.gene_aster) if not os.path.isabs(args.gene_aster) else args.gene_aster
    if not executable or not os.path.isfile(executable):
        raise RuntimeError(f'Cannot find ASTER2 executable: {args.gene_aster}')
    output = os.path.join(args.o.strip(), output_name)
    log = os.path.join(args.o.strip(), output_name + '.log')
    if os.path.isfile(output):
        os.remove(output)
    command = [executable, '-i', trees, '-o', output, '-t', str(max(1, args.p))]
    if args.gene_species_mode == 'multicopy':
        command.extend(['-a', mapping])
    with open(log, 'w') as log_handle:
        subprocess.run(command, check=True, stdout=log_handle, stderr=subprocess.STDOUT)
    if not os.path.isfile(output) or os.path.getsize(output) == 0:
        raise RuntimeError(f'ASTER2 completed without a species tree; inspect {log}')
    with open(output, 'r') as tree_handle:
        tree = next((line.strip() for line in tree_handle if line.strip()), '')
    if not tree.startswith('(') or not tree.endswith(';'):
        raise RuntimeError(f'ASTER2 output is not a Newick tree; inspect {log}')
    provenance = os.path.join(args.o.strip(), 'gene_tree_provenance.tsv')
    with open(provenance, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(('field', 'value'))
        writer.writerow(('mode', args.gene_species_mode))
        writer.writerow(('aster_executable', os.path.realpath(executable)))
        writer.writerow(('command', shlex.join(command)))
        writer.writerow(('gene_trees', os.path.realpath(trees)))
        writer.writerow(('gene_trees_sha256', gene_file_sha256(trees)))
        if args.gene_species_mode == 'multicopy':
            writer.writerow(('leaf_to_species', os.path.realpath(mapping)))
            writer.writerow(('leaf_to_species_sha256', gene_file_sha256(mapping)))
        writer.writerow(('species_tree', os.path.realpath(output)))
        writer.writerow(('species_tree_sha256', gene_file_sha256(output)))


def execute_tasks(args, samples):
    """照命令顺序调度整条流程，哪步出岔子都稳当收口。"""
    if 'te' not in args.command and 'gene-tree' not in args.command and 'gene-resolve' not in args.command and 'gene-annotate' not in args.command and not os.path.isdir(args.r):
        print(f"Reference directory '{args.r}' does not exist")
        return 2

    commands = frozenset(args.command)

    do_profile = 'profiling' in commands
    do_filter = 'filter' in commands
    do_refilter = 'refilter' in commands
    do_assemble = 'assemble' in commands
    do_population = 'population' in commands
    do_te = 'te' in commands
    do_consensus = 'consensus' in commands
    do_trim = 'trim' in commands
    do_combine = 'combine' in commands
    do_tree = 'tree' in commands
    do_stats = 'stats' in commands
    do_gene = bool(getattr(args, 'is_gene_workflow', False))
    do_gene_resolve = 'gene-resolve' in commands
    do_gene_annotate = 'gene-annotate' in commands
    do_gene_tree = 'gene-tree' in commands

    try:
        if do_gene_annotate:
            run_gene_annotate(args)
            return 0
        if do_gene_resolve:
            run_gene_resolve(args)
            return 0
        if do_gene_tree:
            run_gene_tree(args)
            return 0
        if getattr(args, 'is_mito_workflow', False) and getattr(args, 'mito_adaptive_stop', False):
            run_mito_adaptive(args, samples)
            return 0
        if do_profile:
            if len(commands) != 1:
                raise RuntimeError('profiling is a complete marker workflow and cannot be combined with other subcommands')
            do_filter_assemble(args, samples, True, False, True)
        elif do_filter or do_refilter or do_assemble:
            do_filter_assemble(args, samples, do_filter, do_refilter, do_assemble)

            if do_assemble and args.assembly_mode == 'uce':
                write_uce_outputs(args, samples)

            if do_assemble and getattr(args, 'is_mito_workflow', False):
                run_mito_finalize(args, samples)

            if do_assemble and do_gene:
                run_gene_finalize(args, samples)

        if do_te:
            if len(commands) != 1:
                raise RuntimeError('te is a complete workflow and cannot be combined with other subcommands')
            run_te(args)
            return 0

        if do_population:
            if do_assemble and args.assembly_mode != 'uce':
                raise RuntimeError('population with assemble requires --assembly-mode uce')

            run_population(args)

        if do_consensus:
            generate_consensus(args, samples)

        if do_trim:
            if not args.trim_source:
                args.trim_source = 'consensus' if do_consensus else 'assembly'

            blast_trim(args, samples)

        if do_combine:
            if not args.combine_source:
                if do_trim:
                    args.combine_source = 'trimmed'
                elif do_consensus:
                    args.combine_source = 'consensus'
                else:
                    args.combine_source = 'assembly'

            combine_genes(args, samples)

        if do_tree:
            if args.tree_method == 'coalescent':
                build_coalescent_tree(args)
            else:
                build_concatenation_tree(args)

        if do_stats:
            run_stats(args, samples)

    # 文件和外部工具报错都兜住，别甩一屏 traceback。
    except (OSError, ValueError, csv.Error, RuntimeError, subprocess.SubprocessError) as e:
        print(f'Error: {e}')
        return 1

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description="GeneMiner2-UCE extracts phylogenetic marker loci for UCE workflows.",
                                     epilog=HELP_EPILOG)
    parser.add_argument('command',
                        choices=('filter', 'refilter', 'assemble', 'mito', 'gene', 'gene-annotate', 'gene-resolve', 'gene-tree', 'profiling', 'population', 'te', 'consensus', 'trim', 'combine', 'tree', 'stats'),
                        help='One or several of the following actions, separated by space:' + COMMAND_HELP,
                        metavar='command',
                        nargs='*')

    group_io = parser.add_argument_group('input/output parameters')
    group_io.add_argument('-f', help='Sample list file', metavar='FILE', required=False, default='')
    group_io.add_argument('-r', help='Reference directory; optional with mito, which derives it from --mito-genbank', metavar='DIR', required=False, default='')
    group_io.add_argument('--mito-genbank', default='', help='Annotated mitochondrial GenBank reference; enables the mito workflow', metavar='FILE')
    group_io.add_argument('--mito-flank', default=150, type=int, help=argparse.SUPPRESS, metavar='INT')
    group_io.add_argument('--mito-tile-length', default=1200, type=int, help=argparse.SUPPRESS, metavar='INT')
    group_io.add_argument('--mito-tile-step', default=600, type=int, help=argparse.SUPPRESS, metavar='INT')
    group_io.add_argument('--mito-min-overlap', default=41, type=int, help=argparse.SUPPRESS, metavar='INT')
    group_io.add_argument('--mito-min-overlap-identity', default=0.98, type=float, help=argparse.SUPPRESS, metavar='FLOAT')
    group_io.add_argument('--mito-min-junction-support', default=3, type=int, help=argparse.SUPPRESS, metavar='INT')
    group_io.add_argument('--mito-terminal-window', default=500, type=int, help=argparse.SUPPRESS, metavar='INT')
    group_io.add_argument('--mito-link-kmer', default=31, type=int, help=argparse.SUPPRESS, metavar='INT')
    group_io.add_argument('--mito-min-link-hits', default=2, type=int, help=argparse.SUPPRESS, metavar='INT')
    group_io.add_argument('--mito-min-pair-support', default=3, type=int, help=argparse.SUPPRESS, metavar='INT')
    group_io.add_argument('--mito-bridge-kmer', default=31, type=int, help=argparse.SUPPRESS, metavar='INT')
    group_io.add_argument('--mito-bridge-min-depth', default=2, type=int, help=argparse.SUPPRESS, metavar='INT')
    group_io.add_argument('--mito-max-bridge', default=1000, type=int, help=argparse.SUPPRESS, metavar='INT')
    group_io.add_argument('--mito-initial-reads', default=10, type=int, help=argparse.SUPPRESS, metavar='INT')
    group_io.add_argument('--mito-max-reads', default=320, type=int, help='Maximum adaptive input limit in approximately 1.05M paired-read blocks (default = 320)', metavar='INT')
    group_io.add_argument('--no-mito-adaptive-stop', dest='mito_adaptive_stop', action='store_false', default=True, help='Disable default staged mito adaptive stopping')
    group_io.add_argument('-o', help='Output directory', metavar='DIR', required=True)
    group_io.add_argument('-p', default=1, help='Number of parallel processes', metavar='INT', type=int)
    group_io.add_argument('--gene-protein-reference', default='', help='Optional per-family protein FASTA directory; enables gene annotation', metavar='DIR')
    group_io.add_argument('--gene-miniprot', default='miniprot', help='miniprot executable for --gene-protein-reference', metavar='FILE')
    group_io.add_argument('--gene-input', default='', help='Gene workflow input directory (annotation for gene-resolve; resolve output for gene-tree)', metavar='DIR')
    group_io.add_argument('--gene-mafft', default='mafft', help='MAFFT executable for gene-resolve', metavar='FILE')
    group_io.add_argument('--gene-iqtree', default='iqtree', help='IQ-TREE executable for gene-resolve', metavar='FILE')
    group_io.add_argument('--gene-min-taxa', default=4, type=int, help='Minimum taxa per resolved SC-OG (default = 4)', metavar='INT')
    group_io.add_argument('--gene-outgroup', default='', help='Optional TSV/CSV whose first column lists outgroup samples', metavar='FILE')
    group_io.add_argument('--gene-ufboot', default=0, type=int, help='Optional IQ-TREE UFBoot replicates; use 0 or >=1000 (default = 0)', metavar='INT')
    group_io.add_argument('--gene-taper', default='', help='Optional TAPER correction_multi.jl script for AA-alignment masking', metavar='FILE')
    group_io.add_argument('--gene-julia', default='julia', help='Julia executable used with --gene-taper', metavar='FILE')
    group_io.add_argument('--gene-species-mode', choices=('strict', 'multicopy'), default='strict', help='gene-tree route: strict uses resolved SC-OGs; multicopy uses all family trees', metavar='MODE')
    group_io.add_argument('--gene-aster', default='astral', help='ASTER2 astral executable for both gene-tree routes', metavar='FILE')

    group_filter = parser.add_argument_group('arguments for filtering')
    group_filter.add_argument('-kf', default=31, help='Filter k-mer size', metavar='INT', type=int)
    group_filter.add_argument('-s', '--step-size', default=4, help='Filter step size', metavar='INT', type=int)
    group_filter.add_argument('--max-reads', default=0, help='Million reads to process per file', metavar='INT', type=int)
    group_filter.add_argument('--reuse-reference-cache', action='store_true', default=False, help='Reuse a fingerprinted reference k-mer cache instead of rebuilding it every run')
    group_filter.add_argument('--reference-cache-dir', default=None, help='Directory for --reuse-reference-cache files (default = output/.gm2_reference_cache)', metavar='DIR')

    group_refilter = parser.add_argument_group('arguments for futher filtering')
    group_refilter.add_argument('--depth-low-water-mark', default=50, help='If depth is lower than this value, try to find more reads with relaxed criteria', metavar='INT', type=int)
    group_refilter.add_argument('--depth-limit', default=768, help='Maximum depth processed during re-filtering', metavar='INT', type=int)
    group_refilter.add_argument('--file-size-limit', default=6, help='Maximum file size during re-filtering', metavar='INT', type=int)

    group_assembly = parser.add_argument_group('arguments for assembly')
    group_assembly.add_argument('-ka', default=0, help='Assembly k-mer size (default = auto)', metavar='INT', type=int)
    group_assembly.add_argument('--min-ka', default=21, help='Minimum auto-estimated assembly k-mer size', metavar='INT', type=int)
    group_assembly.add_argument('--max-ka', default=51, help='Maximum auto-estimated assembly k-mer size', metavar='INT', type=int)
    group_assembly.add_argument('-e', '--error-threshold', default=2, help='Error threshold', metavar='INT', type=int)
    group_assembly.add_argument('-sb', '--soft-boundary', default='auto', help='Soft boundary (default = auto)', metavar='{INT,auto,unlimited}', type=str)
    group_assembly.add_argument('-i', '--search-depth', default=4096, help='Search depth', metavar='INT', type=int)
    group_assembly.add_argument('--min-coverage', default=0, help='Minimum read depth required for contig generation', metavar='INT', type=int)
    group_assembly.add_argument('--assembler-implementation', choices=('auto', 'uce-rust', 'original', 'original-rust'), default='auto', help='Assembler implementation: auto uses original-rust in original mode and uce-rust in UCE mode; original-rust is the deterministic Rust compatibility implementation for original mode; uce-rust selects the UCE-oriented Rust assembler; original and original-rust are unavailable in UCE mode')
    group_assembly.add_argument('--assembler-read-chunk-size', default=8192, help='Reads per bounded Rust assembler batch (default = 8192)', metavar='INT', type=int)
    group_assembly.add_argument('--assembly-mode', choices=('original', 'uce'), default='original', help='Assembly mode: original performs conventional reference-guided recovery; uce preserves UCE flanks')
    group_assembly.add_argument('--assembler-graph-format', choices=('none', 'gfa', 'dot', 'both'), default='none', help='Write compact per-locus Rust assembly graphs (default = none)')
    group_profile = parser.add_argument_group('arguments for marker profiling')
    group_profile.add_argument('--profile-kmer-size', default=21, help='Profiling: k-mer size for both recruitment and Themisto (odd integer 15-31; default = 21)', metavar='INT', type=int)
    group_profile.add_argument('--profile-pseudoalign-threshold', default=0.80, help='Profiling: Themisto pseudoalignment threshold (default = 0.80)', metavar='FLOAT', type=float)
    group_profile.add_argument('--profile-relevant-kmer-fraction', default=0.50, help='Profiling: minimum fraction of query k-mers found in the reference index (default = 0.50)', metavar='FLOAT', type=float)
    group_profile.add_argument('--profile-group-map', default='', help='Profiling: optional TSV mapping reference ID to reporting group', metavar='FILE')
    group_profile.add_argument('--profile-decoy', default='', help='Profiling: optional non-target decoy FASTA', metavar='FILE')
    group_profile.add_argument('--profile-index-dir', default='', help='Profiling: cache directory for split references and the Themisto index', metavar='DIR')
    group_profile.add_argument('--profile-index-memory-gb', default=2, help='Profiling: Themisto index-build memory limit in GiB (default = 2)', metavar='INT', type=int)
    group_profile.add_argument('--profile-themisto', default='', help='Profiling: Themisto executable path; by default use PATH', metavar='FILE')
    group_profile.add_argument('--profile-force-rebuild', action='store_true', default=False, help='Profiling: rebuild the cached Themisto reference index')
    group_assembly.add_argument('--uce-side-candidates', default=8, help='One-sided branch candidates to combine in UCE mode (default = 8)', metavar='INT', type=int)
    group_assembly.add_argument('--uce-path-strategy', choices=('search', 'backbone'), default='backbone', help='UCE path handling: backbone commits one bounded-lookahead path without backtracking; search preserves legacy branch enumeration (default = backbone)')
    group_assembly.add_argument('--uce-backbone-lookahead', default=24, help='Greedy look-ahead steps used to choose a UCE backbone edge at each bubble (default = 24)', metavar='INT', type=int)
    group_assembly.add_argument('--uce-max-contig-length', default=0, help='Maximum UCE contig length kept before scoring; use 0 to disable (default = 0)', metavar='INT', type=int)
    group_assembly.add_argument('--uce-min-read-density', default=0.003, help='Minimum uniquely placed read_count/length for long UCE contigs before scoring (default = 0.003)', metavar='FLOAT', type=float)
    group_assembly.add_argument('--uce-density-check-min-length', default=1000, help='Minimum contig length where the UCE read-density guardrail applies (default = 1000)', metavar='INT', type=int)
    group_assembly.add_argument('--uce-max-depth-cv', default=0, help='Optional maximum k-mer depth coefficient of variation for UCE contigs; 0 disables (default = 0)', metavar='FLOAT', type=float)
    group_assembly.add_argument('--uce-max-depth-ratio', default=0, help='Optional maximum max/median k-mer depth ratio for UCE contigs; 0 disables (default = 0)', metavar='FLOAT', type=float)
    group_assembly.add_argument('--uce-rescue-reads', action='store_true', default=False, help='UCE mode only: after the first assembly, recruit raw reads once using preliminary contigs plus original references, then re-filter and re-assemble')
    group_assembly.add_argument('--uce-rescue-min-contig-length', default=60, help='Minimum preliminary contig length used as a UCE raw-read rescue reference (default = 60)', metavar='INT', type=int)
    group_assembly.add_argument('--uce-rescue-min-density-ratio', default=0.5, help='Minimum rescue/before read-density ratio kept after UCE raw-read rescue (default = 0.5)', metavar='FLOAT', type=float)

    group_te = parser.add_argument_group('arguments for reference-free repeatome analysis')
    group_te.add_argument('--te-stage', choices=('all', 'discover', 'curate', 'annotate', 'quantify'), default='all', help='Repeatome stage: discover, curate, annotate, or quantify (default = all)')
    group_te.add_argument('--te-kmer', default=25, type=int, help='Canonical k-mer length for repeat seeds (default = 25)', metavar='INT')
    group_te.add_argument('--te-min-kmer-count', default=8, type=int, help='Minimum sampled k-mer count retained in the repeat catalog (default = 8)', metavar='INT')
    group_te.add_argument('--te-catalog-pairs', default=10000, type=int, help='Read-pair quota per taxon for catalog discovery (default = 10000)', metavar='INT')
    group_te.add_argument('--te-read-ledger', default='', help='Optional TSV of high-confidence read IDs to exclude: sample_id read_id', metavar='FILE')
    group_te.add_argument('--te-library', default='', help='Optional classified TE-library FASTA used during the annotate stage (headers: name#Class/Subclass)', metavar='FILE')
    group_te.add_argument('--te-annotate-min-fragment', default=80, type=int, help='Minimum read-supported fragment length for homology annotation (default = 80)', metavar='INT')
    group_te.add_argument('--te-annotate-max-fragment', default=800, type=int, help='Maximum bounded annotation-fragment length (default = 800)', metavar='INT')
    group_te.add_argument('--te-annotate-min-support', default=5, type=int, help='Minimum unique read-pair support for a confident annotation (default = 5)', metavar='INT')
    group_te.add_argument('--te-annotate-min-identity', default=0.80, type=float, help='Minimum gapless identity for a confident TE-library annotation (default = 0.80)', metavar='FLOAT')
    group_te.add_argument('--te-annotate-min-coverage', default=0.60, type=float, help='Minimum fragment coverage for a confident TE-library annotation (default = 0.60)', metavar='FLOAT')
    group_te.add_argument('--te-annotate-min-delta', default=0.10, type=float, help='Minimum best-versus-second-class score margin (default = 0.10)', metavar='FLOAT')

    group_population = parser.add_argument_group('arguments for population SNP analysis')
    group_population.add_argument('--engine', choices=('pseudoref', 'panref', 'panrefv2'), default='pseudoref', help='Population reference engine: accepted-contig pseudoref, legacy panref, or streaming two-pass PanRefV2')
    group_population.add_argument('--population-panrefv2-include-low-confidence', action='store_true', default=False, help='Include short or low-support PanRefV2 loci in the mapping FASTA')
    group_population.add_argument('--population-reference-strategy', choices=('sqcl-longest', 'supported'), default='sqcl-longest', help='Public-reference representative selection: SqCL-like longest accepted contig or support-first (default = sqcl-longest)')
    group_population.add_argument('--population-reference-fasta', default=None, help='Use a fixed external cohort FASTA instead of building a reference from accepted contigs', metavar='FILE')
    group_population.add_argument('--population-min-mapq', default=20, help='Minimum mapping quality for joint calling (default = 20)', metavar='INT', type=int)
    group_population.add_argument('--population-min-baseq', default=20, help='Minimum base quality for joint calling (default = 20)', metavar='INT', type=int)
    group_population.add_argument('--population-min-dp', default=5, help='Set genotypes below this depth to missing (default = 5)', metavar='INT', type=int)
    group_population.add_argument('--population-min-gq', default=20, help='Set genotypes below this quality to missing (default = 20)', metavar='INT', type=int)
    group_population.add_argument('--population-min-qual', default=20.0, help='Minimum site QUAL (default = 20)', metavar='FLOAT', type=float)
    group_population.add_argument('--population-min-call-rate', default=0.8, help='Minimum non-missing genotype fraction (default = 0.8)', metavar='FLOAT', type=float)
    group_population.add_argument('--population-min-mac', default=2, help='Minimum minor allele count (default = 2)', metavar='INT', type=int)
    group_population.add_argument('--population-ld-window', default=50, help='SNPs per LD-pruning window (default = 50)', metavar='INT', type=int)
    group_population.add_argument('--population-ld-step', default=5, help='SNPs shifted per LD-pruning window (default = 5)', metavar='INT', type=int)
    group_population.add_argument('--population-ld-r2', default=0.2, help='LD-pruning r-squared threshold (default = 0.2)', metavar='FLOAT', type=float)
    group_population.add_argument('--population-admixture-k-min', default=2, help='Minimum ADMIXTURE K (default = 2)', metavar='INT', type=int)
    group_population.add_argument('--population-admixture-k-max', default=6, help='Maximum ADMIXTURE K (default = 6)', metavar='INT', type=int)
    group_population.add_argument('--population-admixture-cv', default=10, help='ADMIXTURE cross-validation folds (default = 10)', metavar='INT', type=int)
    group_population.add_argument('--population-start-at', choices=('reference', 'mapping', 'calling', 'selection'), default='reference', help='Start at this population stage, reusing validated existing outputs when later than reference (default = reference)')
    group_population.add_argument('--population-stop-after', choices=('reference', 'mapping', 'calling', 'selection'), default='selection', help='Stop after this population stage (default = selection)')
    group_population.add_argument('--population-skip-mark-duplicates', action='store_true', default=False, help='Skip samtools duplicate marking')
    group_population.add_argument('--population-skip-plink', action='store_true', default=False, help='Do not export PLINK files or run PCA, LD pruning or ADMIXTURE')
    group_population.add_argument('--population-skip-admixture', action='store_true', default=False, help='Do not run ADMIXTURE on the primary one-SNP-per-UCE panel')
    group_population.add_argument('--population-minibwa', default='minibwa', help='minibwa executable (default = minibwa)', metavar='PATH')
    group_population.add_argument('--population-samtools', default='samtools', help='samtools executable (default = samtools)', metavar='PATH')
    group_population.add_argument('--population-bcftools', default='bcftools', help='bcftools executable (default = bcftools)', metavar='PATH')
    group_population.add_argument('--population-plink', default='plink', help='PLINK 1.9 executable (default = plink)', metavar='PATH')
    group_population.add_argument('--population-admixture', default='admixture', help='ADMIXTURE executable (default = admixture)', metavar='PATH')

    group_consensus = parser.add_argument_group('argument for consensus generation')
    group_consensus.add_argument('-c', '--consensus-threshold', default='0.75', help='Consensus threshold (default = 0.75)', metavar='FLOAT', type=float)

    group_trim = parser.add_argument_group('arguments for sequence trimming')
    group_trim.add_argument('-ts', '--trim-source', choices=('assembly', 'consensus'), default=None, help='Whether to trim the primary assembly or the consensus sequence (default = output of last step, assembly if no other command given)')
    group_trim.add_argument('-tm', '--trim-mode', choices=('all', 'longest', 'terminal', 'isoform'), default='terminal', help='Trim mode (default = terminal)', type=str)
    group_trim.add_argument('-tr', '--trim-retention', default=0, help='Retention length threshold (default = 0.0)', metavar='FLOAT', type=float)

    group_combine = parser.add_argument_group('arguments for sequence alignment and clustering')
    group_combine.add_argument('-cs', '--combine-source', choices=('assembly', 'consensus', 'trimmed'), default=None, help='Whether to combine the primary assembly, the consensus sequences or the trimmed sequences (default = output of last step, assembly if no other command given)')
    group_combine.add_argument('-cd', '--clean-difference', default=1, help='Maximum acceptable pairwise difference in an alignment (default = 1.0)', metavar='FLOAT', type=float)
    group_combine.add_argument('-cn', '--clean-sequences', default=0, help='Number of sequences required in an alignment (default = 0)', metavar='INT', type=int)
    group_combine.add_argument('--msa-program', choices=('clustalo', 'mafft'), default='mafft', help='Program for multiple sequence alignment', type=str)
    group_combine.add_argument('--msa-threads', default=1, help='Threads used by each multiple-sequence-alignment job (default = 1)', metavar='INT', type=int)
    group_combine.add_argument('--alignment-filter', choices=('trimal', 'alifilter', 'none'), default=None, help='Program for filtering aligned loci before tree reconstruction (default = trimal)', type=str)
    group_combine.add_argument('--filter-processes', default=None, help='Maximum number of concurrent alignment filtering jobs (default = -p)', metavar='INT', type=int)
    group_combine.add_argument('--alifilter-model', default=None, help='AliFilter model specification or model.json path when --alignment-filter alifilter is used', metavar='MODEL', type=str)
    group_combine.add_argument('--strict-combine-errors', action='store_true', default=False, help='Stop combine if any locus fails during alignment, cleanup, or alignment filtering')
    group_combine.add_argument('--no-alignment', action='store_true', default=False, help='Do not perform multiple sequence alignment')
    group_combine.add_argument('--no-trimal', action='store_true', default=False, help='Do not run alignment filtering (deprecated; use --alignment-filter none)')

    group_tree = parser.add_argument_group('argument for tree inference')
    group_tree.add_argument('-m', '--tree-method', choices=('coalescent', 'concatenation'), default='coalescent', help='Multi-gene tree reconstruction method (default = coalescent)')
    group_tree.add_argument('-b', '--bootstrap', default=1000, help='Number of bootstrap replicates', metavar='INT', type=int)
    group_tree.add_argument('--phylo-program', choices=('raxmlng', 'iqtree', 'fasttree', 'veryfasttree'), default='fasttree', help='Program for phylogenetic tree reconstruction', type=str)

    group_stats = parser.add_argument_group('arguments for UCE statistics')
    group_stats.add_argument('--stats-no-heatmap', action='store_true', default=False, help='Do not create UCE statistics heatmaps')
    group_stats.add_argument('--stats-count-input-reads', action='store_true', default=False, help='Count input FASTQ reads for InputReads and PctFiltered statistics; can be slow for large datasets')

    parser.add_argument('--min-depth', help=argparse.SUPPRESS)
    parser.add_argument('--max-depth', help=argparse.SUPPRESS)

    args = parser.parse_args()
    args.is_mito_workflow = False
    args.is_gene_workflow = False

    if 'gene-annotate' in args.command:
        if len(args.command) != 1 or not args.gene_input or not args.gene_protein_reference:
            parser.error('gene-annotate requires --gene-input and --gene-protein-reference only')
    elif 'gene-resolve' in args.command:
        if len(args.command) != 1:
            parser.error('gene-resolve is a complete workflow and cannot be combined with other subcommands')
        if not args.gene_input:
            parser.error('gene-resolve requires --gene-input')
        if args.gene_min_taxa < 2:
            parser.error('--gene-min-taxa must be at least 2')
        if args.gene_ufboot != 0 and args.gene_ufboot < 1000:
            parser.error('--gene-ufboot must be 0 or at least 1000')
    elif 'gene-tree' in args.command:
        if len(args.command) != 1:
            parser.error('gene-tree is a complete workflow and cannot be combined with other subcommands')
        if not args.gene_input:
            parser.error('gene-tree requires --gene-input')
    elif 'gene' in args.command:
        if len(args.command) != 1:
            parser.error('gene is a complete workflow and cannot be combined with other subcommands')
        if args.assembly_mode != 'original':
            parser.error('gene requires --assembly-mode original')
        if args.assembler_implementation not in ('auto', 'original-rust'):
            parser.error('gene requires the original-rust assembler')
        args.is_gene_workflow = True
        args.command = ('filter', 'refilter', 'assemble')
        args.assembler_implementation = 'original-rust'
    elif 'mito' in args.command:
        args.is_mito_workflow = True
        if len(args.command) != 1:
            parser.error("mito is a complete workflow and cannot be combined with other subcommands")
        if not args.mito_genbank:
            parser.error("mito requires --mito-genbank")
        if args.mito_tile_length < 1 or args.mito_tile_step < 1 or args.mito_tile_step > args.mito_tile_length:
            parser.error("mito requires 0 < --mito-tile-step <= --mito-tile-length")
        if args.mito_min_overlap < 1 or not 0 < args.mito_min_overlap_identity <= 1:
            parser.error("mito overlap parameters require positive overlap and identity in (0, 1]")
        if args.mito_min_junction_support < 1:
            parser.error("--mito-min-junction-support must be positive")
        if args.mito_terminal_window < 1 or args.mito_min_link_hits < 1 or args.mito_min_pair_support < 1:
            parser.error("mito terminal-window, link-hit and pair-support values must be positive")
        if not 1 <= args.mito_link_kmer <= 63 or not 1 <= args.mito_bridge_kmer <= 63:
            parser.error("mito link and bridge k-mers must be between 1 and 63")
        if args.mito_bridge_min_depth < 1 or args.mito_max_bridge < 1:
            parser.error("mito bridge depth and maximum bridge must be positive")
        if args.mito_initial_reads < 1:
            parser.error('--mito-initial-reads must be positive')
        if args.mito_max_reads < args.mito_initial_reads:
            parser.error('--mito-max-reads must be at least --mito-initial-reads')
        args.command = ('filter', 'refilter', 'assemble')
        args.assembly_mode = 'uce'
        if args.ka == 0:
            args.ka = 31
        args.uce_min_read_density = 0
        args.search_depth = max(args.search_depth, 30000)
    elif args.mito_genbank:
        parser.error('--mito-genbank is only valid with the mito subcommand')
    elif not args.r and 'te' not in args.command and 'gene-tree' not in args.command and 'gene-resolve' not in args.command and 'gene-annotate' not in args.command:
        parser.error("-r is required unless the mito or te subcommand is used")

    if not args.f and 'gene-tree' not in args.command and 'gene-resolve' not in args.command and 'gene-annotate' not in args.command:
        parser.error('-f is required unless gene-resolve is used')

    if args.reference_cache_dir and not args.reuse_reference_cache:
        parser.error('--reference-cache-dir requires --reuse-reference-cache')

    args.uce_side_candidates = max(args.uce_side_candidates, 3)
    if 'profiling' in args.command:
        if args.profile_group_map and not os.path.isfile(args.profile_group_map):
            parser.error('--profile-group-map must be a readable TSV file')
        if args.profile_kmer_size < 15 or args.profile_kmer_size > 31 or args.profile_kmer_size % 2 == 0:
            parser.error('--profile-kmer-size must be an odd integer from 15 to 31')
        args.kf = args.ka = args.min_ka = args.max_ka = args.profile_kmer_size
        if not 0 < args.profile_pseudoalign_threshold <= 1:
            parser.error('--profile-pseudoalign-threshold must be in (0, 1]')
        if not 0 <= args.profile_relevant_kmer_fraction <= 1:
            parser.error('--profile-relevant-kmer-fraction must be in [0, 1]')
        if args.profile_index_memory_gb < 1:
            parser.error('--profile-index-memory-gb must be at least 1')
    args.uce_backbone_lookahead = max(args.uce_backbone_lookahead, 1)
    args.uce_max_contig_length = max(args.uce_max_contig_length, 0)
    args.uce_density_check_min_length = max(args.uce_density_check_min_length, 1)
    args.uce_rescue_min_contig_length = max(args.uce_rescue_min_contig_length, args.kf)

    if args.uce_min_read_density < 0:
        parser.error('--uce-min-read-density must be greater than or equal to 0')

    if args.uce_max_depth_cv < 0:
        parser.error('--uce-max-depth-cv must be greater than or equal to 0')

    if args.uce_max_depth_ratio < 0:
        parser.error('--uce-max-depth-ratio must be greater than or equal to 0')

    if args.uce_rescue_min_density_ratio <= 0:
        parser.error('--uce-rescue-min-density-ratio must be greater than 0')

    if args.population_min_mapq < 0 or args.population_min_baseq < 0:
        parser.error('population mapping and base quality thresholds must be non-negative')

    if args.population_min_dp < 0 or args.population_min_gq < 0 or args.population_min_mac < 1:
        parser.error('population DP/GQ thresholds must be non-negative and MAC must be at least 1')

    if not 0 < args.population_min_call_rate <= 1:
        parser.error('--population-min-call-rate must be in (0, 1]')

    if args.population_ld_window < 1 or args.population_ld_step < 1 or not 0 < args.population_ld_r2 < 1:
        parser.error('population LD window/step must be positive and LD r-squared must be in (0, 1)')

    if args.population_admixture_k_min < 2 or args.population_admixture_k_max < args.population_admixture_k_min or args.population_admixture_cv < 2:
        parser.error('population ADMIXTURE requires 2 <= K-min <= K-max and at least 2 CV folds')

    if args.no_trimal and args.alignment_filter not in (None, 'none'):
        parser.error('--no-trimal cannot be combined with --alignment-filter trimal or --alignment-filter alifilter')

    if get_alifilter_model(args) and get_alignment_filter(args) != 'alifilter':
        parser.error('--alifilter-model requires --alignment-filter alifilter')

    if args.p < 1:
        parser.error('-p must be at least 1')

    if args.msa_threads < 1:
        parser.error('--msa-threads must be at least 1')

    if args.msa_threads > args.p:
        parser.error('--msa-threads cannot be greater than -p')

    if args.filter_processes is not None and args.filter_processes < 1:
        parser.error('--filter-processes must be at least 1')

    if not args.command:
        if args.assembly_mode == 'uce':
            args.command = ('filter', 'refilter', 'assemble', 'combine', 'tree')
        else:
            args.command = ('filter', 'refilter', 'assemble', 'trim', 'combine', 'tree')

    if args.min_depth is not None:
        print('  Option --min-depth has been removed. Please use --depth-low-water-mark, --error-threshold or --min-coverage instead.')
        print(DEPTH_DEPRECATION_EXPLAINER)
        sys.exit(2)

    if args.max_depth is not None:
        print('  Option --max-depth has been removed. Please use --depth-limit instead.')
        print(DEPTH_DEPRECATION_EXPLAINER)
        sys.exit(2)

    try:
        if args.is_mito_workflow:
            prepare_mito_reference(args)
        materialize_profile_reference(args)
        prepare_profile_cache_key(args)
    except RuntimeError as e:
        parser.error(str(e))

    if 'gene-annotate' in args.command or 'gene-resolve' in args.command or 'gene-tree' in args.command:
        print(f'Running tasks: {args.command[0]}')
        sys.exit(execute_tasks(args, {}))

    samples = prepare_workdir(args)
    if samples:
        print(f'Running tasks: {", ".join(args.command)}')
        print()
        sys.exit(execute_tasks(args, samples))
    else:
        print('Sample list is empty or invalid, exiting')
        sys.exit(2)
