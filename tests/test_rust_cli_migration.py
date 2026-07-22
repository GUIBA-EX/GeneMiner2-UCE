from pathlib import Path
import csv
import os
import random
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "rust" / "geneminer2_cli" / "Cargo.toml"


class RustCliMigrationTests(unittest.TestCase):
    def test_rust_is_the_default_engine_without_environment_override(self):
        proc = subprocess.run(
            ["cargo", "run", "--quiet", "--manifest-path", str(MANIFEST), "--", "--help"],
            cwd=ROOT,
            text=True,
            env={key: value for key, value in os.environ.items() if key != "GENEMINER2_ENGINE"},
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("no Python runtime is required", proc.stdout)

    def test_output_comparator_accepts_identical_directories(self):
        fixture = ROOT / "tests" / "__pycache__"
        proc = subprocess.run(
            ["cargo", "run", "--quiet", "--manifest-path", str(MANIFEST), "--bin", "compare_outputs", "--", str(fixture), str(fixture)],
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Compatibility check passed", proc.stdout)

    def _fake_component(self, directory, name):
        target = directory / name
        target.write_text('#!/bin/sh\nprintf "%s\n" "$@" > "$GM2_CAPTURE"\n')
        target.chmod(0o755)
        return target

    def test_native_te_forwards_legacy_defaults_to_rust_backend(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            components, capture = tmp / 'components', tmp / 'te.args'
            components.mkdir()
            self._fake_component(components, 'main_repeat')
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--',
                 'te', '-f', 'taxa.tsv', '-o', str(tmp / 'out'), '-p', '3',
                 '--te-stage', 'discover', '--te-library', 'library.fasta'],
                cwd=ROOT, text=True,
                env={**os.environ, 'GENEMINER2_ENGINE': 'rust',
                     'GM2_COMPONENT_DIR': str(components), 'GM2_CAPTURE': str(capture)},
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            received = capture.read_text().splitlines()
            self.assertEqual(received[:6], ['--samples', 'taxa.tsv', '--output', str(tmp / 'out'), '--stage', 'discover'])
            self.assertIn(str(components / 'MainFilterNew'), received)
            self.assertEqual(received[received.index('--threads') + 1], '3')
            self.assertEqual(received[received.index('--te-library') + 1], 'library.fasta')

    def test_native_population_forwards_panref_options_to_rust_backend(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            components, capture = tmp / 'components', tmp / 'population.args'
            components.mkdir()
            self._fake_component(components, 'main_population')
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--',
                 'population', '-f', 'samples.tsv', '-r', 'baits', '-o', str(tmp / 'out'),
                 '-p', '2', '--engine', 'panrefv2',
                 '--population-panrefv2-include-low-confidence',
                 '--population-skip-plink', '--population-stop-after', 'mapping'],
                cwd=ROOT, text=True,
                env={**os.environ, 'GENEMINER2_ENGINE': 'rust',
                     'GM2_COMPONENT_DIR': str(components), 'GM2_CAPTURE': str(capture)},
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            received = capture.read_text().splitlines()
            self.assertEqual(received[received.index('--engine') + 1], 'panrefv2')
            self.assertEqual(received[received.index('--panref-baits') + 1], 'baits')
            self.assertEqual(received[received.index('--threads') + 1], '2')
            self.assertEqual(received[received.index('--stop-after') + 1], 'mapping')
            self.assertIn('--panrefv2-include-low-confidence', received)
            self.assertIn('--skip-plink', received)

    def test_native_stats_runs_rust_backend(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            references, output = tmp / 'references', tmp / 'out'
            references.mkdir()
            (references / 'uce_demo.fasta').write_text('>uce_demo\nACGTACGT\n')
            reads = tmp / 'reads.fq'
            reads.write_text('@read\nACGTACGT\n+\nIIIIIIII\n')
            sample_list = tmp / 'samples.tsv'
            sample_list.write_text(f'demo\t{reads}\t{reads}\n')
            (output / '1_Demo').mkdir(parents=True)
            (output / '1_Demo' / 'uce_assembly_summary.csv').write_text(
                'locus,accepted,selected_contig_length\nuce_demo,1,8\n'
            )
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--',
                 'stats', '-f', str(sample_list), '-r', str(references), '-o', str(output),
                 '--stats-no-heatmap'],
                cwd=ROOT, text=True,
                env={**os.environ, 'GENEMINER2_ENGINE': 'rust',
                     'GM2_COMPONENT_DIR': str(ROOT / 'cli' / 'bin')},
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue((output / 'uce_stats.tsv').is_file())

    def test_native_consensus_preserves_legacy_task_layout(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            components = tmp / 'components'
            components.mkdir()
            minimap = components / 'minimap2'
            minimap.write_text(
                '#!/bin/sh\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-o" ]; then out="$2"; shift 2; else shift; fi; done\n'
                'printf "@HD\tVN:1.6\n" > "$out"\n'
            )
            minimap.chmod(0o755)
            consensus = components / 'build_consensus'
            consensus.write_text(
                '#!/bin/sh\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-o" ]; then out="$2"; shift 2; else shift; fi; done\n'
                'printf ">consensus\nACGT\n" > "$out/consensus.fasta"\n'
            )
            consensus.chmod(0o755)
            references, output = tmp / 'references', tmp / 'out'
            references.mkdir()
            (references / 'gene.fasta').write_text('>gene\nACGT\n')
            reads = tmp / 'reads.fastq'
            reads.write_text('@read\nACGT\n+\nIIII\n')
            sample_list = tmp / 'samples.tsv'
            sample_list.write_text(f'demo\t{reads}\n')
            sample = output / '1_Demo'
            (sample / 'results').mkdir(parents=True)
            (sample / 'filtered').mkdir()
            (sample / 'results' / 'gene.fasta').write_text('>gene\nACGT\n')
            (sample / 'filtered' / 'gene.fq').write_text(reads.read_text())
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--',
                 'consensus', '-f', str(sample_list), '-r', str(references), '-o', str(output), '-p', '2'],
                cwd=ROOT, text=True,
                env={**os.environ, 'GENEMINER2_ENGINE': 'rust',
                     'GM2_COMPONENT_DIR': str(components),
                     'GM2_MINIMAP2': str(minimap)},
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue((sample / 'consensus' / 'consensus.fasta').is_file())
            self.assertFalse((sample / 'consensus' / 'gene.sam').exists())

    def test_native_trim_preserves_blast_database_and_task_layout(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            components = tmp / 'components'
            components.mkdir()
            makeblastdb = components / 'makeblastdb'
            makeblastdb.write_text(
                '#!/bin/sh\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-out" ]; then touch "$2.nhr"; shift 2; else shift; fi; done\n'
            )
            makeblastdb.chmod(0o755)
            trim = components / 'build_trimed'
            trim.write_text(
                '#!/bin/sh\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-o" ]; then out="$2"; shift 2; else shift; fi; done\n'
                'printf ">trimmed\nACGT\n" > "$out"\n'
            )
            trim.chmod(0o755)
            references, output = tmp / 'references', tmp / 'out'
            references.mkdir()
            (references / 'gene.fasta').write_text('>gene\nACGT\n')
            sample_list = tmp / 'samples.tsv'
            sample_list.write_text('demo\treads.fq\n')
            sample = output / '1_Demo'
            (sample / 'results').mkdir(parents=True)
            (sample / 'results' / 'gene.fasta').write_text('>gene\nACGT\n')
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--',
                 'trim', '-f', str(sample_list), '-r', str(references), '-o', str(output),
                 '--trim-mode', 'longest', '--trim-retention', '0.5'],
                cwd=ROOT, text=True,
                env={**os.environ, 'GENEMINER2_ENGINE': 'rust',
                     'GM2_COMPONENT_DIR': str(components),
                     'GM2_MAKEBLASTDB': str(makeblastdb), 'GM2_BLASTN': 'blastn'},
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue((output / 'blast_db' / 'gene.nhr').is_file())
            self.assertTrue((sample / 'blast' / 'gene.fasta').is_file())

    def test_native_combine_merges_first_records_and_honors_uce_acceptance(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            references, output = tmp / 'references', tmp / 'out'
            references.mkdir()
            (references / 'gene.fasta').write_text('>gene\nACGT\n')
            sample_list = tmp / 'samples.tsv'
            sample_list.write_text('one\tr1.fq\ntwo\tr2.fq\n')
            one, two = output / '1_One', output / '2_Two'
            for sample, sequence, accepted in ((one, 'AAAA', '1'), (two, 'CCCC', '0')):
                (sample / 'results').mkdir(parents=True)
                if sample == one:
                    (sample / 'results' / 'gene.fasta').write_text(f'>first\n{sequence}\n>second\nGGGG\n')
                (sample / 'uce_assembly_summary.csv').write_text(f'locus,accepted\ngene,{accepted}\n')
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--',
                 'combine', '--no-alignment', '--assembly-mode', 'uce',
                 '-f', str(sample_list), '-r', str(references), '-o', str(output)],
                cwd=ROOT, text=True,
                env={**os.environ, 'GENEMINER2_ENGINE': 'rust',
                     'GM2_COMPONENT_DIR': str(ROOT / 'cli' / 'bin')},
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            merged = (output / 'combined_results' / 'gene.fasta').read_text()
            self.assertEqual(merged, '>1_One\nAAAA\n')

    def test_native_tree_builds_coalescent_tree_from_locus_alignments(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            tools = tmp / 'tools'
            tools.mkdir()
            fasttree = tools / 'FastTree'
            fasttree.write_text(
                '#!/bin/sh\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-out" ]; then out="$2"; shift 2; else shift; fi; done\n'
                'printf "(A,B);\n" > "$out"\n'
            )
            fasttree.chmod(0o755)
            astral = tools / 'astral'
            astral.write_text(
                '#!/bin/sh\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-o" ]; then out="$2"; shift 2; else shift; fi; done\n'
                'printf "(A,B);\n" > "$out"\n'
            )
            astral.chmod(0o755)
            references, output = tmp / 'references', tmp / 'out'
            references.mkdir()
            (references / 'gene.fasta').write_text('>gene\nACGT\n')
            aligned = output / 'combined_results' / 'aligned'
            aligned.mkdir(parents=True)
            (aligned / 'gene.fasta').write_text('>A\nACGT\n>B\nACGT\n')
            samples = tmp / 'samples.tsv'
            samples.write_text('a\tr1.fq\n')
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--',
                 'tree', '-f', str(samples), '-r', str(references), '-o', str(output),
                 '--alignment-filter', 'none', '--phylo-program', 'fasttree'],
                cwd=ROOT, text=True,
                env={**os.environ, 'GENEMINER2_ENGINE': 'rust',
                     'GM2_COMPONENT_DIR': str(ROOT / 'cli' / 'bin'),
                     'GM2_FASTTREE': str(fasttree), 'GM2_ASTRAL': str(astral)},
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual((output / 'Coalescent.tree').read_text(), '(A,B);\n')
            self.assertEqual((output / 'combined_genes.trees').read_text(), '(A,B);\n')

    def test_native_gene_annotate_forwards_to_gene_workflow_without_samples(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            components, capture = tmp / 'components', tmp / 'args.txt'
            components.mkdir()
            workflow = components / 'gene_workflow'
            workflow.write_text('#!/bin/sh\nprintf "%s\n" "$@" > "$GM2_CAPTURE"\n')
            workflow.chmod(0o755)
            annotation, proteins = tmp / 'gene', tmp / 'proteins'
            annotation.mkdir(); proteins.mkdir()
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--',
                 'gene-annotate', '--gene-input', str(annotation),
                 '--gene-protein-reference', str(proteins), '-o', str(tmp / 'out'), '-p', '2'],
                cwd=ROOT, text=True,
                env={**os.environ, 'GENEMINER2_ENGINE': 'rust', 'GM2_COMPONENT_DIR': str(components),
                     'GM2_CAPTURE': str(capture)}, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            received = capture.read_text().splitlines()
            self.assertEqual(received[0], 'annotate')
            self.assertEqual(received[received.index('--threads') + 1], '2')

    def test_native_gene_resolve_forwards_all_resolution_options(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            components, capture = tmp / 'components', tmp / 'resolve.args'
            components.mkdir()
            workflow = components / 'gene_workflow'
            workflow.write_text('#!/bin/sh\nprintf "%s\n" "$@" > "$GM2_CAPTURE"\n')
            workflow.chmod(0o755)
            annotation, outgroup, taper = tmp / 'annotation', tmp / 'outgroups.tsv', tmp / 'taper.jl'
            annotation.mkdir(); outgroup.write_text('sample\n'); taper.write_text('# script\n')
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--',
                 'gene-resolve', '--gene-input', str(annotation), '-o', str(tmp / 'out'), '-p', '2',
                 '--gene-min-taxa', '5', '--gene-outgroup', str(outgroup), '--gene-ufboot', '1000',
                 '--gene-taper', str(taper), '--gene-julia', 'julia-bin'],
                cwd=ROOT, text=True,
                env={**os.environ, 'GENEMINER2_ENGINE': 'rust', 'GM2_COMPONENT_DIR': str(components), 'GM2_CAPTURE': str(capture)},
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            received = capture.read_text().splitlines()
            self.assertEqual(received[0], 'resolve')
            self.assertEqual(received[received.index('--threads') + 1], '2')
            self.assertEqual(received[received.index('--outgroup') + 1], str(outgroup))
            self.assertEqual(received[received.index('--ufboot') + 1], '1000')

    def test_native_combine_alignment_and_filter_chain(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); components = tmp / 'components'; components.mkdir()
            for name, script in {
                'fix_alignment': '#!/bin/sh\nexit 0\n',
                'merge_seq': '#!/bin/sh\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-output" ]; then out="$2"; shift 2; else shift; fi; done\nprintf ">merged\nACGT\n" > "$out"\n',
            }.items():
                path = components / name; path.write_text(script); path.chmod(0o755)
            mafft = components / 'mafft'; mafft.write_text('#!/bin/sh\nlast=""; for x in "$@"; do last="$x"; done; cat "$last"\n'); mafft.chmod(0o755)
            trimal = components / 'trimal'; trimal.write_text('#!/bin/sh\nwhile [ "$#" -gt 0 ]; do case "$1" in -in) inp="$2"; shift 2;; -out) out="$2"; shift 2;; *) shift;; esac; done; cp "$inp" "$out"\n'); trimal.chmod(0o755)
            refs, out = tmp / 'refs', tmp / 'out'; refs.mkdir(); (refs / 'gene.fasta').write_text('>gene\nACGT\n')
            samples = tmp / 'samples.tsv'; samples.write_text('a\tr1.fq\nb\tr2.fq\n')
            for number, seq in ((1, 'ACGT'), (2, 'AGGT')):
                target = out / f'{number}_{"A" if number == 1 else "B"}' / 'results'; target.mkdir(parents=True); (target / 'gene.fasta').write_text(f'>x\n{seq}\n')
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--', 'combine',
                 '-f', str(samples), '-r', str(refs), '-o', str(out), '-p', '2', '--msa-threads', '1', '--filter-processes', '1', '--alignment-filter', 'trimal'],
                cwd=ROOT, text=True,
                env={**os.environ, 'GENEMINER2_ENGINE': 'rust', 'GM2_COMPONENT_DIR': str(components), 'GM2_MAFFT': str(mafft), 'GM2_TRIMAL': str(trimal)},
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue((out / 'combined_results' / 'aligned' / 'gene.fasta').is_file())
            self.assertTrue((out / 'combined_trimed' / 'gene.fasta').is_file())
            self.assertTrue((out / 'combined_results.fasta').is_file())
            self.assertTrue((out / 'combined_trimed.fasta').is_file())

    def test_native_tree_builds_concatenation_tree(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); tools = tmp / 'tools'; tools.mkdir()
            fasttree = tools / 'FastTree'; fasttree.write_text('#!/bin/sh\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-out" ]; then out="$2"; shift 2; else shift; fi; done\nprintf "(A,B);\n" > "$out"\n'); fasttree.chmod(0o755)
            refs, out = tmp / 'refs', tmp / 'out'; refs.mkdir(); (refs / 'gene.fasta').write_text('>gene\nACGT\n')
            (out).mkdir(); (out / 'combined_results.fasta').write_text('>A\nACGT\n>B\nACGT\n')
            samples = tmp / 'samples.tsv'; samples.write_text('a\tr1.fq\n')
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--', 'tree',
                 '-f', str(samples), '-r', str(refs), '-o', str(out), '--tree-method', 'concatenation', '--alignment-filter', 'none', '--phylo-program', 'fasttree'],
                cwd=ROOT, text=True, env={**os.environ, 'GENEMINER2_ENGINE': 'rust', 'GM2_FASTTREE': str(fasttree)},
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual((out / 'Concatenation.tree').read_text(), '(A,B);\n')

    def test_native_gene_tree_writes_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            resolved = tmp / 'resolved'
            trees = resolved / 'astral_input'
            trees.mkdir(parents=True)
            (trees / 'resolved_1to1.trees').write_text('(A,B);\n')
            aster = tmp / 'astral'
            aster.write_text(
                '#!/bin/sh\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-o" ]; then out="$2"; shift 2; else shift; fi; done\n'
                'printf "(A,B);\n" > "$out"\n'
            )
            aster.chmod(0o755)
            output = tmp / 'tree'
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--',
                 'gene-tree', '--gene-input', str(resolved), '-o', str(output), '--gene-aster', str(aster)],
                cwd=ROOT, text=True, env={**os.environ, 'GENEMINER2_ENGINE': 'rust'},
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual((output / 'gene_strict_aster.tree').read_text(), '(A,B);\n')
            provenance = (output / 'gene_tree_provenance.tsv').read_text()
            self.assertIn('gene_trees_sha256', provenance)
            self.assertIn('species_tree_sha256', provenance)

    def test_native_profiling_recruits_then_quantifies_single_marker(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            components = tmp / 'components'
            components.mkdir()
            uce_filter = components / 'uce_filter'
            uce_filter.write_text(
                '#!/bin/sh\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-o" ]; then out="$2"; shift 2; else shift; fi; done\n'
                'mkdir -p "$out/filtered"\nprintf "@read\nACGT\n+\nIIII\n" > "$out/filtered/marker.fq"\n'
            )
            uce_filter.chmod(0o755)
            marker_profile = components / 'marker_profile'
            marker_profile.write_text(
                '#!/bin/sh\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "--output" ]; then out="$2"; shift 2; else shift; fi; done\n'
                'mkdir -p "$out"\nprintf "reference\tsupport\nmarker\t1\n" > "$out/marker_reference_support.tsv"\n'
            )
            marker_profile.chmod(0o755)
            reference = tmp / 'marker.fasta'
            reference.write_text('>marker\nACGT\n')
            reads = tmp / 'reads.fq'
            reads.write_text('@read\nACGT\n+\nIIII\n')
            samples = tmp / 'samples.tsv'
            samples.write_text(f'demo\t{reads}\n')
            output = tmp / 'out'
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--',
                 'profiling', '-f', str(samples), '-r', str(reference), '-o', str(output),
                 '--profile-themisto', 'themisto'],
                cwd=ROOT, text=True,
                env={**os.environ, 'GENEMINER2_ENGINE': 'rust', 'GM2_COMPONENT_DIR': str(components)},
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue((output / '1_Demo' / 'marker_profile' / 'marker_reference_support.tsv').is_file())
            self.assertTrue((output / '.marker_profile_reference' / 'marker.fasta').is_file())

    def test_native_mito_runs_bait_collapse_refilter_assembly_and_finalize(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            components = tmp / 'components'; components.mkdir()
            mainfilter = components / 'MainFilterNew'
            mainfilter.write_text(
                '#!/bin/sh\nwhile [ "$#" -gt 0 ]; do case "$1" in -o) out="$2"; shift 2;; -lkd) dict="$2"; shift 2;; -m) mode="$2"; shift 2;; *) shift;; esac; done\n'
                'mkdir -p "$out"\nif [ "$mode" = "2" ]; then : > "$dict"; else mkdir -p "$out/filtered_pe"; printf "@r/1\nACGT\n+\nIIII\n" > "$out/filtered_pe/mitochondrion_1.fq"; printf "@r/2\nACGT\n+\nIIII\n" > "$out/filtered_pe/mitochondrion_2.fq"; fi\n'
            ); mainfilter.chmod(0o755)
            refilter = components / 'main_refilter_new'
            refilter.write_text(
                '#!/bin/sh\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-o" ]; then out="$2"; shift 2; else shift; fi; done\n'
                'mkdir -p "$out"\nprintf "@r\nACGT\n+\nIIII\n" > "$out/mitochondrion.fq"\n'
            ); refilter.chmod(0o755)
            assembler = components / 'main_assembler-rust'
            assembler.write_text(
                '#!/bin/sh\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-o" ]; then out="$2"; shift 2; else shift; fi; done\n'
                'mkdir -p "$out/contigs_all" "$out/assembly_graphs"\nprintf ">mito\nACGTACGT\n" > "$out/contigs_all/mitochondrion.fasta"\nprintf "H\tVN:Z:1.0\n" > "$out/assembly_graphs/mitochondrion.gfa"\n'
            ); assembler.chmod(0o755)
            mito = components / 'mito_workflow'
            mito.write_text(
                '#!/bin/sh\ncmd="$1"; shift\ncase "$cmd" in\n'
                'prepare-reference) while [ "$#" -gt 0 ]; do if [ "$1" = "--out-dir" ]; then out="$2"; shift 2; else shift; fi; done; mkdir -p "$out/metadata"; printf ">ref\nACGT\n" > "$out/metadata/mitochondrial_reference.fasta"; printf "gene\n" > "$out/metadata/mitochondrial_genes.tsv"; printf ">bait\nACGT\n" > "$out/mitochondrion.fasta";;\n'
                'collapse-baits) while [ "$#" -gt 0 ]; do if [ "$1" = "--out-dir" ]; then out="$2"; shift 2; else shift; fi; done; mkdir -p "$out"; printf "@r/1\nACGT\n+\nIIII\n" > "$out/mitochondrion_1.fq"; printf "@r/2\nACGT\n+\nIIII\n" > "$out/mitochondrion_2.fq";;\n'
                'finalize) while [ "$#" -gt 0 ]; do if [ "$1" = "--out-dir" ]; then out="$2"; shift 2; else shift; fi; done; mkdir -p "$out"; printf "status\tcircular\n" > "$out/mitochondrial_assembly_summary.tsv"; printf ">mito\nACGTACGT\n" > "$out/mitochondrial_assembly.fasta";; esac\n'
            ); mito.chmod(0o755)
            gb = tmp / 'mito.gb'; gb.write_text('LOCUS TEST 4 bp DNA circular\nORIGIN\n        1 acgt\n//\n')
            reads = tmp / 'r1.fq'; reads.write_text('@r\nACGT\n+\nIIII\n')
            samples = tmp / 'samples.tsv'; samples.write_text(f'demo\t{reads}\t{reads}\n')
            output = tmp / 'out'
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--',
                 'mito', '-f', str(samples), '-o', str(output), '--mito-genbank', str(gb), '--mito-initial-reads', '1', '--mito-max-reads', '2'],
                cwd=ROOT, text=True,
                env={**os.environ, 'GENEMINER2_ENGINE': 'rust', 'GM2_COMPONENT_DIR': str(components)},
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue((output / '.gm2_mito_reference' / 'metadata' / 'mitochondrial_reference.fasta').is_file())
            self.assertTrue((output / '1_Demo' / 'filtered' / 'mitochondrion.fq').is_file())
            self.assertTrue((output / '.mito_adaptive' / '1m' / '1_Demo' / 'mito_rescue_round_1' / 'assembly_refs' / 'mitochondrion.fasta').is_file())
            self.assertTrue((output / '1_Demo' / 'mito' / 'mitochondrial_assembly_summary.tsv').is_file())

    def test_native_legacy_uce_filter_path_recovers_a_paired_synthetic_locus(self):
        rng = random.Random(20260723)
        truth = ''.join(rng.choice('ACGT') for _ in range(700))
        complement = str.maketrans('ACGT', 'TGCA')
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            references = tmp / 'references'; references.mkdir()
            (references / 'uce_demo.fasta').write_text(f'>uce_demo\n{truth}\n')
            r1, r2 = tmp / 'r1.fq', tmp / 'r2.fq'
            with r1.open('w') as first, r2.open('w') as second:
                for index, start in enumerate(range(0, 551, 10)):
                    read = truth[start:start + 150]
                    first.write(f'@read{index}/1\n{read}\n+\n' + 'I' * 150 + '\n')
                    second.write(f'@read{index}/2\n{read.translate(complement)[::-1]}\n+\n' + 'I' * 150 + '\n')
            sample_list = tmp / 'samples.tsv'; sample_list.write_text(f'demo\t{r1}\t{r2}\n')
            output = tmp / 'out'
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--',
                 'filter', 'refilter', 'assemble', '--assembly-mode', 'uce', '--legacy-uce-filter',
                 '-f', str(sample_list), '-r', str(references), '-o', str(output),
                 '-kf', '31', '-ka', '31', '--min-ka', '31', '--max-ka', '31'],
                cwd=ROOT, text=True,
                env={**os.environ, 'GENEMINER2_ENGINE': 'rust', 'GM2_COMPONENT_DIR': str(ROOT / 'cli' / 'bin')},
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=90,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            with (output / '1_Demo' / 'uce_assembly_summary.csv').open(newline='') as handle:
                row = next(csv.DictReader(handle))
            self.assertEqual(row['accepted'], '1')
            self.assertTrue((output / '1_Demo' / 'results' / 'uce_demo.fasta').is_file())

    def test_native_uce_rescue_runs_terminal_second_round_and_writes_reports(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); components = tmp / 'components'; components.mkdir()
            uce_filter = components / 'uce_filter'
            uce_filter.write_text(
                '#!/bin/sh\nwhile [ "$#" -gt 0 ]; do case "$1" in -o) out="$2"; shift 2;; *) shift;; esac; done\n'
                'mkdir -p "$out/filtered"; printf "@r\nACGT\n+\nIIII\n" > "$out/filtered/gene.fq"; printf "gene,1\n" > "$out/ref_reads_count_dict.txt"\n'
            ); uce_filter.chmod(0o755)
            assembler = components / 'main_assembler-rust'
            assembler.write_text(
                '#!/bin/sh\nwhile [ "$#" -gt 0 ]; do case "$1" in -o) out="$2"; shift 2;; -r) ref="$2"; shift 2;; *) shift;; esac; done\n'
                'case "$ref" in *uce_rescue_round_2*) len=180; reads=4;; *uce_rescue_round_1*) len=150; reads=3;; *) len=100; reads=1;; esac\n'
                'mkdir -p "$out/results" "$out/contigs_all"; seq=$(printf "%*s" "$len" "" | tr " " A); printf ">gene\n%s\n" "$seq" > "$out/results/gene.fasta"; cp "$out/results/gene.fasta" "$out/contigs_all/gene.fasta"; printf "locus,status,accepted,selected_contig_length,unique_read_count,read_count\ngene,success,1,%s,%s,%s\n" "$len" "$reads" "$reads" > "$out/uce_assembly_summary.csv"; printf "gene,success,%s,\n" "$len" > "$out/result_dict.txt"\n'
            ); assembler.chmod(0o755)
            refs = tmp / 'refs'; refs.mkdir(); (refs / 'gene.fasta').write_text('>gene\n' + 'A' * 100 + '\n')
            reads = tmp / 'reads.fq'; reads.write_text('@r\nACGT\n+\nIIII\n')
            samples = tmp / 'samples.tsv'; samples.write_text(f'demo\t{reads}\t{reads}\n')
            out = tmp / 'out'
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--', 'filter', 'refilter', 'assemble',
                 '--assembly-mode', 'uce', '--uce-rescue-reads', '--uce-rescue-rounds', '2', '-f', str(samples), '-r', str(refs), '-o', str(out)],
                cwd=ROOT, text=True, env={**os.environ, 'GENEMINER2_ENGINE': 'rust', 'GM2_COMPONENT_DIR': str(components)},
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue((out / '1_Demo' / 'uce_rescue_round_2' / 'terminal_baits' / 'gene.fasta').is_file())
            self.assertIn(',2,gene,terminal-only,', (out / '1_Demo' / 'uce_rescue_rounds.csv').read_text())
            self.assertIn('gene,success,1,180,4,4', (out / '1_Demo' / 'uce_assembly_summary.csv').read_text())
            self.assertTrue((out / '1_Demo' / 'uce_rescue_summary.csv').is_file())

    def test_native_uce_rescue_restores_only_locus_lost_by_rescue(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); components = tmp / 'components'; components.mkdir()
            uce_filter = components / 'uce_filter'
            uce_filter.write_text('#!/bin/sh\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-o" ]; then out="$2"; shift 2; else shift; fi; done\nmkdir -p "$out/filtered"; : > "$out/ref_reads_count_dict.txt"\n'); uce_filter.chmod(0o755)
            assembler = components / 'main_assembler-rust'
            assembler.write_text(
                '#!/bin/sh\nwhile [ "$#" -gt 0 ]; do case "$1" in -o) out="$2"; shift 2;; -r) ref="$2"; shift 2;; *) shift;; esac; done\n'
                'mkdir -p "$out/results" "$out/contigs_all"; case "$ref" in *uce_rescue_round_1*) printf ">a\n%s\n" "$(printf "%*s" 150 "" | tr " " A)" > "$out/results/a.fasta"; cp "$out/results/a.fasta" "$out/contigs_all/a.fasta"; printf "locus,status,accepted,selected_contig_length,unique_read_count,read_count\na,success,1,150,3,3\nb,failed,0,0,0,0\n" > "$out/uce_assembly_summary.csv";; *) for g in a b; do printf ">$g\n%s\n" "$(printf "%*s" 100 "" | tr " " A)" > "$out/results/$g.fasta"; cp "$out/results/$g.fasta" "$out/contigs_all/$g.fasta"; done; printf "locus,status,accepted,selected_contig_length,unique_read_count,read_count\na,success,1,100,1,1\nb,success,1,100,1,1\n" > "$out/uce_assembly_summary.csv";; esac; : > "$out/result_dict.txt"\n'
            ); assembler.chmod(0o755)
            refs = tmp / 'refs'; refs.mkdir(); (refs / 'a.fasta').write_text('>a\n' + 'A' * 100 + '\n'); (refs / 'b.fasta').write_text('>b\n' + 'A' * 100 + '\n')
            reads = tmp / 'reads.fq'; reads.write_text('@r\nACGT\n+\nIIII\n'); samples = tmp / 'samples.tsv'; samples.write_text(f'demo\t{reads}\t{reads}\n'); out = tmp / 'out'
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--', 'filter', 'refilter', 'assemble', '--assembly-mode', 'uce', '--uce-rescue-reads', '--uce-rescue-rounds', '1', '-f', str(samples), '-r', str(refs), '-o', str(out)],
                cwd=ROOT, text=True, env={**os.environ, 'GENEMINER2_ENGINE': 'rust', 'GM2_COMPONENT_DIR': str(components)}, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            summary = (out / '1_Demo' / 'uce_assembly_summary.csv').read_text()
            self.assertIn('a,success,1,150,3,3', summary)
            self.assertIn('b,success,1,100,1,1', summary)
            self.assertEqual(len(''.join((out / '1_Demo' / 'results' / 'b.fasta').read_text().splitlines()[1:])), 100)

    def test_native_uce_pipeline_recovers_a_paired_synthetic_locus(self):
        rng = random.Random(20260722)
        truth = ''.join(rng.choice('ACGT') for _ in range(700))
        complement = str.maketrans('ACGT', 'TGCA')
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            references = tmp / 'references'
            references.mkdir()
            (references / 'uce_demo.fasta').write_text(f'>uce_demo\n{truth}\n')
            r1, r2 = tmp / 'r1.fq', tmp / 'r2.fq'
            with r1.open('w') as first, r2.open('w') as second:
                for index, start in enumerate(range(0, 551, 10)):
                    read = truth[start:start + 150]
                    first.write(f'@read{index}/1\n{read}\n+\n' + 'I' * 150 + '\n')
                    second.write(f'@read{index}/2\n{read.translate(complement)[::-1]}\n+\n' + 'I' * 150 + '\n')
            sample_list = tmp / 'samples.tsv'
            sample_list.write_text(f'demo\t{r1}\t{r2}\n')
            proc = subprocess.run(
                ["cargo", "run", "--quiet", "--manifest-path", str(MANIFEST), "--", "filter", "refilter", "assemble", "--assembly-mode", "uce", "-f", str(sample_list), "-r", str(references), "-o", str(tmp / 'out'), "-kf", "31", "-ka", "31", "--min-ka", "31", "--max-ka", "31", "--uce-rescue-reads"],
                cwd=ROOT,
                text=True,
                env={**os.environ, "GENEMINER2_ENGINE": "rust", "GM2_COMPONENT_DIR": str(ROOT / 'cli' / 'bin')},
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=90,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            with (tmp / 'out' / '1_Demo' / 'uce_assembly_summary.csv').open(newline='') as handle:
                row = next(csv.DictReader(handle))
            self.assertEqual(row['accepted'], '1')
            self.assertGreaterEqual(int(row['selected_contig_length']), 650)
            self.assertTrue((tmp / 'out' / '1_Demo' / 'uce_rescue_round_1' / 'assembly_refs' / 'uce_demo.fasta').is_file())

    def test_native_original_pipeline_runs_without_gene_postprocessing(self):
        rng = random.Random(91); truth = ''.join(rng.choice('ACGT') for _ in range(500)); complement = str.maketrans('ACGT', 'TGCA')
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); references = tmp / 'references'; references.mkdir(); (references / 'locus.fasta').write_text(f'>locus\n{truth}\n')
            r1, r2 = tmp / 'r1.fq', tmp / 'r2.fq'
            with r1.open('w') as first, r2.open('w') as second:
                for index, start in enumerate(range(0, 351, 10)):
                    sequence = truth[start:start + 150]
                    first.write(f'@read{index}/1\n{sequence}\n+\n' + 'I' * len(sequence) + '\n')
                    second.write(f'@read{index}/2\n{sequence.translate(complement)[::-1]}\n+\n' + 'I' * len(sequence) + '\n')
            samples = tmp / 'samples.tsv'; samples.write_text(f'demo\t{r1}\t{r2}\n'); output = tmp / 'out'
            proc = subprocess.run(['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--', 'filter', 'refilter', 'assemble', '-f', str(samples), '-r', str(references), '-o', str(output), '-kf', '31', '-ka', '31', '--min-ka', '31', '--max-ka', '31', '--reuse-reference-cache', '--workflow-profile'], cwd=ROOT, text=True, env={**os.environ, 'GENEMINER2_ENGINE': 'rust', 'GM2_COMPONENT_DIR': str(ROOT / 'cli' / 'bin')}, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=90)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue((output / '1_Demo' / 'results' / 'locus.fasta').is_file()); self.assertTrue((output / '.gm2_reference_cache').is_dir()); self.assertTrue((output / 'workflow_profile.tsv').is_file()); self.assertFalse((output / 'gene').exists())
            profile = (output / 'workflow_profile.tsv').read_text()
            self.assertIn('__reference__\t0\tmainfilter_index\t', profile)
            self.assertIn('1_Demo\t0\tfilter\t', profile)
            self.assertIn('1_Demo\t0\trefilter\t', profile)
            self.assertIn('1_Demo\t0\tassemble\t', profile)

    def test_native_cleanup_runs_only_after_successful_complete_pipeline(self):
        rng = random.Random(92); truth = ''.join(rng.choice('ACGT') for _ in range(500)); complement = str.maketrans('ACGT', 'TGCA')
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); references = tmp / 'references'; references.mkdir(); (references / 'locus.fasta').write_text(f'>locus\n{truth}\n')
            r1, r2 = tmp / 'r1.fq', tmp / 'r2.fq'
            with r1.open('w') as first, r2.open('w') as second:
                for index, start in enumerate(range(0, 351, 10)):
                    sequence = truth[start:start + 150]
                    first.write(f'@r{index}/1\n{sequence}\n+\n' + 'I' * len(sequence) + '\n')
                    second.write(f'@r{index}/2\n{sequence.translate(complement)[::-1]}\n+\n' + 'I' * len(sequence) + '\n')
            samples = tmp / 'samples.tsv'; samples.write_text(f'demo\t{r1}\t{r2}\n'); output = tmp / 'out'
            proc = subprocess.run(['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--', 'filter', 'refilter', 'assemble', '--cleanup-intermediates', '-f', str(samples), '-r', str(references), '-o', str(output), '-kf', '31', '-ka', '31', '--min-ka', '31', '--max-ka', '31'], cwd=ROOT, text=True, env={**os.environ, 'GENEMINER2_ENGINE': 'rust', 'GM2_COMPONENT_DIR': str(ROOT / 'cli' / 'bin')}, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=90)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertFalse((output / '1_Demo' / 'filtered').exists()); self.assertFalse((output / '1_Demo' / 'filtered_pe').exists()); self.assertTrue((output / '1_Demo' / 'results' / 'locus.fasta').is_file()); self.assertIn('reproducible filtered reads', (output / 'cleanup_manifest.tsv').read_text())

    def test_native_gene_pipeline_runs_filter_to_cohort(self):
        rng = random.Random(77)
        truth = ''.join(rng.choice('ACGT') for _ in range(700))
        complement = str.maketrans('ACGT', 'TGCA')
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            references = tmp / 'references'
            references.mkdir()
            (references / 'gene_demo.fasta').write_text(f'>gene_demo\n{truth}\n')
            (tmp / 'out' / '99_Stale').mkdir(parents=True)
            r1, r2 = tmp / 'r1.fq', tmp / 'r2.fq'
            with r1.open('w') as first, r2.open('w') as second:
                for index, start in enumerate(range(0, 551, 10)):
                    read = truth[start:start + 150]
                    first.write(f'@read{index}/1\n{read}\n+\n' + 'I' * 150 + '\n')
                    second.write(f'@read{index}/2\n{read.translate(complement)[::-1]}\n+\n' + 'I' * 150 + '\n')
            sample_list = tmp / 'samples.tsv'
            sample_list.write_text(f'demo\t{r1}\t{r2}\n')
            proc = subprocess.run(
                ["cargo", "run", "--quiet", "--manifest-path", str(MANIFEST), "--", "gene", "-f", str(sample_list), "-r", str(references), "-o", str(tmp / 'out'), "-kf", "31", "-ka", "31", "--min-ka", "31", "--max-ka", "31"],
                cwd=ROOT,
                text=True,
                env={**os.environ, "GENEMINER2_ENGINE": "rust", "GM2_COMPONENT_DIR": str(ROOT / 'cli' / 'bin')},
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=120,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue((tmp / 'out' / '1_Demo' / 'results' / 'gene_demo.fasta').is_file())
            self.assertTrue((tmp / 'out' / 'gene' / 'family_summary.tsv').is_file())
            self.assertNotIn('99_Stale', (tmp / 'out' / 'gene' / 'family_count_matrix.tsv').read_text())

    def test_native_gene_pipeline_accepts_explicit_paired_reads(self):
        rng = random.Random(78); truth = ''.join(rng.choice('ACGT') for _ in range(700)); complement = str.maketrans('ACGT', 'TGCA')
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); references = tmp / 'references'; references.mkdir(); (references / 'gene_demo.fasta').write_text(f'>gene_demo\n{truth}\n')
            r1, r2 = tmp / 'r1.fq', tmp / 'r2.fq'
            with r1.open('w') as first, r2.open('w') as second:
                for index, start in enumerate(range(0, 551, 10)):
                    read = truth[start:start + 150]
                    first.write(f'@read{index}/1\n{read}\n+\n' + 'I' * 150 + '\n'); second.write(f'@read{index}/2\n{read.translate(complement)[::-1]}\n+\n' + 'I' * 150 + '\n')
            samples = tmp / 'samples.tsv'; samples.write_text(f'demo\t{r1}\t{r2}\n')
            proc = subprocess.run(['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--', 'gene', '-f', str(samples), '-r', str(references), '-o', str(tmp / 'out'), '-kf', '31', '-ka', '31', '--min-ka', '31', '--max-ka', '31'], cwd=ROOT, text=True, env={**os.environ, 'GENEMINER2_ENGINE': 'rust', 'GM2_COMPONENT_DIR': str(ROOT / 'cli' / 'bin')}, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=120)
            self.assertEqual(proc.returncode, 0, proc.stderr); self.assertTrue((tmp / 'out' / '1_Demo' / 'results' / 'gene_demo.fasta').is_file())

    def test_native_rejects_invalid_combine_parallelism_even_without_alignment(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            references = tmp / 'references'; references.mkdir()
            (references / 'locus.fasta').write_text('>locus\nACGT\n')
            samples = tmp / 'samples.tsv'; samples.write_text('sample\treads.fq\n')
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--',
                 'combine', '--no-alignment', '--filter-processes', '0',
                 '-f', str(samples), '-r', str(references), '-o', str(tmp / 'out')],
                cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn('--filter-processes must be at least 1', proc.stderr)

    def test_native_workflow_profile_is_written_for_failed_sample(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); components = tmp / 'components'; components.mkdir()
            filter_bin = components / 'MainFilterNew'
            filter_bin.write_text('#!/bin/sh\nfor arg in "$@"; do [ "$arg" = "-q1" ] && exit 1; done\nexit 0\n')
            filter_bin.chmod(0o755)
            references = tmp / 'references'; references.mkdir()
            (references / 'locus.fasta').write_text('>locus\nACGT\n')
            reads = tmp / 'reads.fq'; reads.write_text('@r\nACGT\n+\nIIII\n')
            samples = tmp / 'samples.tsv'; samples.write_text(f'sample\t{reads}\n')
            output = tmp / 'out'
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--',
                 'filter', '--workflow-profile', '-f', str(samples), '-r', str(references), '-o', str(output)],
                cwd=ROOT, text=True,
                env={**os.environ, 'GM2_COMPONENT_DIR': str(components)},
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertNotEqual(proc.returncode, 0)
            profile = (output / 'workflow_profile.tsv').read_text()
            self.assertIn('1_Sample\t0\tfilter\t', profile)
            self.assertIn('\tfailed\n', profile)

    def test_native_uce_rejects_incompatible_assembler_implementation(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            references = tmp / 'references'; references.mkdir()
            (references / 'uce.fasta').write_text('>uce\nACGT\n')
            samples = tmp / 'samples.tsv'; samples.write_text('demo\treads_1.fq\treads_2.fq\n')
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--',
                 'filter', 'assemble', '--assembly-mode', 'uce', '--assembler-implementation', 'original',
                 '-f', str(samples), '-r', str(references), '-o', str(tmp / 'out')],
                cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn('UCE assembly requires --assembler-implementation auto or uce-rust', proc.stderr)

    def test_native_uce_refilter_alone_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            reads, references = tmp / 'reads.fq', tmp / 'references'
            reads.write_text('@read/1\nACGTACGTACGTACGTACGTACGTACGTACG\n+\n' + 'I' * 32 + '\n')
            references.mkdir()
            (references / 'uce.fasta').write_text('>uce\nACGTACGTACGTACGTACGTACGTACGTACG\n')
            sample_list = tmp / 'samples.tsv'
            sample_list.write_text(f'demo\t{reads}\t{reads}\n')
            proc = subprocess.run(
                ["cargo", "run", "--quiet", "--manifest-path", str(MANIFEST), "--", "refilter", "--assembly-mode", "uce", "-f", str(sample_list), "-r", str(references), "-o", str(tmp / 'out')],
                cwd=ROOT,
                text=True,
                env={**os.environ, "GENEMINER2_ENGINE": "rust", "GM2_COMPONENT_DIR": str(ROOT / 'cli' / 'bin')},
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=30,
            )
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn('UCE refilter is fused', proc.stderr)

    def test_native_two_column_sample_list_uses_legacy_duplicate_mate_convention(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); components = tmp / 'components'; components.mkdir(); capture = tmp / 'uce.args'
            uce_filter = components / 'uce_filter'
            uce_filter.write_text('#!/bin/sh\nprintf "%s\n" "$@" > "$GM2_CAPTURE"\nmkdir -p "$4/filtered"\n')
            uce_filter.chmod(0o755)
            reads, references = tmp / 'reads.fq', tmp / 'references'
            reads.write_text('@read\nACGTACGTACGTACGTACGTACGTACGTACG\n+\n' + 'I' * 32 + '\n')
            references.mkdir(); (references / 'uce.fasta').write_text('>uce\nACGTACGTACGTACGTACGTACGTACGTACG\n')
            sample_list = tmp / 'samples.tsv'; sample_list.write_text(f'demo\t{reads}\n')
            proc = subprocess.run(
                ['cargo', 'run', '--quiet', '--manifest-path', str(MANIFEST), '--', 'filter', '--assembly-mode', 'uce', '-f', str(sample_list), '-r', str(references), '-o', str(tmp / 'out')],
                cwd=ROOT, text=True, env={**os.environ, 'GENEMINER2_ENGINE': 'rust', 'GM2_COMPONENT_DIR': str(components), 'GM2_CAPTURE': str(capture)},
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            received = capture.read_text().splitlines()
            self.assertEqual(received[received.index('-q1') + 1], str(reads))
            self.assertEqual(received[received.index('-q2') + 1], str(reads))


if __name__ == "__main__":
    unittest.main()
