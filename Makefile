CONSENSUS_BIN := cli/bin/build_consensus
CONSENSUS_RUST_MANIFEST := rust/build_consensus/Cargo.toml
CONSENSUS_RUST_SOURCES := $(CONSENSUS_RUST_MANIFEST) $(wildcard rust/build_consensus/src/*.rs)
RUST_ASSEMBLER_BIN := cli/bin/main_assembler-rust
ORIGINAL_RUST_ASSEMBLER_BIN := cli/bin/main_assembler-original-rust
ORIGINAL_RUST_ASSEMBLER_MANIFEST := rust/main_assembler_original/Cargo.toml
ORIGINAL_RUST_ASSEMBLER_SOURCES := $(ORIGINAL_RUST_ASSEMBLER_MANIFEST) $(wildcard rust/main_assembler_original/src/*.rs)
ASSEMBLER_RUST_MANIFEST := rust/main_assembler/Cargo.toml
ASSEMBLER_RUST_SOURCES := $(ASSEMBLER_RUST_MANIFEST) $(wildcard rust/main_assembler/src/*.rs)
REFILTER_BIN := cli/bin/main_refilter_new
UCE_FILTER_BIN := cli/bin/uce_filter
UCE_FILTER_MANIFEST := rust/uce_filter/Cargo.toml
UCE_FILTER_SOURCES := $(UCE_FILTER_MANIFEST) $(wildcard rust/uce_filter/src/*.rs) $(wildcard rust/uce_filter_core/src/*.rs) rust/uce_filter_core/Cargo.toml
FILTER_RUST_MANIFEST := rust/main_filter_new/Cargo.toml
FILTER_RUST_SOURCES := $(FILTER_RUST_MANIFEST) rust/main_filter_new/src/main.rs
POPULATION_BIN := cli/bin/main_population
POPULATION_RUST_MANIFEST := rust/main_population/Cargo.toml
POPULATION_RUST_SOURCES := $(POPULATION_RUST_MANIFEST) $(wildcard rust/main_population/src/*.rs) $(wildcard rust/main_population/src/panref/*.rs)
TOOLS_RUST_MANIFEST := rust/gm2_tools/Cargo.toml
TOOLS_RUST_SOURCES := $(TOOLS_RUST_MANIFEST) $(wildcard rust/gm2_tools/src/*.rs) $(wildcard rust/gm2_tools/src/bin/*.rs)
ALIGNMENT_CLEAN_BIN := cli/bin/fix_alignment
MERGE_SEQ_BIN := cli/bin/merge_seq
BUILD_TRIMED_BIN := cli/bin/build_trimed
GM2_STATS_BIN := cli/bin/gm2_stats
MITO_WORKFLOW_BIN := cli/bin/mito_workflow
GENE_WORKFLOW_BIN := cli/bin/gene_workflow
RAD_WORKFLOW_BIN := cli/bin/rad_workflow
MARKER_PROFILE_BIN := cli/bin/marker_profile
REPEAT_BIN := cli/bin/main_repeat
RUST_CLI_BIN := cli/bin/geneminer2-rust
RUST_CLI_MANIFEST := rust/geneminer2_cli/Cargo.toml
RUST_CLI_SOURCES := $(RUST_CLI_MANIFEST) $(wildcard rust/geneminer2_cli/src/*.rs)
REPEAT_RUST_MANIFEST := rust/main_repeat/Cargo.toml
REPEAT_RUST_SOURCES := $(REPEAT_RUST_MANIFEST) $(wildcard rust/main_repeat/src/*.rs)
FILTER_HAXE_SOURCES := $(wildcard scripts/filter/*.h scripts/filter/*.hpp scripts/filter/*.hx)

.PHONY: build clean distclean haxe-filter rust-assembler

build: $(CONSENSUS_BIN) cli/bin/MainFilterNew $(REFILTER_BIN) $(UCE_FILTER_BIN) $(ORIGINAL_RUST_ASSEMBLER_BIN) $(RUST_ASSEMBLER_BIN) $(POPULATION_BIN) $(ALIGNMENT_CLEAN_BIN) $(MERGE_SEQ_BIN) $(BUILD_TRIMED_BIN) $(GM2_STATS_BIN) $(MARKER_PROFILE_BIN) $(MITO_WORKFLOW_BIN) $(GENE_WORKFLOW_BIN) $(RAD_WORKFLOW_BIN) $(REPEAT_BIN) $(RUST_CLI_BIN)
	cd cli && ln -sfn -r bin/geneminer2-rust geneminer2

clean:
	rm -f -r scripts/filter/bin
	rm -f -r rust/main_assembler_original/target
	rm -f -r rust/main_filter_new/target
	rm -f -r rust/build_consensus/target
	rm -f -r rust/main_refilter_new/target
	rm -f -r rust/uce_filter/target
	rm -f -r rust/uce_filter_core/target
	rm -f -r rust/main_assembler/target
	rm -f -r rust/main_population/target
	rm -f -r rust/main_repeat/target
	rm -f -r rust/gm2_tools/target

distclean: clean
	rm -f -r cli/bin
	rm -f cli/geneminer2

cli/bin:
	mkdir -p cli/bin

$(CONSENSUS_BIN): $(CONSENSUS_RUST_SOURCES) | cli/bin
	cargo build --release --manifest-path $(CONSENSUS_RUST_MANIFEST)
	install rust/build_consensus/target/release/build_consensus $(CONSENSUS_BIN)

cli/bin/MainFilterNew: $(FILTER_RUST_SOURCES) $(FILTER_HAXE_SOURCES) | cli/bin
	if command -v cargo >/dev/null 2>&1; then \
		cargo build --release --manifest-path $(FILTER_RUST_MANIFEST); \
		install rust/main_filter_new/target/release/MainFilterNew cli/bin/MainFilterNew; \
	else \
		(cd scripts/filter && haxe -cpp bin -dce full -D analyzer-optimize -D HXCPP_GC_BIG_BLOCKS -D HXCPP_GC_MOVING -D HXCPP_M64 -D HXCPP_OPTIMIZE_LINK -D HXCPP_SINGLE_THREADED_APP -D HXCPP_VISIT_ALLOCS -main MainFilterNew.hx); \
		install scripts/filter/bin/MainFilterNew cli/bin/MainFilterNew; \
	fi

haxe-filter: cli/bin/MainFilterNew-haxe

cli/bin/MainFilterNew-haxe: $(FILTER_HAXE_SOURCES) | cli/bin
	cd scripts/filter && haxe -cpp bin -dce full -D analyzer-optimize -D HXCPP_GC_BIG_BLOCKS -D HXCPP_GC_MOVING -D HXCPP_M64 -D HXCPP_OPTIMIZE_LINK -D HXCPP_SINGLE_THREADED_APP -D HXCPP_VISIT_ALLOCS -main MainFilterNew.hx
	install scripts/filter/bin/MainFilterNew cli/bin/MainFilterNew-haxe
$(ORIGINAL_RUST_ASSEMBLER_BIN): $(ORIGINAL_RUST_ASSEMBLER_SOURCES) | cli/bin
	command -v cargo >/dev/null 2>&1 || { echo "Cargo is required for the Rust assembler" >&2; exit 1; }
	cargo build --release --manifest-path $(ORIGINAL_RUST_ASSEMBLER_MANIFEST)
	install rust/main_assembler_original/target/release/main_assembler_original $(ORIGINAL_RUST_ASSEMBLER_BIN)


$(REFILTER_BIN): rust/main_refilter_new/Cargo.toml rust/main_refilter_new/src/main.rs | cli/bin
	cargo build --release --manifest-path rust/main_refilter_new/Cargo.toml
	install -D -t cli/bin rust/main_refilter_new/target/release/main_refilter_new

$(UCE_FILTER_BIN): $(UCE_FILTER_SOURCES) | cli/bin
	cargo build --release --manifest-path $(UCE_FILTER_MANIFEST)
	install rust/uce_filter/target/release/uce_filter $(UCE_FILTER_BIN)

rust-assembler:
	command -v cargo >/dev/null 2>&1 || { echo "Cargo is required for the Rust assembler" >&2; exit 1; }
	$(MAKE) $(RUST_ASSEMBLER_BIN)

$(RUST_ASSEMBLER_BIN): $(ASSEMBLER_RUST_SOURCES) | cli/bin
	command -v cargo >/dev/null 2>&1 || { echo "Cargo is required for the Rust assembler" >&2; exit 1; }
	cargo build --release --manifest-path $(ASSEMBLER_RUST_MANIFEST)
	install rust/main_assembler/target/release/main_assembler $(RUST_ASSEMBLER_BIN)


$(POPULATION_BIN): $(POPULATION_RUST_SOURCES) | cli/bin
	cargo build --release --manifest-path $(POPULATION_RUST_MANIFEST)
	install rust/main_population/target/release/main_population $(POPULATION_BIN)

$(ALIGNMENT_CLEAN_BIN): $(TOOLS_RUST_SOURCES) | cli/bin
	cargo build --release --manifest-path $(TOOLS_RUST_MANIFEST) --bin fix_alignment
	install rust/gm2_tools/target/release/fix_alignment $(ALIGNMENT_CLEAN_BIN)

$(MERGE_SEQ_BIN): $(TOOLS_RUST_SOURCES) | cli/bin
	cargo build --release --manifest-path $(TOOLS_RUST_MANIFEST) --bin merge_seq
	install rust/gm2_tools/target/release/merge_seq $(MERGE_SEQ_BIN)

$(BUILD_TRIMED_BIN): $(TOOLS_RUST_SOURCES) | cli/bin
	cargo build --release --manifest-path $(TOOLS_RUST_MANIFEST) --bin build_trimed
	install rust/gm2_tools/target/release/build_trimed $(BUILD_TRIMED_BIN)

$(GM2_STATS_BIN): $(TOOLS_RUST_SOURCES) | cli/bin
	cargo build --release --manifest-path $(TOOLS_RUST_MANIFEST) --bin gm2_stats
	install rust/gm2_tools/target/release/gm2_stats $(GM2_STATS_BIN)

$(MITO_WORKFLOW_BIN): $(TOOLS_RUST_SOURCES) | cli/bin
	cargo build --release --manifest-path $(TOOLS_RUST_MANIFEST) --bin mito_workflow
	install rust/gm2_tools/target/release/mito_workflow $(MITO_WORKFLOW_BIN)

$(GENE_WORKFLOW_BIN): $(TOOLS_RUST_SOURCES) | cli/bin
	cargo build --release --manifest-path $(TOOLS_RUST_MANIFEST) --bin gene_workflow
	install rust/gm2_tools/target/release/gene_workflow $(GENE_WORKFLOW_BIN)

$(RAD_WORKFLOW_BIN): $(TOOLS_RUST_SOURCES) | cli/bin
	cargo build --release --manifest-path $(TOOLS_RUST_MANIFEST) --bin rad_workflow
	install rust/gm2_tools/target/release/rad_workflow $(RAD_WORKFLOW_BIN)

$(MARKER_PROFILE_BIN): rust/marker_profile/Cargo.toml rust/marker_profile/src/main.rs | cli/bin
	cargo build --release --manifest-path rust/marker_profile/Cargo.toml
	install rust/marker_profile/target/release/marker_profile $(MARKER_PROFILE_BIN)

$(REPEAT_BIN): $(REPEAT_RUST_SOURCES) | cli/bin
	cargo build --release --manifest-path $(REPEAT_RUST_MANIFEST)
	install rust/main_repeat/target/release/main_repeat $(REPEAT_BIN)

$(RUST_CLI_BIN): $(RUST_CLI_SOURCES) | cli/bin
	cargo build --release --manifest-path $(RUST_CLI_MANIFEST)
	install rust/geneminer2_cli/target/release/geneminer2_cli $(RUST_CLI_BIN)
