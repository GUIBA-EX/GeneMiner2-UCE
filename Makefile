PY_SRC := build_consensus merge_seq unix_command
PY_BIN := $(patsubst %,cli/bin/%,$(PY_SRC))
RUST_ASSEMBLER_BIN := cli/bin/main_assembler-rust
ORIGINAL_ASSEMBLER_BIN := cli/bin/main_assembler-original
ORIGINAL_ASSEMBLER_SOURCE := scripts/main_assembler_original.py
ASSEMBLER_RUST_MANIFEST := rust/main_assembler/Cargo.toml
ASSEMBLER_RUST_SOURCES := $(ASSEMBLER_RUST_MANIFEST) $(wildcard rust/main_assembler/src/*.rs)
REFILTER_BIN := cli/bin/main_refilter_new
FILTER_RUST_MANIFEST := rust/main_filter_new/Cargo.toml
FILTER_RUST_SOURCES := $(FILTER_RUST_MANIFEST) rust/main_filter_new/src/main.rs
POPULATION_BIN := cli/bin/main_population
POPULATION_RUST_MANIFEST := rust/main_population/Cargo.toml
POPULATION_RUST_SOURCES := $(POPULATION_RUST_MANIFEST) rust/main_population/src/main.rs
FILTER_HAXE_SOURCES := $(wildcard scripts/filter/*.h scripts/filter/*.hpp scripts/filter/*.hx)

.PHONY: build clean cython distclean haxe-filter rust-assembler

build: cli/bin/MainFilterNew $(REFILTER_BIN) $(ORIGINAL_ASSEMBLER_BIN) $(POPULATION_BIN) $(PY_BIN)
	for target in $(PY_SRC); do cp -L -r -t cli/bin --reflink=auto --update=none scripts/dist/$$target/_internal; done
	if command -v cargo >/dev/null 2>&1; then $(MAKE) $(RUST_ASSEMBLER_BIN); fi
	cd cli && ln -f -r -s bin/unix_command geneminer2

clean:
	rm -f -r scripts/filter/bin
	rm -f -r scripts/build
	rm -f -r scripts/dist
	rm -f -r rust/main_filter_new/target
	rm -f -r rust/main_refilter_new/target
	rm -f -r rust/main_assembler/target
	rm -f -r rust/main_population/target

distclean: clean
	for target in $(PY_SRC); do rm -f scripts/$$target.spec; done
	rm -f scripts/main_refilter_new.spec
	rm -f scripts/main_assembler_original.spec
	rm -f -r cli/bin
	rm -f cli/geneminer2

cli/bin:
	mkdir -p cli/bin

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

cython:
	cd scripts && cythonize -i main_refilter_ext.pyx

$(REFILTER_BIN): scripts/main_refilter_new.py rust/main_refilter_new/Cargo.toml rust/main_refilter_new/src/main.rs | cython
	if command -v cargo >/dev/null 2>&1; then \
		cargo build --release --manifest-path rust/main_refilter_new/Cargo.toml; \
		install -D -t cli/bin rust/main_refilter_new/target/release/main_refilter_new; \
	else \
		(cd scripts && pyinstaller -D -y --optimize 2 main_refilter_new.py); \
		install -D -t cli/bin scripts/dist/main_refilter_new/main_refilter_new; \
		cp -L -r -t cli/bin --reflink=auto --update=none scripts/dist/main_refilter_new/_internal; \
	fi

$(ORIGINAL_ASSEMBLER_BIN): $(ORIGINAL_ASSEMBLER_SOURCE) | cython cli/bin
	cd scripts && pyinstaller -D -y --optimize 2 main_assembler_original.py
	install scripts/dist/main_assembler_original/main_assembler_original $(ORIGINAL_ASSEMBLER_BIN)
	cp -L -r -t cli/bin --reflink=auto --update=none scripts/dist/main_assembler_original/_internal

rust-assembler:
	command -v cargo >/dev/null 2>&1 || { echo "Cargo is required for the optional Rust assembler" >&2; exit 1; }
	$(MAKE) $(RUST_ASSEMBLER_BIN)

$(RUST_ASSEMBLER_BIN): $(ASSEMBLER_RUST_SOURCES) | cli/bin
	cargo build --release --manifest-path $(ASSEMBLER_RUST_MANIFEST)
	install rust/main_assembler/target/release/main_assembler $(RUST_ASSEMBLER_BIN)

cli/bin/unix_command: scripts/gm2_stats.py

$(POPULATION_BIN): $(POPULATION_RUST_SOURCES) | cli/bin
	cargo build --release --manifest-path $(POPULATION_RUST_MANIFEST)
	install rust/main_population/target/release/main_population $(POPULATION_BIN)

$(PY_BIN): cli/bin/%: scripts/%.py | cython
	cd scripts && pyinstaller -D -y --optimize 2 $(notdir $<)
	install -D -t cli/bin scripts/dist/$(notdir $@)/$(notdir $@)
