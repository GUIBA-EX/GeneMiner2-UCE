//! Fast, non-cryptographic hash tables used by the assembler.
//!
//! AHash keeps k-mer and graph lookups cache-friendly without changing any
//! sequence, graph, or path-selection semantics.

pub use ahash::{AHashMap as HashMap, AHashSet as HashSet};
