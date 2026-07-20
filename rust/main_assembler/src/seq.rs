use crate::hash::HashMap;

// DNA 四个字母压成 2-bit；碰到含糊字符就停，别往图里掺水。
pub fn base_bits(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

pub fn bits_base(bits: u8) -> u8 {
    b"ACGT"[(bits & 3) as usize]
}

pub fn kmer_mask(k: usize) -> u128 {
    if k >= 64 {
        u128::MAX
    } else {
        (1_u128 << (2 * k)) - 1
    }
}

// 一个 k-mer 塞进 u128，k 不超过 63 时既快又省地儿。
pub fn encode_kmer(sequence: &[u8]) -> Option<u128> {
    let mut encoded = 0_u128;
    for &base in sequence {
        encoded = (encoded << 2) | u128::from(base_bits(base)?);
    }
    Some(encoded)
}

pub fn decode_kmer(mut encoded: u128, k: usize) -> Vec<u8> {
    let mut sequence = vec![b'A'; k];
    for base in sequence.iter_mut().rev() {
        *base = bits_base(encoded as u8);
        encoded >>= 2;
    }
    sequence
}

// 不解码成字符串，直接在位上翻正反链，热循环里能省不少功夫。
pub fn reverse_complement_kmer(mut encoded: u128, k: usize) -> u128 {
    let mut reverse = 0_u128;
    for _ in 0..k {
        reverse = (reverse << 2) | ((encoded as u8 ^ 3) & 3) as u128;
        encoded >>= 2;
    }
    reverse
}

pub fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

// 含糊碱基会切断连续区段，k-mer 只能在干净区段里滚。
pub fn valid_runs(sequence: &[u8]) -> Vec<(usize, &[u8])> {
    let mut runs = Vec::new();
    let mut start = None;
    for (index, &base) in sequence.iter().enumerate() {
        if base_bits(base).is_some() {
            start.get_or_insert(index);
        } else if let Some(run_start) = start.take() {
            runs.push((run_start, &sequence[run_start..index]));
        }
    }
    if let Some(run_start) = start {
        runs.push((run_start, &sequence[run_start..]));
    }
    runs
}

// 滚动更新正反两个编码，别每挪一步都重新切片算一遍。
pub fn for_each_kmer<F>(run: &[u8], k: usize, mut visit: F)
where
    F: FnMut(usize, u128, u128),
{
    if k == 0 || k > 63 || run.len() < k {
        return;
    }
    let mask = kmer_mask(k);
    let suffix_mask = kmer_mask(k - 1);
    let mut forward = match encode_kmer(&run[..k]) {
        Some(value) => value,
        None => return,
    };
    let mut reverse = reverse_complement_kmer(forward, k);
    visit(0, forward, reverse);

    for start in 1..=run.len() - k {
        let next = base_bits(run[start + k - 1]).expect("valid run");
        forward = ((forward & suffix_mask) << 2) | u128::from(next);
        reverse = (u128::from(next ^ 3) << (2 * (k - 1))) | (reverse >> 2);
        debug_assert_eq!(forward & !mask, 0);
        visit(start, forward, reverse);
    }
}

pub fn median(values: &[i64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 1 {
        sorted[mid] as f64
    } else {
        (sorted[mid - 1] as f64 + sorted[mid] as f64) / 2.0
    }
}

pub fn quartiles(values: &[i64]) -> (f64, f64, f64, i64) {
    if values.is_empty() {
        return (0.0, 0.0, 0.0, 1);
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let mid = sorted.len() / 2;
    let (left, right, q2) = if sorted.len() % 2 == 1 {
        (&sorted[..mid], &sorted[mid + 1..], sorted[mid] as f64)
    } else {
        (
            &sorted[..mid],
            &sorted[mid..],
            (sorted[mid - 1] as f64 + sorted[mid] as f64) / 2.0,
        )
    };
    (
        median(left),
        q2,
        median(right),
        sorted.last().copied().unwrap_or(0) + 1,
    )
}

pub fn calculate_auto_k(
    reference_sequences: &[Vec<u8>],
    read_slices: &HashMap<Vec<u8>, u64>,
    slice_len: usize,
    mut k_min: usize,
    k_max: usize,
    error_limit: u32,
) -> usize {
    if k_min.is_multiple_of(2) {
        k_min += 1;
    }
    if slice_len <= k_min || k_max < k_min || k_max > 63 {
        return k_min;
    }

    let mut counts: HashMap<u128, u32> = HashMap::new();
    for sequence in read_slices.keys() {
        for (_, run) in valid_runs(sequence) {
            for_each_kmer(run, k_min, |_, forward, reverse| {
                *counts.entry(forward).or_default() += 1;
                *counts.entry(reverse).or_default() += 1;
            });
        }
    }
    counts.retain(|_, depth| *depth > error_limit);

    let range = k_max - k_min + 1;
    let mut run_stats = vec![0_u64; range];
    for sequence in reference_sequences {
        for (_, run) in valid_runs(sequence) {
            if run.len() < k_min {
                continue;
            }
            let mut lengths = vec![0_usize];
            for_each_kmer(run, k_min, |_, forward, _| {
                let last = lengths.len() - 1;
                if counts.contains_key(&forward) {
                    lengths[last] += 1;
                    if lengths[last] >= range {
                        lengths.push(range / 2);
                    }
                } else if lengths[last] != 0 {
                    lengths.push(0);
                }
            });

            let mut frequencies: HashMap<usize, u64> = HashMap::new();
            for length in lengths {
                *frequencies.entry(length).or_default() += 1;
            }
            for (length, frequency) in frequencies {
                if length == 0 {
                    continue;
                }
                let mut index = length - 1;
                index -= index % 2;
                if index >= run_stats.len() {
                    index = run_stats.len() - 1;
                    index -= index % 2;
                }
                run_stats[index] += frequency;
                for offset in (2..=index).step_by(2) {
                    run_stats[index - offset] += frequency;
                }
            }
        }
    }

    let upper = match run_stats.iter().rposition(|count| *count > 0) {
        Some(index) => k_min + index,
        None => return k_min,
    };
    let lower = upper.div_ceil(2);
    let candidates: Vec<(usize, f64)> = run_stats
        .iter()
        .enumerate()
        .map(|(index, count)| {
            let k = k_min + index;
            let coverage = if slice_len > k {
                *count as f64 * k as f64 / (slice_len - k + 1) as f64
            } else {
                0.0
            };
            (k, coverage)
        })
        .filter(|(k, _)| *k > lower && *k <= upper)
        .collect();
    if candidates.is_empty() {
        return k_min;
    }
    let cutoff = candidates
        .iter()
        .map(|(_, coverage)| *coverage)
        .fold(0.0_f64, f64::max)
        / 2.0;
    candidates
        .iter()
        .rev()
        .find(|(_, coverage)| *coverage > cutoff)
        .map(|(k, _)| *k)
        .unwrap_or(k_min)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn u128_supports_k51_roundtrip() {
        let sequence: Vec<u8> = (0..51).map(|index| b"ACGT"[index % 4]).collect();
        assert_eq!(sequence.len(), 51);
        let encoded = encode_kmer(&sequence).unwrap();
        assert_eq!(decode_kmer(encoded, 51), sequence);
    }

    #[test]
    fn rolling_forward_and_reverse_match_direct_encoding() {
        let sequence = b"AACCGTTAGC";
        for_each_kmer(sequence, 5, |start, forward, reverse| {
            let direct = encode_kmer(&sequence[start..start + 5]).unwrap();
            assert_eq!(forward, direct);
            assert_eq!(reverse, reverse_complement_kmer(direct, 5));
        });
    }

    #[test]
    fn ambiguity_splits_runs() {
        let runs = valid_runs(b"AACNNGTT");
        assert_eq!(runs.len(), 2);
        assert_eq!(runs[0], (0, &b"AAC"[..]));
        assert_eq!(runs[1], (5, &b"GTT"[..]));
    }

    #[test]
    fn reverse_complement_preserves_ambiguity_like_python() {
        assert_eq!(reverse_complement(b"ACNUT"), b"AUNGT");
    }
}
