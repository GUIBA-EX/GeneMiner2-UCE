const RUN_LEN_CONST: f64 = 0.577_215_664_9_f64 / std::f64::consts::LN_2 - 1.5;
const THR_P95_2T: f64 = 1.96;
const THR_1E5_1T: f64 = 3.74;
const TOLERANCE: f64 = 1e-5;

pub fn collect_runs_stats(orientations: &[u8]) -> [usize; 13] {
    let mut results = [0_usize; 13];
    if orientations.is_empty() {
        return results;
    }
    let mut curr_dir = 0_usize;
    let mut curr_len = 0_usize;
    results[12] = orientations.len();
    for &value in orientations {
        let orient = value as usize;
        if orient != curr_dir {
            results[curr_dir] = results[curr_dir].max(curr_len);
            results[curr_dir + 4] += 1;
            curr_dir = orient;
            curr_len = 0;
        }
        if curr_dir != 0 {
            curr_len += 1;
            results[curr_dir + 8] += 1;
        }
    }
    results[curr_dir] = results[curr_dir].max(curr_len);
    results[curr_dir + 4] += 1;
    results
}

fn is_close(a: f64, b: f64, abs_tol: f64) -> bool {
    let rel_tol = 1e-9;
    (a - b).abs() <= abs_tol.max(rel_tol * a.abs().max(b.abs()))
}

pub fn infer_orientation(stats: [usize; 13]) -> u8 {
    let fwd_l = stats[1] as f64;
    let rev_l = stats[2] as f64;
    let fwd_r = stats[5] as f64;
    let rev_r = stats[6] as f64;
    let fwd_n = stats[9] as f64;
    let rev_n = stats[10] as f64;
    let amb_n = stats[11] as f64;
    let tot_n = stats[12] as f64;
    if fwd_n <= 1.0 {
        if rev_n <= 1.0 {
            return if amb_n <= 1.0 { 0 } else { 3 };
        }
        return 2;
    } else if rev_n <= 1.0 {
        return 1;
    }
    let npr = 2.0 * fwd_n * rev_n;
    let nht = fwd_n + rev_n;
    let erc = npr / nht + 1.0;
    let vrn = npr * (npr - nht) / (nht * nht * (nht - 1.0));
    if (fwd_r + rev_r - erc) / vrn.sqrt() > -THR_1E5_1T {
        let ntt = if fwd_r > rev_r {
            fwd_n + fwd_r - rev_r
        } else {
            fwd_n + rev_r - fwd_r
        };
        let rex = fwd_n / ntt;
        if is_close(rex, 1.0, TOLERANCE)
            || (1.0 - rex) / (rex * (1.0 - rex) / ntt).sqrt() < THR_P95_2T
        {
            return 0;
        }
    }
    let erl = (tot_n.log2() + RUN_LEN_CONST).max(0.0) + 4.0;
    let mut orient = ((fwd_l > erl) as u8) + ((rev_l > erl) as u8) * 2;
    if orient != 3 {
        return orient;
    }
    let lpf = (1.0 / (1.0 - fwd_l + fwd_n))
        .ln()
        .mul_add(1.0 / fwd_l, 0.0)
        .exp();
    let lpr = (1.0 / (1.0 - rev_l + rev_n))
        .ln()
        .mul_add(1.0 / rev_l, 0.0)
        .exp();
    let fpz = is_close(lpf, 0.0, TOLERANCE);
    let rpz = is_close(lpr, 0.0, TOLERANCE);
    if fpz {
        orient = 2 - 2 * (rpz as u8);
    } else if rpz {
        orient = 1;
    } else if (is_close(lpf, 1.0, TOLERANCE) && is_close(lpr, 1.0, TOLERANCE))
        || (lpf - lpr).abs()
            / (lpf.powi(2) * (1.0 - lpf) / fwd_n + lpr.powi(2) * (1.0 - lpr) / rev_n).sqrt()
            < THR_P95_2T
    {
        orient = 0;
    }
    orient
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_statistics_preserve_position_and_ambiguity() {
        let stats = collect_runs_stats(&[0, 1, 1, 0, 2, 2, 2, 3, 0]);
        assert_eq!(stats[1], 2);
        assert_eq!(stats[2], 3);
        assert_eq!(stats[9], 2);
        assert_eq!(stats[10], 3);
        assert_eq!(stats[11], 1);
        assert_eq!(stats[12], 9);
    }

    #[test]
    fn empty_orientations_are_unmapped() {
        assert_eq!(infer_orientation(collect_runs_stats(&[])), 0);
    }
}
