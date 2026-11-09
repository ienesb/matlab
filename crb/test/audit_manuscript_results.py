"""Read-only audit of historical and completed fresh-seed receiver results.

Run with the workspace SciPy directory on PYTHONPATH. Prints JSON to stdout;
never rewrites MAT files and never treats checkpoint files as completed data.
Bootstrap resampling preserves the common-trial pairing of each comparison.
"""
import json
from pathlib import Path
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
FILES = [
    "mcrb_rx_eva_reduced_n100_n100.mat",
    "mcrb_rx_eva_reduced_fast_n60_fast.mat",
    "mcrb_rx_eva_reduced_n100_finalization20260914_seed9101.mat",
    "mcrb_rx_eva_reduced_fast_n100_finalization20260914_seed9101.mat",
]
RNG = np.random.default_rng(20260914)


def ratio_ci(numerator, denominator):
    if not len(numerator):
        return {"n": 0, "ratio": None, "ci95": None}
    ix = RNG.integers(0, len(numerator), size=(5000, len(numerator)))
    vals = np.sqrt(np.median(numerator[ix], axis=1) /
                   np.median(denominator[ix], axis=1))
    return {"n": len(numerator),
            "ratio": float(np.sqrt(np.median(numerator) / np.median(denominator))),
            "ci95": np.quantile(vals, [.025, .975]).tolist()}


def audit(path):
    r = loadmat(path, simplify_cells=True)
    t = r["trial"]
    tr = t["tr"]
    names = list(r["trn"])
    snrs = np.atleast_1d(r["snr_db"])
    wfs = list(r["waveforms"])
    control_idx = [q for q, n in enumerate(names) if n.endswith("_bsc")]
    real_idx = [q for q, n in enumerate(names)
                if not n.endswith("_bsc") and not n.endswith("_pil")]
    valid = ((t["pcvalid"] > 0) & (t["exitflag"] > 0)
             & np.isfinite(t["kappa"]) & (t["kappa"] <= r["opts"]["kappa_thresh"])
             & np.all(np.isfinite(tr[:, real_idx, :, :]), axis=1))
    explicit_control_pcvalid = "pcvalid_bsc" in t
    control_pcvalid = (t["pcvalid_bsc"] > 0 if explicit_control_pcvalid
                       else t["kappa_bsc"] < 1e12)
    vb = (control_pcvalid & (t["exitflag_bsc"] > 0)
          & np.isfinite(t["kappa_bsc"])
          & (t["kappa_bsc"] <= r["opts"]["kappa_thresh"])
          & np.all(np.isfinite(tr[:, control_idx, :, :]), axis=1))
    paired = valid & vb
    out = {"file": path.name, "trials": int(r["nT"]), "seed": int(r["opts"]["seed"]),
           "seconds": float(r["seconds"]), "snr_db": snrs.tolist(),
           "control_pcvalid_saved": explicit_control_pcvalid,
           "negative_trace_count": int(np.count_nonzero(tr < 0)),
           "min_finite_trace": float(np.min(tr[np.isfinite(tr)])),
           "max_drop_real": (1-valid.mean(axis=0)).max(axis=0).tolist(),
           "max_drop_control": (1-vb.mean(axis=0)).max(axis=0).tolist(),
           "max_drop_paired": (1-paired.mean(axis=0)).max(axis=0).tolist(),
           "invalid_control_on_valid_real_draws": np.sum(valid & ~vb, axis=(0, 1)).tolist(),
           "max_score_floor_over_noise": np.nanmax(t["F_over_G"], axis=(0, 1)).tolist(),
           "points": []}
    for snr in [4, 8]:
        found = np.flatnonzero(snrs == snr)
        if not len(found):
            continue
        si = found[0]
        point = {"snr_db": snr, "waveforms": {}, "otfs_over_ofdm": {}}
        for iw, wf in enumerate(wfs):
            vr = valid[:, si, iw]
            vp = paired[:, si, iw]
            arm = {"n_valid": int(vr.sum()), "n_control_valid": int(vb[:, si, iw].sum()),
                   "n_paired": int(vp.sum()), "ber": float(r["ber"][si, iw])}
            for blk, scale in [("tau", 1e9), ("nu", 1)]:
                q = names.index("tr_" + blk + "_lb")
                qc = names.index("tr_" + blk + "_lb_bsc")
                arm[blk + "_lb"] = float(scale*np.sqrt(np.median(tr[vr, q, si, iw])))
                arm[blk + "_real_over_control"] = ratio_ci(tr[vp, q, si, iw], tr[vp, qc, si, iw])
            point["waveforms"][wf] = arm
        io, it = wfs.index("cp-ofdm"), wfs.index("cp-otfs")
        vw = valid[:, si, io] & valid[:, si, it]
        for blk in ["tau", "nu"]:
            q = names.index("tr_" + blk + "_lb")
            point["otfs_over_ofdm"][blk] = ratio_ci(tr[vw, q, si, it], tr[vw, q, si, io])
        out["points"].append(point)
    return out


if __name__ == "__main__":
    reports = [audit(ROOT / "results" / fn) for fn in FILES if (ROOT / "results" / fn).is_file()]
    print(json.dumps({"bootstrap_replicates": 5000, "bootstrap_seed": 20260914,
                      "ratio_definition": "sqrt(median numerator trace / median denominator trace), common valid trials",
                      "datasets": reports}, indent=2, allow_nan=False))
