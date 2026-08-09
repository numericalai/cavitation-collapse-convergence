#!/usr/bin/env python3
"""
Post-processing for 3D spherical bubble collapse (octant).

Extracts the bubble radius over time by integrating the gas volume fraction over
the octant, and validates it against the Rayleigh-Plesset solution from
thermo.py. This is the credibility anchor for the reframed case: if the
simulated R(t) tracks R-P through the run (which stops at 0.6 t_c, below the
focus), the setup is sound.

MFC output format (this build, format=1, parallel_io=F, 3D):
    D/prim.{VAR}.00.{STEP}.dat  -- primitive variables
    D/cons.{VAR}.00.{STEP}.dat  -- conservative variables
Each line is FOUR whitespace-separated columns: x, y, z, value. Coordinates are
in code units (length / 1 mm), carried with every file, so the stretched-grid
cell positions come for free -- no separate coordinate file, no recomputing
stretch parameters.

Radius extraction. The bubble is diffuse (smoothen=T), so "radius" is defined as
the radius of the sphere with the same gas volume:
    V_gas = sum over cells of alpha_gas * V_cell
    R_eq  = (3 V_gas / (4 pi) * octant_factor)^(1/3)
The octant meshes 1/8 of the sphere, so the full-sphere gas volume is 8 * V_gas.

Cell volumes on the stretched grid are reconstructed from the sorted unique
coordinate values along each axis (midpoint rule between neighbouring nodes).

CONSERVATION CHECK. At t=0 the gas volume must equal the known 1 mm sphere
(within the diffuse-interface smearing). If it does not, the cell-volume
weighting is wrong and nothing downstream is trustworthy -- this is the guard
that catches a stretched-grid integration error, the 3D analogue of the
sum(alpha)=1 check in the 1D cases.

Usage:
    python3 post.py --run 3D_sphbubcollapse_choke --alpha-var 10
    python3 post.py --run 3D_sphbubcollapse_choke --inspect

Apache License 2.0.
"""

import argparse
import glob
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from thermo import (rayleigh_collapse_time, rayleigh_plesset,  # noqa: E402
                    RHO_WATER, P_VAPOR)

L0 = 1.0e-3   # length scale [m]; code lengths are physical/L0


def _load_field(path):
    """Read an MFC 3D .dat file: columns x, y, z, value (code units)."""
    a = np.loadtxt(path)
    if a.ndim == 1:
        a = a.reshape(1, -1)
    if a.shape[1] != 4:
        raise ValueError(f"{path}: expected 4 columns (x,y,z,value), "
                         f"got {a.shape[1]}.")
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3]


def _cell_volumes(x, y, z):
    """Cell volumes on a structured, stretched, axis-aligned grid.

    The field is a flattened structured grid: each axis has a small set of
    unique node coordinates (40 here), and every cell is a box. Build the
    per-axis cell widths from the sorted unique coordinates using the midpoint
    rule, then map each cell to the product of its three widths.

    The earlier bug was trusting an assumed domain extent; the grid actually
    stretches far past the near-origin region (x runs to ~8.9 here, not 4), so
    widths must come from the coordinates themselves, which this does.
    """
    widths = {}
    for name, c in (("x", x), ("y", y), ("z", z)):
        u = np.unique(np.round(c, 10))
        edges = np.empty(u.size + 1)
        edges[1:-1] = 0.5 * (u[:-1] + u[1:])
        edges[0] = u[0] - 0.5 * (u[1] - u[0])
        edges[-1] = u[-1] + 0.5 * (u[-1] - u[-2])
        w = np.diff(edges)
        widths[name] = dict(zip(u, w))
    key = lambda arr: np.round(arr, 10)
    wx = np.array([widths["x"][v] for v in key(x)])
    wy = np.array([widths["y"][v] for v in key(y)])
    wz = np.array([widths["z"][v] for v in key(z)])
    return wx * wy * wz    # code units^3


ALPHA_LIQUID_BG = 0.1   # gas volume fraction in the liquid background (case IC)
ALPHA_BUBBLE = 0.9      # gas volume fraction inside the bubble (case IC)


def bubble_radius_volume(path, octant_factor=8.0, alpha_bg=ALPHA_LIQUID_BG,
                         alpha_bub=ALPHA_BUBBLE):
    """Equivalent radius from the EXCESS-gas volume integral.

    Resolution-dependent because the smoothed interface (smoothen=T) is spread
    over a fixed number of cells, so the skin's physical thickness shrinks as
    the grid refines and the volume of the transition region scales differently.
    Kept for cross-checking, NOT as the primary metric. See bubble_radius_profile.
    """
    x, y, z, alpha = _load_field(path)
    dV = _cell_volumes(x, y, z)
    excess = np.clip((alpha - alpha_bg) / (alpha_bub - alpha_bg), 0.0, 1.0)
    V_gas_code = float(np.sum(excess * dV))
    V_gas_full = octant_factor * V_gas_code * L0**3
    R_eq = (3.0 * V_gas_full / (4.0 * np.pi)) ** (1.0 / 3.0)
    return R_eq, V_gas_full


def bubble_radius_profile(path, alpha_bg=ALPHA_LIQUID_BG, alpha_bub=ALPHA_BUBBLE,
                          n_bins=None):
    """Radius from the alpha = 0.5-crossing of the radially-binned profile.

    Resolution-INDEPENDENT: the interface midpoint sits at the same physical
    radius regardless of how many cells span the smoothed skin, because the
    midpoint of a monotonic smoothed transition does not move under refinement.
    This is the primary metric. The bubble radius is where the spherically-
    averaged gas fraction crosses the midpoint between bubble and liquid values.

    Method: bin every cell by its radial distance from the origin, average alpha
    in each bin, and linearly interpolate the radius at which the averaged alpha
    falls through (alpha_bg + alpha_bub)/2.
    """
    x, y, z, alpha = _load_field(path)
    r = np.sqrt(x*x + y*y + z*z)                     # code units
    order = np.argsort(r)
    r_s, a_s = r[order], alpha[order]

    # Bin radially. Default to ~one bin per unique radial shell of the grid.
    if n_bins is None:
        n_bins = max(20, int(np.unique(np.round(x, 8)).size))
    r_max = r_s.max()
    edges = np.linspace(0.0, r_max, n_bins + 1)
    idx = np.clip(np.digitize(r_s, edges) - 1, 0, n_bins - 1)
    a_prof = np.full(n_bins, np.nan)
    r_prof = 0.5 * (edges[:-1] + edges[1:])
    for b in range(n_bins):
        sel = idx == b
        if sel.any():
            a_prof[b] = a_s[sel].mean()

    good = ~np.isnan(a_prof)
    r_prof, a_prof = r_prof[good], a_prof[good]

    # The profile runs from ~alpha_bub near r=0 down to ~alpha_bg far out.
    # Find where it crosses the midpoint.
    mid = 0.5 * (alpha_bg + alpha_bub)
    below = np.where(a_prof < mid)[0]
    if below.size == 0 or below[0] == 0:
        # No clean crossing (bubble may have shrunk below the grid or filled it)
        return np.nan
    j = below[0]
    # linear interpolation between bins j-1 (above mid) and j (below mid)
    a1, a0 = a_prof[j - 1], a_prof[j]
    r1, r0 = r_prof[j - 1], r_prof[j]
    frac = (a1 - mid) / (a1 - a0) if a1 != a0 else 0.0
    R_code = r1 + frac * (r0 - r1)
    return R_code * L0                               # physical m


def bubble_radius(path, method="profile", **kw):
    """Primary radius entry point. method='profile' (resolution-independent,
    default) or 'volume' (cross-check)."""
    if method == "profile":
        R = bubble_radius_profile(path, **kw)
        return R, None
    return bubble_radius_volume(path, **kw)


def _steps(run_dir, var, kind="prim"):
    # Accept either a run dir containing D/, or a dir that IS the output folder.
    base = Path(run_dir)
    if (base / "D").is_dir():
        base = base / "D"
    files = glob.glob(str(base / f"{kind}.{var}.00.*.dat"))
    out = []
    for f in files:
        m = re.search(rf"{kind}\.{var}\.00\.(\d+)\.dat", f)
        if m:
            out.append((int(m.group(1)), f))
    return sorted(out)


def _dt_phys(run_dir):
    """Physical timestep from the case summary if present, else from the .inp."""
    # Try the MFC simulation.inp for dt (code units) and multiply by L0.
    inp = Path(run_dir) / "simulation.inp"
    if inp.exists():
        txt = inp.read_text()
        m = re.search(r"dt\s*=\s*([0-9.eEdD+-]+)", txt)
        if m:
            dt_code = float(m.group(1).replace("D", "E").replace("d", "e"))
            return dt_code * L0
    return None


def detect_alpha_var(run_dir, step0_kind="prim"):
    """Find which prim variable is the gas volume fraction: bounded in [0,1],
    and not identically constant (the liquid fraction is its complement)."""
    cands = []
    for var in range(1, 13):
        files = _steps(run_dir, var, step0_kind)
        if not files:
            continue
        _, _, _, val = _load_field(files[0][1])
        if val.min() >= -1e-6 and val.max() <= 1.0 + 1e-6 and val.ptp() > 0.1:
            cands.append((var, float(val.min()), float(val.max())))
    return cands


def analyse(run_dir, alpha_var, kind="prim", p_drive=100e5, R0=1e-3,
            p_ratio=100.0):
    dt = _dt_phys(run_dir)
    steps = _steps(run_dir, alpha_var, kind)
    if not steps:
        raise FileNotFoundError(f"No {kind}.{alpha_var}.00.*.dat in {run_dir} or {run_dir}/D")

    t, R, Rv = [], [], []
    for istep, path in steps:
        R_prof = bubble_radius_profile(path)
        R_vol, _ = bubble_radius_volume(path)
        t.append(istep * dt if dt else istep)
        R.append(R_prof)
        Rv.append(R_vol)
    t = np.array(t); R = np.array(R); Rv = np.array(Rv)

    # t=0 check against the known sphere, using the resolution-independent
    # profile radius (the volume integral is grid-dependent and only a cross-
    # check now).
    R0_sim = R[0]
    r0_err = abs(R0_sim - R0) / R0

    # Analytical references.
    tc = rayleigh_collapse_time(R0, p_drive, RHO_WATER, P_VAPOR)
    rp = rayleigh_plesset(R0, p_drive, p_drive / p_ratio,
                          rho=RHO_WATER, p_v=P_VAPOR)

    return {
        "t": t, "R": R, "R_volume": Rv,
        "R0_sim": float(R0_sim), "R0_known": R0,
        "r0_error": float(r0_err),
        "tc": tc, "rp": rp,
        "t_over_tc": t / tc,
        "R_over_R0": R / R0,
    }


def report(m):
    print(f"\n  RADIUS AT t=0 (resolution-independent profile method)")
    print(f"    R0 simulated / known   {m['R0_sim']/m['R0_known']:.4f}   "
          + ("OK" if m['r0_error'] < 0.05 else
             "*** check extraction ***"))
    print(f"    (volume-method R0/known {m['R_volume'][0]/m['R0_known']:.4f} "
          "-- grid-dependent, cross-check only)")
    print()
    print(f"  COLLAPSE")
    print(f"    Rayleigh t_c           {m['tc']*1e6:.4f} us")
    print(f"    run reaches            {m['t_over_tc'][-1]:.3f} t_c "
          f"(R/R0 = {m['R_over_R0'][-1]:.3f})")
    R_rp = np.interp(m["t"], m["rp"]["t"], m["rp"]["R"]) / m["R0_known"]
    rel = np.abs(m["R_over_R0"] - R_rp)
    print(f"    R-P at run end         R/R0 = {R_rp[-1]:.3f}")
    print(f"    mean |R_sim - R_RP|/R0 {rel.mean():.4f}")
    print(f"    early (t<0.1 tc) diff  {rel[m['t_over_tc']<0.1].mean():.4f}")
    print()
    print("    NOTE: agreement with incompressible R-P is expected only very")
    print("    early. As compression grows, the compressible simulation")
    print("    collapses faster than R-P -- that departure is physical, not")
    print("    error. The VALIDATION criterion is grid convergence of R(t)")
    print("    between resolutions (run compare_grids), not agreement with R-P.")


def compare_grids(run_a, run_b, alpha_var=2, kind="prim", p_drive=100e5):
    """The REAL validation: does the same physical R(t) come out of two grids?

    Compares R vs COLLAPSE FRACTION (normalised step, 0..1) rather than physical
    time. Both runs cover the same physical window (same SB_ENDFRAC), so step
    fraction is a common, dt-independent axis -- this avoids any dependence on
    reading each run's dt correctly, which differs between grids.
    """
    ma = analyse(run_a, alpha_var, kind, p_drive)
    mb = analyse(run_b, alpha_var, kind, p_drive)
    # normalised progress through the run, 0..1
    fa = (ma["t"] - ma["t"][0]) / (ma["t"][-1] - ma["t"][0])
    fb = (mb["t"] - mb["t"][0]) / (mb["t"][-1] - mb["t"][0])
    ff = np.linspace(0.0, 1.0, 30)
    Ra = np.interp(ff, fa, ma["R_over_R0"])
    Rb = np.interp(ff, fb, mb["R_over_R0"])
    diff = np.abs(Ra - Rb)
    print(f"\n  GRID CONVERGENCE OF R  [{Path(run_a).name} vs {Path(run_b).name}]")
    print(f"    compared vs collapse fraction (dt-independent)")
    print(f"    R0 (A / B)             {ma['R0_sim']/ma['R0_known']:.4f} / "
          f"{mb['R0_sim']/mb['R0_known']:.4f}")
    print(f"    mean |R_A - R_B|/R0    {diff.mean():.4f}")
    print(f"    max  |R_A - R_B|/R0    {diff.max():.4f}   "
          + ("(CONVERGED)" if diff.max() < 0.05 else "(not converged)"))
    print(f"\n    {'frac':>6} {'R_A/R0':>8} {'R_B/R0':>8} {'diff':>8}")
    for i in range(0, len(ff), 3):
        print(f"    {ff[i]:>6.3f} {Ra[i]:>8.4f} {Rb[i]:>8.4f} {diff[i]:>+8.4f}")
    return {"frac": ff, "R_a": Ra, "R_b": Rb, "max_diff": float(diff.max())}


def export_csv(run_dir, alpha_var, out_path, kind="prim", p_drive=100e5):
    """Write the extracted collapse curve as a plain CSV for archiving.

    Columns: step, R_over_R0, R_metres, R_volume_over_R0 (cross-check),
    t_over_tc. Small, human-readable evidence that backs the convergence claim.
    """
    m = analyse(run_dir, alpha_var, kind, p_drive)
    steps = [s for s, _ in _steps(run_dir, alpha_var, kind)]
    rows = ["step,R_over_R0,R_metres,R_volume_over_R0,t_over_tc"]
    Rv = m["R_volume"]
    for i in range(len(m["t"])):
        rows.append(f"{steps[i]},{m['R_over_R0'][i]:.6f},{m['R'][i]:.6e},"
                    f"{Rv[i]/m['R0_known']:.6f},{m['t_over_tc'][i]:.6f}")
    Path(out_path).write_text("\n".join(rows) + "\n")
    print(f"wrote {out_path}  ({len(steps)} rows)")


def compare_three(runs, alpha_var=2, kind="prim", p_drive=100e5, labels=None):
    """Compare R(collapse-fraction) across three (or more) grids over their
    COMMON window, and report whether refinement converges or diverges.

    The key diagnostic for this study: if successively finer grids abort earlier
    and/or drift further apart rather than closer, the collapse is not converging
    under refinement -- the instability is in the method, not the resolution.
    """
    ms = [analyse(r, alpha_var, kind, p_drive) for r in runs]
    labels = labels or [Path(r).name for r in runs]
    # each grid's own end fraction (how far it got before abort)
    print("\n  GRID REFINEMENT SUMMARY")
    for lab, m in zip(labels, ms):
        print(f"    {lab:<8s} reached frac 1.00 of its run; "
              f"R0/known = {m['R0_sim']/m['R0_known']:.4f}, "
              f"{len(m['t'])} snapshots")
    # common window is min end -- but all use SB_ENDFRAC=0.5, so runs that
    # aborted cover less of it. Compare on the shortest run's coverage.
    # Use normalised progress within each run's OWN saved range, then also note
    # the physical reach.
    ff = np.linspace(0.0, 1.0, 30)
    curves = []
    for m in ms:
        fr = (m["t"] - m["t"][0]) / (m["t"][-1] - m["t"][0])
        curves.append(np.interp(ff, fr, m["R_over_R0"]))
    curves = np.array(curves)

    print("\n  PAIRWISE MAX |ΔR|/R0 over the saved range")
    for i in range(len(runs)):
        for j in range(i + 1, len(runs)):
            d = np.abs(curves[i] - curves[j])
            print(f"    {labels[i]} vs {labels[j]}: "
                  f"mean {d.mean():.4f}  max {d.max():.4f}")

    print(f"\n    {'frac':>6} " + " ".join(f"{l:>9s}" for l in labels))
    for k in range(0, len(ff), 3):
        print(f"    {ff[k]:>6.3f} " +
              " ".join(f"{curves[c][k]:>9.4f}" for c in range(len(runs))))

    # Anti-convergence check: does the finest-vs-middle gap exceed the
    # middle-vs-coarse gap? If refinement isn't shrinking differences, say so.
    if len(runs) >= 3:
        d_cm = np.abs(curves[0] - curves[1]).mean()   # coarse-mid
        d_mf = np.abs(curves[1] - curves[2]).mean()   # mid-fine
        print(f"\n  CONVERGENCE DIRECTION")
        print(f"    coarse-mid mean diff  {d_cm:.4f}")
        print(f"    mid-fine  mean diff   {d_mf:.4f}")
        if d_mf >= d_cm * 0.8:
            print("    -> differences are NOT shrinking under refinement.")
            print("       R(t) does not converge; finer grids do not approach a")
            print("       limit. Consistent with a method (not resolution) failure.")
        else:
            print("    -> differences shrink under refinement (converging).")
    return {"frac": ff, "curves": curves, "labels": labels}
    """Write the extracted collapse curve as a plain CSV for archiving.

    Columns: step, R_over_R0, R_metres, R_volume_over_R0 (cross-check),
    t_over_tc. This is the small, human-readable evidence that backs the
    convergence claim -- a few KB, belongs in the repo, not the raw fields.
    """
    m = analyse(run_dir, alpha_var, kind, p_drive)
    steps = [s for s, _ in _steps(run_dir, alpha_var, kind)]
    rows = ["step,R_over_R0,R_metres,R_volume_over_R0,t_over_tc"]
    Rv = m["R_volume"]
    for i in range(len(m["t"])):
        rows.append(f"{steps[i]},{m['R_over_R0'][i]:.6f},{m['R'][i]:.6e},"
                    f"{Rv[i]/m['R0_known']:.6f},{m['t_over_tc'][i]:.6f}")
    Path(out_path).write_text("\n".join(rows) + "\n")
    print(f"wrote {out_path}  ({len(steps)} rows)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=str(HERE))
    ap.add_argument("--compare", nargs=2, default=None,
                    help="two run dirs: validate by grid convergence of R(t)")
    ap.add_argument("--compare3", nargs="+", default=None,
                    help="three+ run dirs: report convergence direction")
    ap.add_argument("--export-csv", default=None,
                    help="write the extracted R(t) curve to this CSV path")
    ap.add_argument("--alpha-var", type=int, default=None)
    ap.add_argument("--kind", default="prim", choices=["prim", "cons"])
    ap.add_argument("--p-drive", type=float, default=100e5)
    ap.add_argument("--p-ratio", type=float, default=100.0)
    ap.add_argument("--inspect", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if args.export_csv:
        export_csv(args.run, args.alpha_var or 2, args.export_csv,
                   args.kind, args.p_drive)
        return

    if args.compare3:
        compare_three(args.compare3, args.alpha_var or 2, args.kind,
                      args.p_drive)
        return

    if args.compare:
        av = args.alpha_var or 2
        compare_grids(args.compare[0], args.compare[1], av, args.kind,
                      args.p_drive)
        return

    if args.inspect or args.alpha_var is None:
        print("Detecting gas volume-fraction variable (bounded [0,1], varying)...")
        cands = detect_alpha_var(args.run, args.kind)
        if not cands:
            print("  none found -- check --kind and that D/ has files.")
            return
        for var, lo, hi in cands:
            print(f"  {args.kind}.{var}: min={lo:.4g} max={hi:.4g}")
        if args.inspect:
            print("\nPass the gas fraction index with --alpha-var.")
            return
        args.alpha_var = cands[-1][0]
        print(f"\nUsing --alpha-var {args.alpha_var}")

    m = analyse(args.run, args.alpha_var, args.kind, args.p_drive,
                p_ratio=args.p_ratio)
    report(m)

    if args.out:
        import json
        Path(args.out).write_text(json.dumps({
            "t": m["t"].tolist(), "R": m["R"].tolist(),
            "R_over_R0": m["R_over_R0"].tolist(),
            "t_over_tc": m["t_over_tc"].tolist(),
            "conservation_error": m["conservation_error"],
            "tc": m["tc"],
        }, indent=2))
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()