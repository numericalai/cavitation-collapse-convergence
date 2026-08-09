#!/usr/bin/env python3
"""
MFC case -- 3D spherical bubble collapse, first-octant symmetry.

Built directly on the shipped 3D_sphbubcollapse tutorial so that Phase 0
reproduces it exactly (SB_DRIVING=1 SB_PRATIO=100 gives the shipped numbers).
The parametrisation exposes the one axis that matters for choke-trim
cavitation: the driving pressure, i.e. the differential across the trim.

Nondimensionalisation follows the template: lengths are scaled by L0 = 1 mm,
so R0 = 1 in code units. Pressures are in SI. Post-processing multiplies code
lengths by L0 and code times by the template's time unit before any comparison
to the SI Rayleigh prediction in thermo.py.

Run:
    mfc.sh run case.py -t pre_process simulation
    SB_PDRIVE=100e5 mfc.sh run case.py -t pre_process simulation
    SB_SUMMARY=1 python3 case.py

Apache License 2.0.
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo import (WATER, GAS, RHO_WATER, P_VAPOR,
                    rayleigh_collapse_time, rayleigh_plesset)

L0 = 1.0e-3                       # length scale [m]; R0 = 1 in code units


def _env(name, default, cast=float):
    return cast(os.environ.get(name, default))


# --- parameters -------------------------------------------------------------
# The driving pressure is THE sweep axis. The shipped template uses p_out = 1 bar
# with a 100x interior deficit; real choke trim cavitates at tens to hundreds of
# bar differential. SB_PDRIVE sets the external (liquid) pressure directly.
P_DRIVE = _env("SB_PDRIVE", 1.0e5)          # external liquid pressure [Pa]
P_RATIO = _env("SB_PRATIO", 100.0)          # p_out / p_bubble
P_BUBBLE = P_DRIVE / P_RATIO                 # bubble interior pressure [Pa]

NX = _env("SB_NX", 99, int)                  # cells per octant edge (m=n=p=NX)
R0_PHYS = _env("SB_R0", 1.0e-3)              # initial bubble radius [m]
DOMAIN = _env("SB_DOMAIN", 4.0)              # domain edge in R0 units
CFL_DT = _env("SB_DT", 0.0)                 # timestep [s]; 0 = auto-size below
N_SAVE = _env("SB_NSAVE", 60, int)   # snapshots over the run; ~30 needed for
                                     # the convergence compare, 60 is ample.

R0 = R0_PHYS / L0                            # = 1.0 by default

# Auto-size the timestep when SB_DT is not given. The reframed runs stop BEFORE
# the focus, so the relevant sound speed is the early-collapse liquid value
# (~1700 m/s), not the enormous focus value. That permits a far larger timestep
# than the template's fixed 2e-10 and is what makes these runs cheap. A run that
# deliberately approaches the focus must set SB_DT explicitly.
if CFL_DT <= 0.0:
    c_early = WATER.c(P_DRIVE, RHO_WATER)   # liquid sound speed at driving pres.
    dx_min_phys = (DOMAIN * R0) / NX / 4.0 * L0   # smallest cell under 4x stretch
    CFL_TARGET = _env("SB_CFL", 0.4)
    CFL_DT = CFL_TARGET * dx_min_phys / c_early

dt_code = CFL_DT / L0                        # template scales dt by 1/L0 too

# End time. The reframed scope stops the run BELOW the focus, where the diffuse-
# interface collapse is stable and Rayleigh-Plesset is valid. The focus (R -> 0,
# ~20x compression) is where the six-equation model goes unstable and is
# explicitly OUT OF SCOPE. Default 0.6 t_c: compression only ~1.2x, safely below
# the observed ~0.66 t_c abort. Raise SB_ENDFRAC deliberately to map the
# stability ceiling, never blindly.
END_FRAC = _env("SB_ENDFRAC", 0.6)
tc = rayleigh_collapse_time(R0_PHYS, P_DRIVE, RHO_WATER, P_VAPOR)
t_end_phys = END_FRAC * tc
n_steps = int(math.ceil(t_end_phys / CFL_DT))
t_save = max(1, n_steps // N_SAVE)

# Timestep must resolve the acoustic CFL near minimum radius, where the liquid
# is most compressed and the sound speed highest. The template's fixed dt was
# tuned for its 1-bar case; higher driving pressure needs a check. The stiffest
# point is the peak gas pressure at minimum radius.
_rp_peak = rayleigh_plesset(R0_PHYS, P_DRIVE, P_BUBBLE,
                            rho=RHO_WATER, p_v=P_VAPOR)["peak_gas_pressure"]
c_max = WATER.c(_rp_peak, RHO_WATER)
# Smallest cell is at the origin under 4x stretching; estimate it.
dx_min = (DOMAIN * R0) / NX / 4.0
cfl_est = c_max * CFL_DT / (dx_min * L0)

# The 1-bar case genuinely needs ~700k steps to fully collapse, which is why the
# template stops early. For the choke sweep (high driving pressure) the collapse
# is far faster and the step count is modest. Warn rather than silently produce
# a runaway run.
if n_steps > 400000:
    sys.stderr.write(
        f"WARNING: {n_steps} steps to reach {END_FRAC:g} t_c at "
        f"{P_DRIVE/1e5:g} bar. Low driving pressure means slow collapse. "
        "The choke sweep (>= 10 bar) is far cheaper; for a full-collapse run "
        "at 1 bar, expect this cost or raise SB_DT.\n")
if cfl_est > 0.8:
    sys.stderr.write(
        f"WARNING: estimated acoustic CFL {cfl_est:.2f} at minimum radius "
        f"(c_max ~ {c_max:.0f} m/s). Reduce SB_DT for the high-pressure runs.\n")

# --- densities from the template's convention -------------------------------
# The template sets alpha and alpha_rho directly rather than from an EOS. We
# keep its values so Phase 0 is bit-for-bit the tutorial, and only scale the
# pressures. (alpha_rho is partial density = alpha * rho_phase.)
LIQ_ALPHA_RHO_W = 1000.0
LIQ_ALPHA_RHO_G = 0.1
LIQ_ALPHA_W = 0.9
LIQ_ALPHA_G = 0.1

BUB_ALPHA_RHO_W = 100.0
BUB_ALPHA_RHO_G = 0.9
BUB_ALPHA_W = 0.1
BUB_ALPHA_G = 0.9

STR = DOMAIN                                 # domain end in code units

case = {
    "run_time_info": "T",

    "x_domain%beg": 0.0, "x_domain%end": STR,
    "y_domain%beg": 0.0, "y_domain%end": STR,
    "z_domain%beg": 0.0, "z_domain%end": STR,

    # Grid stretching, concentrating cells near the origin (the collapse point).
    "stretch_x": "T", "a_x": 4.0, "x_a": -1.5 * R0, "x_b": 1.5 * R0,
    "stretch_y": "T", "a_y": 4.0, "y_a": -1.5 * R0, "y_b": 1.5 * R0,
    "stretch_z": "T", "a_z": 4.0, "z_a": -1.5 * R0, "z_b": 1.5 * R0,
    "cyl_coord": "F",

    "m": NX, "n": NX, "p": NX,
    "dt": dt_code,
    "t_step_start": 0,
    "t_step_stop": n_steps,
    "t_step_save": t_save,

    "num_patches": 2,
    "model_eqns": 3,
    "alt_soundspeed": "F",
    "num_fluids": 2,
    "mpp_lim": "T",
    "mixture_err": "T",

    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "weno_Re_flux": "F",
    "weno_avg": "F",
    "avg_state": 2,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": 2,
    "wave_speeds": 1,

    # Reflecting on the octant symmetry planes, non-reflecting on the far faces.
    "bc_x%beg": -2, "bc_x%end": -6,
    "bc_y%beg": -2, "bc_y%end": -6,
    "bc_z%beg": -2, "bc_z%end": -6,

    "format": 1,
    "precision": 2,
    # OUTPUT REDUCTION. The analysis reads only prim.2 (gas fraction). Writing
    # conservative variables as well doubles the output for nothing, so switch
    # them off. Combined with a lower SB_NSAVE this keeps even the 149^3 run to
    # a fraction of what the 99^3 produced. prim_vars_wrt stays on because
    # prim.2 is the one variable we need.
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "F",
    "parallel_io": "F",

    # Patch 1: high-pressure liquid background (cube).
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 80.0 * R0,
    "patch_icpp(1)%y_centroid": 80.0 * R0,
    "patch_icpp(1)%z_centroid": 80.0 * R0,
    "patch_icpp(1)%length_x": 160.0 * R0,
    "patch_icpp(1)%length_y": 160.0 * R0,
    "patch_icpp(1)%length_z": 160.0 * R0,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": P_DRIVE,
    "patch_icpp(1)%alpha_rho(1)": LIQ_ALPHA_RHO_W,
    "patch_icpp(1)%alpha_rho(2)": LIQ_ALPHA_RHO_G,
    "patch_icpp(1)%alpha(1)": LIQ_ALPHA_W,
    "patch_icpp(1)%alpha(2)": LIQ_ALPHA_G,

    # Patch 2: spherical gas bubble at the origin.
    "patch_icpp(2)%geometry": 8,
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.5,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%z_centroid": 0.0,
    "patch_icpp(2)%radius": R0,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%vel(3)": 0.0,
    "patch_icpp(2)%pres": P_BUBBLE,
    "patch_icpp(2)%alpha_rho(1)": BUB_ALPHA_RHO_W,
    "patch_icpp(2)%alpha_rho(2)": BUB_ALPHA_RHO_G,
    "patch_icpp(2)%alpha(1)": BUB_ALPHA_W,
    "patch_icpp(2)%alpha(2)": BUB_ALPHA_G,

    "fluid_pp(1)%gamma": WATER.mfc_gamma,
    "fluid_pp(1)%pi_inf": WATER.mfc_pi_inf,
    "fluid_pp(2)%gamma": GAS.mfc_gamma,
    "fluid_pp(2)%pi_inf": GAS.mfc_pi_inf,
}


def summary():
    rp = rayleigh_plesset(R0_PHYS, P_DRIVE, P_BUBBLE, rho=RHO_WATER, p_v=P_VAPOR)
    lines = [
        "3D bubble collapse (octant) -- case summary",
        "-" * 60,
        f"  driving pressure p_out   {P_DRIVE/1e5:>12.4g} bar",
        f"  bubble pressure p_in     {P_BUBBLE/1e5:>12.4g} bar   "
        f"(ratio {P_RATIO:g})",
        f"  initial radius R0        {R0_PHYS*1e3:>12.4g} mm",
        f"  cells per octant edge    {NX:>12d}   ({NX**3:,} total)",
        f"  timestep                 {CFL_DT:>12.4g} s",
        f"  steps                    {n_steps:>12d}   (save every {t_save})",
        "",
        "  ANALYTICAL ANCHORS (SI, from thermo.py)",
        f"    Rayleigh collapse time t_c   {tc*1e6:>10.4f} us",
        f"    R-P minimum radius           {rp['R_min']/R0_PHYS:>10.4f} R0",
        f"    R-P compression ratio        {rp['compression_ratio']:>10.1f}",
        f"    R-P peak gas pressure        {rp['peak_gas_pressure']/1e5:>10.1f} bar",
        "",
        f"  run reaches {END_FRAC:.1f} t_c = {t_end_phys*1e6:.3f} us",
        "",
        f"  peak sound speed est.    {c_max:>10.0f} m/s   "
        f"(acoustic CFL ~ {cfl_est:.2f})",
        f"  smallest cell (origin)   {dx_min*L0*1e6:>10.2f} um   "
        f"({(R0_PHYS)/(dx_min*L0):.1f} cells per R0)",
    ]
    if abs(P_DRIVE - 1e5) < 1 and abs(P_RATIO - 100) < 1:
        lines.append("  (this is the SHIPPED template configuration -- Phase 0)")
    return "\n".join(lines)


if __name__ == "__main__":
    if os.environ.get("SB_SUMMARY"):
        print(summary())
    else:
        print(json.dumps(case))