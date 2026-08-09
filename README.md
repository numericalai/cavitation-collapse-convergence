# Cavitation bubble collapse: a grid-convergence study

**Does a 3D cavitation bubble collapse converge on a uniform grid? No — and the
finest grid is the one that fails first.**

This repository documents an attempt to validate a spherical cavitation-bubble
collapse against analytical theory (Rayleigh–Plesset) using the open-source
compressible multiphase solver [MFC](https://github.com/MFlowCode/MFC), as a
step toward simulating cavitation erosion. The collapse **does not
grid-converge** on a uniform octant mesh — not at 60k cells, not at 970k, and
not at 3.3 million — and refining the grid makes the run fail *sooner*, not
later. Everything needed to reproduce that finding is here.

It's a negative result, published because it's useful: if you're reaching for a
fixed-grid diffuse-interface method to resolve a bubble collapse to its focus,
this is the evidence that it won't get there, and a map of exactly where it
stops.

## The finding

Three resolutions (39³, 99³, 149³) produce three different collapse curves. The
99³ and 149³ start within ~1% of each other but diverge to ~10% through the
steep mid-collapse — refining does **not** bring them inside a usable band. And
the finest grid aborts at 0.35 t_c, *earlier* than the 99³ at 0.56 t_c, because
resolving the collapse more sharply makes the interface instability arrive
sooner.

The one thing that *does* converge is the static initial sphere
(`figures/fig2_r0_converge.png`) — 7.9% radius error at 39³, 2.6% at 99³, 2.2%
at 149³ — which localises the failure to the collapse dynamics, not the
discretisation. The method is sound; it's the moving, self-concentrating
collapse a fixed grid can't keep up with.

## Results at a glance

| grid | cells | initial-radius error | reached | outcome |
|------|-------|----------------------|---------|---------|
| 39³ | 59,319 | 7.9% | 0.50 t_c | completed (too dissipative to go unstable) |
| 99³ | 970,299 | 2.6% | 0.56 t_c | aborts near focus |
| 149³ | 3,307,949 | 2.2% | 0.35 t_c | **aborts earliest** |

Pairwise collapse-curve difference over the common window: 39³-vs-99³ ~10% max,
99³-vs-149³ ~10% max (spiking through the steep collapse). Refinement does not
shrink it below the ~5% you would want for a converged result.

## Repository layout

```
code/
  thermo.py     analytical anchors — Rayleigh collapse time, Rayleigh–Plesset,
                stiffened-gas EOS. `python code/thermo.py` runs 7 self-checks.
  case.py       parametrised MFC case (first-octant symmetry, driven collapse).
                Env vars: SB_NX, SB_PDRIVE, SB_ENDFRAC, SB_NSAVE.
  post.py       radius extraction (resolution-independent alpha=0.5 profile
                method), two- and three-grid comparison, CSV export.
  reap_cons.sh  disk janitor — deletes redundant conservative-variable output
                during a run, for space-constrained machines.
results/
  R_of_t_39.csv, R_of_t_99.csv, R_of_t_149.csv   extracted R(t) per grid
  convergence_summary.csv                        the headline numbers
  DATA_ARCHIVE.md                                raw field data: format + link
figures/
  fig1_nonconvergence.png   three grids diverge; refinement doesn't help
  fig2_r0_converge.png      the static sphere DOES converge (method is sound)
  fig3_focus_wall.png       the focus is unreachable; finer grids abort sooner
post/
  blog_negative_result.md   the write-up
  EDITORIAL_NOTES.md        framing, reviewer objections, and how they're met
```

The raw field data (gas-fraction fields for all three grids, ~29 GB) is too
large for git and is archived separately — see
[`results/DATA_ARCHIVE.md`](results/DATA_ARCHIVE.md) for the download and format.

## Reproduce it

**1. Verify the analytical anchor** (no CFD, runs in seconds):

```bash
python code/thermo.py
```

**2. Run the three grids** (GPU; dual-GPU for the two fine grids — see the post
for machine notes). The runs stop below the collapse focus, where the method is
stable:

```bash
SB_NX=39  SB_PDRIVE=100e5 SB_ENDFRAC=0.5 ./mfc.sh run code/case.py -t pre_process simulation
SB_NX=99  SB_PDRIVE=100e5 SB_ENDFRAC=0.5 ./mfc.sh run code/case.py -t pre_process simulation -n 2
SB_NX=149 SB_PDRIVE=100e5 SB_ENDFRAC=0.5 ./mfc.sh run code/case.py -t pre_process simulation -n 2
```

**3. Extract R(t) and compare** — this is the finding:

```bash
python code/post.py --compare3 D_39 D_99 D_149 --alpha-var 2
```

The extracted curves in `results/` reproduce `figures/fig1`. If you have only
the archived data (not a full MFC install), steps 1 and 3 still run against the
downloaded `D_*` directories.

## What this is and isn't

- **Is:** a quantified boundary — where a fixed-grid diffuse-interface collapse
  stops being trustworthy — backed by three grids and an exact analytical
  anchor.
- **Isn't:** a cavitation-erosion prediction. The collapse focus (peak pressure,
  microjet) is unreachable in this configuration. Reaching it needs adaptive
  mesh refinement or interface sharpening — noted as the path forward, not done
  here.

## Solver

Built on [MFC](https://github.com/MFlowCode/MFC), an open-source compressible
multiphase flow solver. The case derives from its shipped `3D_sphbubcollapse`
example, re-parametrised for a driven collapse at choke-relevant pressures.

## License

Code: MIT. Data and figures: CC-BY. If you use this, please cite the archive
DOI (see `results/DATA_ARCHIVE.md`) and link back to this repository.