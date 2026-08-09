# Data archive — 3D spherical bubble collapse, grid non-convergence study

Raw field data and extracted curves backing the post
*"When the Bubble Won't Sit Still: A Negative Result on Cavitation Collapse."*

This archive exists so the central claim — that the collapse does **not**
grid-converge across 39³, 99³, and 149³, and that the finest grid fails
*earliest* — is independently reproducible, including re-doing the radius
extraction with a different method if you wish.

## What's here

```
D_39/        gas-fraction field output (prim.2 only), 39³  octant  (60k cells)
D_99/        gas-fraction field output (prim.2 only), 99³  octant  (970k cells)
D_149/       gas-fraction field output (prim.2 only), 149³ octant  (3.3M cells)
results/
  R_of_t_39.csv    extracted bubble radius vs time, 39³
  R_of_t_99.csv    extracted bubble radius vs time, 99³
  R_of_t_149.csv   extracted bubble radius vs time, 149³
  convergence_summary.txt   the three-grid comparison output
figures/           the three figures from the post
code/              case.py, thermo.py, post.py (also on GitHub)
```

**Note on the reduced field set.** The full MFC runs wrote all conservative and
primitive variables at every saved step (~686 GB for the 99³ alone). This
archive contains **only `prim.2`, the gas volume fraction** — the *only*
variable the radius extraction and the convergence test read. The result is
fully reproducible from this ~29 GB subset; the other variables are not needed.
For the complete field data (all variables, for a different analysis), contact
the authors.

The `results/` CSVs are the small, human-readable evidence. The `D_*` field
directories are for anyone reproducing the extraction itself.

## Run configuration (identical physics, three grids)

| | 39³ | 99³ | 149³ |
|---|---|---|---|
| cells per octant edge | 39 | 99 | 149 |
| total cells | 59,319 | 970,299 | 3,307,949 |
| driving pressure (liquid) | 100 bar | 100 bar | 100 bar |
| bubble pressure | 1 bar (ratio 100) | 1 bar | 1 bar |
| initial radius R₀ | 1 mm | 1 mm | 1 mm |
| timestep dt | 6.26e-9 s | 2.466e-9 s | ~1.6e-9 s |
| target window | 0 → 0.5 t_c | 0 → 0.5 t_c | 0 → 0.5 t_c |
| **actually reached** | 0.50 t_c (completed) | **0.56 t_c (abort)** | **0.35 t_c (abort)** |
| Rayleigh collapse time t_c | 9.148 µs | 9.148 µs | 9.148 µs |

All three were launched to 0.5 t_c. The 39³ completed; the 99³ and 149³ aborted
as they approached the collapse focus, and — the key point — the *finer* 149³
aborted *earlier* (0.35 t_c) than the 99³ (0.56 t_c). Resolving the collapse
more sharply makes the instability arrive sooner.

## Field file format (read this before parsing D_*)

MFC writes ASCII `.dat` files, one per variable per saved timestep:

```
D_NN/prim.{VAR}.00.{STEP}.dat
```

Each line is **four whitespace-separated columns: x, y, z, value**. Coordinates
are in **code units** (physical length / 1 mm), so the stretched-grid cell
positions come with every file — no separate coordinate file.

**Variable mapping** (this build; MFC's ordering has moved between releases, so
do not assume it — it is confirmed here by the physical values):

| VAR | quantity | value in bubble / liquid |
|---|---|---|
| `prim.2` | **gas volume fraction** | 0.9 / 0.1 |
| `prim.7` | (fluid-2 related) | 0.1 / 0.9 |
| `prim.8` | (fluid-1 volume fraction) | 0.9 / 0.1 |

**Use `prim.2` for the gas fraction.** `prim.8` also sits in [0,1] but is the
water fraction — integrating it gives the whole domain, not the bubble. (We lost
time to exactly this; the value table above is how to tell them apart.)

## Reproducing the result

```bash
# extract R(t) from each grid
python code/post.py --run D_39  --alpha-var 2 --export-csv results/R_of_t_39.csv
python code/post.py --run D_99  --alpha-var 2 --export-csv results/R_of_t_99.csv
python code/post.py --run D_149 --alpha-var 2 --export-csv results/R_of_t_149.csv

# the three-grid convergence test — this is the finding
python code/post.py --compare3 D_39 D_99 D_149 --alpha-var 2
```

Expected: the three curves do **not** collapse onto one another. Over the common
window the 99³-vs-149³ difference reaches ~10% of R₀ through the steep
mid-collapse — refining from 99³ to 149³ does not bring it inside the ~5% you'd
want. The one quantity that *does* converge is the static initial (t=0) radius:
7.9% error at 39³, 2.6% at 99³, 2.2% at 149³, flattening — which localises the
failure to the collapse dynamics, not the discretisation.

## Radius extraction method

Radius is defined as the **α = 0.5 crossing of the radially-binned gas-fraction
profile** — resolution-independent, because the midpoint of a smoothed interface
does not move under refinement. An equivalent-volume integral is also computed
as a cross-check but is *not* the primary metric: it is grid-dependent because
the diffuse interface's physical thickness changes with resolution. If you
re-extract with a different radius definition and get a different convergence
behaviour, that is itself worth reporting.

## Caveats

- This is a **negative result**: the data demonstrates non-convergence, not a
  validated collapse. Do not use these curves as a validated bubble-collapse
  reference.
- The **collapse focus** (R → 0) is not in this data — all runs stop before it,
  and the finest grid stops earliest (0.35 t_c).
- Reaching a converged result requires adaptive mesh refinement at the collapse
  region; that is the conclusion, and it is not addressed here.

## Citation

If you use this data, please cite the archive DOI (Zenodo) and link the
accompanying code repository.

Apache License 2.0 (code) / CC-BY (data).