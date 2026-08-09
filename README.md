---
title: "When the Bubble Won't Sit Still: A Negative Result on Cavitation Collapse"
description: We tried to validate a collapsing cavitation bubble against theory. Three grids gave three different answers — and the finest, at 3.3M cells, was the one that failed first. Here is why that is the finding, not the failure.
---

Case Study — Methods

# When the Bubble Won't Sit Still

## What we set out to do

Cavitation erosion — the pitting that eats pump impellers, choke trim, and valve
seats — comes from collapsing vapour bubbles. A low-pressure cavity forms, gets
swept into a high-pressure region, and implodes. At the final instant of that
implosion the collapse concentrates ambient pressure into a spike of thousands
of bar and, near a wall, fires a microjet at the surface. That spike and that
jet are what remove metal.

We wanted to simulate it. The plan was ordinary and, we thought, safe: take a
well-established compressible multiphase solver, collapse a single bubble, and
validate the collapse against the exact analytical solution — the
Rayleigh–Plesset equation — before trusting anything harder. Reproduce the known
answer, then push into the regime the analytics can't reach.

We never got to trust it. The validation itself failed, and it failed in an
instructive way. This post is about that failure, because a negative result
that's properly characterised is more useful than a positive one that isn't.

## The anchor that worked

The setup: a 1 mm bubble at low pressure, suddenly exposed to 100 bar of liquid,
collapsing under its own pressure deficit. In three dimensions, exploiting the
symmetry of a spherical collapse to mesh one octant. Water modelled as a
stiffened-gas liquid, the bubble as gas.

The analytical anchor is solid. The Rayleigh collapse time — the exact time for
an empty cavity to implode — is a closed-form expression with no free
parameters. Our thermodynamics reproduce it to five significant figures, and the
full Rayleigh–Plesset history (which adds the gas cushioning that arrests the
collapse) integrates cleanly. Before running any CFD, we had a trustworthy curve
to check against.

Then we ran the simulation, and the trouble started.

## Problem one: the focus is unreachable

A collapsing bubble is a genuinely violent thing. As the radius shrinks toward
zero, the gas inside compresses, the pressure diverges, and the local speed of
sound climbs toward five figures. On a fixed computational grid this is a
nightmare: the timestep is limited by the fastest wave in the domain, and near
the collapse focus that wave becomes almost arbitrarily fast.

Our runs aborted at roughly 56% of the collapse time, with the stability
indicator jumping to values that signalled the solution had gone unphysical, not
merely under-resolved. Crucially, shrinking the timestep did not help — because
the instability isn't a timestep problem. It's that a diffuse, smoothed
gas–liquid interface, compressed twentyfold on a fixed grid, produces a state
the model cannot represent. The final stage of the collapse — the exact stage
that matters for erosion — is out of reach.

[FIGURE 3: focus wall — Rayleigh–Plesset, both abort points, unreachable focus]

That alone reframed the project. If the peak pressure and the jet live at a
focus we can't reach, we can't deliver them. So we retreated to a more modest
goal: forget the focus, just validate the *early and middle* collapse — the part
that's stable — against Rayleigh–Plesset. Surely that much would hold.

## Problem two: it doesn't converge

Here is the test that decides everything, and it's the one people most often
skip. A simulation result is only trustworthy if it stops changing when you
refine the grid. So we ran the same collapse at three resolutions — 39³, 99³,
and 149³, from sixty thousand cells to 3.3 million — and asked a simple
question: do they agree on the bubble radius over time?

They don't.

[FIGURE 1: three grids diverging + pairwise gap panel]

The three curves take three visibly different paths through the collapse. The
coarse 39³ sits high throughout. The 99³ and 149³ start close — within about 1%
at the initial instant — but diverge through the steep part of the collapse,
where the finest grid drops as much as 10% below the middle one before they
cross again. Refining from 99³ to 149³ does not shrink the disagreement to
within the ~5% you would want; through the mid-collapse it makes it *worse*.

And the finer the grid, the sooner the run dies. The 99³ aborts near 0.56 t_c;
the 149³ — with three times the cells — aborts at 0.35 t_c, *earlier*. This is
the opposite of what a converging method does. Resolving the collapse more
sharply makes the instability arrive sooner.

To be precise about the claim, because it's the kind that invites a "you just
had it set up wrong" reply: this is not an artefact of one CFL number or one
reconstruction switch, and it isn't simply "still under-resolved" in the ordinary
sense where more cells would fix it — the finest grid is the one that fails
first. It is specific to *this* configuration: a fixed, uniform-in-the-collapse
octant grid with a diffuse interface. The same diffuse-interface family, given
interface sharpening or adaptive refinement that follows the collapse, does reach
the focus. What we're showing is that the fixed-grid configuration — the one a
reasonable person would reach for first — cannot, and adding cells uniformly
doesn't rescue it.

### The one thing that does converge

There is an important exception, and it's what tells us this is a resolution
problem rather than a broken solver. The *initial*, stationary bubble — a 1 mm
sphere sitting still, the easiest possible thing to represent — converges
cleanly:

[FIGURE 2: R0 error vs grid, converging and flattening]

Its radius error falls from 7.9% at 39³ to 2.6% at 99³ to 2.2% at 149³, and the
trend is flattening toward the ~1% of a well-resolved sphere. So the
discretisation itself is sound: give it something static and it converges. It's
the *collapse* — the moving, steepening, self-concentrating part — that refuses
to. The method works until the physics starts to concentrate, and then a fixed
grid can't keep up with it.

## Why this happens

Both problems have the same root. The solver represents the bubble surface as a
*diffuse* interface — smeared over a few cells rather than tracked sharply. That
choice buys robustness for many problems, but it has a fatal interaction with a
collapse: as the bubble shrinks, the interface region occupies a larger and
larger *fraction* of the bubble, and the smearing that was negligible at full
size dominates near the focus. Refining the grid narrows the smear in absolute
terms but the collapse also concentrates the physics into an ever-smaller region,
and on a fixed grid the resolution simply can't keep pace with the collapse.

This is a known hard problem in compressible multiphase CFD. Resolving a bubble
collapse to its focus needs either interface sharpening, adaptive mesh refinement
that follows the collapse and concentrates cells where the action is, or a
dedicated cavitation model. A fixed stretched grid at feasible resolution — even
a million cells — isn't enough.

## What this is worth

It would have been easy to stop at one grid. The 99³ result, taken alone, looks
plausible: a bubble collapsing on roughly the right timescale, tracking the
analytical curve to within a few percent over the stable window. Published with a
single figure and no convergence study, it would have passed a casual read. It
would also have been wrong — or more precisely, unsupported, which in engineering
is the same thing. It took a coarser grid to show the 99³ wasn't obviously right,
and a finer one to show it wasn't right at all.

Concretely: at the same instant in the collapse, the 99³ grid reads a bubble
radius that the 149³ puts up to 10% lower — a ~25% disagreement in bubble
*volume*, and a correspondingly wrong collapse rate. Anyone who took the 99³
curve alone as the answer would have misreported the collapse dynamics by that
margin and had no way to know it. The convergence test is what separates a result
from a picture. Three grids, one question — *does the answer stop moving?* — and
the honest answer here was no; the third grid is what made that unarguable.

So we're not publishing a cavitation-erosion prediction, because we don't have
one we can stand behind. What we have instead is a clear, quantified boundary:

- The **analytical anchor** (Rayleigh collapse time, Rayleigh–Plesset) is exact
  and reproduced to five figures. That part is trustworthy.
- The **collapse focus** — where erosion physics lives — is unreachable on a
  fixed grid with a diffuse interface, and no timestep reduction fixes it.
- The **collapse dynamics** short of the focus **do not converge** across 39³,
  99³, and 149³ — the finest grid disagrees with the middle one by up to 10%
  through the steep collapse, and aborts *earlier* than the coarser run rather
  than later. Even the stable window isn't validated at 3.3M cells.
- The **static initial sphere does converge**, which localises the problem to
  the collapse itself, not the discretisation.
- Reaching a trustworthy result requires **adaptive mesh refinement** at the
  collapse region — a real piece of development, not a parameter change.

That's a map of exactly where the hard part is and what it would take to get
through it. For anyone scoping cavitation-erosion simulation, that map is worth
more than an unconverged number dressed up as an answer.

## The general point

There's a failure mode in computational engineering where a simulation produces a
plausible picture, the picture matches intuition, and everyone moves on. The
picture might be right. But without a convergence study you don't know, and
"looks reasonable" is not a validation. The discipline that catches this is
boring and cheap next to the run itself: do it at more than one resolution — and
if you can, more than two, because here it was the *third* grid that turned "the
coarse one is probably just too coarse" into "refinement is actively making this
worse." When the answer won't stop moving, that *is* the result — and reporting
it honestly is how the field learns where its tools actually stop.

---

*Reproducibility: the analytical anchor, case setup, and the radius-extraction
and grid-comparison post-processing are open source. The Rayleigh–Plesset
comparison runs in seconds; the 3D collapse runs on GPU. The non-convergence
above is reproducible from the three provided grids (39³, 99³, 149³); the
central claim — that refinement makes the disagreement worse — is the 99³-vs-149³
comparison specifically.*