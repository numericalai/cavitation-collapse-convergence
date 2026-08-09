# Cavitation bubble collapse — a grid-convergence study (negative result)

The collapse does not grid-converge on a uniform octant mesh, even at 3.3M
cells, and the finest grid fails earliest. See post/blog_negative_result.md for
the full write-up, figures/ for the evidence, and results/ for the extracted
curves. Raw field data is archived separately (results/DATA_ARCHIVE.md).

Reproduce:
    python code/thermo.py                                  # verify the anchor
    python code/post.py --compare3 D_39 D_99 D_149 --alpha-var 2

Built on MFC (https://github.com/MFlowCode/MFC), which is MIT-licensed. The code
in this repository is Apache-2.0; data/figures are CC-BY.
