# Whole-Route Planning of a Memory Bus Between Two BGAs

A paper, in draft for *IEEE Transactions on Computer-Aided Design of Integrated
Circuits and Systems*, on routing the bus between a processor and its DRAM
fully automatically — the densest wiring on most processor boards, and still
largely routed by hand.

📄 **[main.pdf](main.pdf)** (a draft)

The router plans every lane of a two-layer BGA-to-BGA bus completely before
any copper is laid, and carries the plan to finished, design-rule-clean copper.
On two public boards, each routed by its designer, it completes every rung of
a ladder of up to 51 nets connected and DRC-clean, never with more vias than
the human and always with less copper:

| board | nets | vias: router / human | copper (mm): router / human |
| --- | --- | --- | --- |
| Allwinner H3 → DDR3 | 51 | **80** / 88 | **1220** / 1337 |
| Zynq-7020 → DDR3 | 44 | **58** / 103 | **1347** / 1617 |

## The idea

**The ends decide the braid.** Once each net's escape from both arrays is
fixed, the two orders of the escapes fix which pairs of nets must cross — the
inversions of a permutation. What is left free is where each crossing happens,
in what order, and where each lane changes layer. These choices interact, since
a via placed to serve one crossing serves or spoils the next, so the router
decides all of them at once:

1. **An exact constraint-programming model** (CP-SAT) chooses every crossing
   position, the crossing order and every layer change jointly, and, for the ends
   it is given, proves the fewest layer changes within the model. Its consistency
   rule over triples of lanes
   admits exactly the orders a wiring diagram of the permutation can realize.
2. **A linear program** places the lanes in millimeters. The constraints it
   has to break go back to the model as cuts, and an audit's findings go back as
   negotiated congestion prices, until the plan passes every rule of the
   production router.
3. **The escapes at both arrays** are chosen by an estimate of the same
   objective, so the ends themselves are not proved best. A
   differential pair is planned as one lane with a polarity at each end, and
   gets an automatic via crossover where the two ends disagree.

## The physics in it

The models read naturally as physics, and the readings were of practical use:

| routing | physics |
| --- | --- |
| a lane along the route coordinate | a worldline in one space and one time dimension |
| a crossing | a two-body scattering event |
| the triple rule on crossing order | the Yang–Baxter relation (third Reidemeister move) |
| a lane's layer | an Ising spin |
| crossing lanes must be on different layers | an antiferromagnetic bond |
| a via | a domain wall along the worldline |
| why vias are forced at all | frustration |
| the geometry LP | statics: strings at unit tension, struts between neighbors, walls at pads |
| its dual variables | contact forces, with complementary slackness as the contact law |
| an elastic constraint | a rigid-plastic strut that yields |

The forces are also useful for reading a plan: in the first pass on the 51-net
rung only 376 of the hard constraints carry a force, and they name the rules
that hold each lane where it is. A constraint the program could not keep is a
strut that yields, and those are what go back to the solve as cuts.

## Files

| file | what it is |
| --- | --- |
| `main.tex` | the paper (IEEEtran journal class); Figs. 1, 2 and 4 are TikZ inline |
| `refs.bib` | the bibliography |
| `main.pdf` | the built paper |
| `figs/lp_cartoon.tex` | Fig. 3, a toy of the geometry LP solved, with its real offsets and dual forces; generated |
| `figs/k51_pair.png` | Fig. 5, the 51-net H3 rung beside the human's board; generated |
| `gen/lp_cartoon.py` | solves the toy LP and writes `figs/lp_cartoon.tex` |
| `gen/fig_board_pair.py` | draws two boards' bus copper side by side, identically, and writes Fig. 5 |
| `gen/make_arxiv.py` | builds the arXiv upload, `arxiv/` and `arxiv_submission.tar.gz` (neither committed) |
| `figs/k51_own_ends.png`, `figs/k51_whole_route.png`, `figs/zynq_k44_whole.png`, `figs/k51_forces.png` | dark-theme renders of the H3 rung on our ends (the earlier Fig. 5) and on the human's ends, of the Zynq rung, and of the LP's binding rules; not currently used in the paper |

## Building

```bash
tectonic -X compile --keep-intermediates main.tex   # or pdflatex, bibtex, then pdflatex twice
python3 gen/lp_cartoon.py                           # regenerates Fig. 3 (numpy, scipy)
python3 gen/make_arxiv.py                           # the arXiv upload; needs the main.bbl the build kept
```

arXiv compiles with pdfLaTeX and does not run BibTeX, so the upload carries
`main.bbl`. It publishes the source, so the script strips the comment block of
build notes and open items from the top of `main.tex`.

`gen/fig_board_pair.py` needs the router's `py_router` package, found through
the `KRT` environment variable, and takes the two `.kicad_pcb` boards and a
list of the run's nets:

```bash
KRT=/path/to/KiCadRoutingTools python3 gen/fig_board_pair.py \
    OURS.kicad_pcb HUMAN.kicad_pcb NETS.lines OUT.png "left label" "right label"
```

## Where the numbers come from

This repository's convention is that a document's numbers come from a script
beside it. Most of this paper's do not: they are results of runs of the router,
which is open source with the benches and recipes that reproduce them:
[KiCadRoutingTools](https://github.com/drandyhaas/KiCadRoutingTools). Only
Fig. 3 is computed here.

## Status

A draft. The author's open items are listed in a comment at the top of
`main.tex`.
