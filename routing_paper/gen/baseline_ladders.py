#!/usr/bin/env python3
"""baseline_ladders.py -- Table I and Fig. 5: KiCadRoutingTools' awx/baseline_bench.py on the paper's two ladders.

  KRT=/path/to/KiCadRoutingTools python3 gen/baseline_ladders.py ARGS...   (ARGS as baseline_bench.py takes them)

baseline_bench.py runs a rung by the toolkit's production router (krt) and by Freerouting (fr), grades those boards,
the whole route's (ours) and the designer's (human) the same way, and writes the table and the drawings. Its own
LADDERS name the article's FIRST Zynq build, which no longer exists; this runs it on the paper's ladders instead:

  h3    awx/fb_t2q_pairs.kicad_pcb (U1 -> DU1), rungs K15 K28 K35 K41 K51, as baseline_bench.py has them
  zynq  awx/tmp/zynq/zynqF.kicad_pcb (U1 -> U2), built from the committed awx/boards/zynq_ad9364.kicad_pcb.gz by the
        awx README's "Building the zynq article"; its ladder's checkpoints are 10 18 26 33 38 43 48, and the paper
        takes 18 to 43

KRT must be a checkout of the commit the paper names (\\krtcommit in main.tex); a different one is run, with a warning.
"""
import os
import subprocess
import sys

COMMIT = '2e319092'                    # main.tex's \krtcommit

KRT = os.environ.get('KRT')
if not KRT or not os.path.isfile(os.path.join(KRT, 'awx', 'baseline_bench.py')):
    sys.exit('baseline_ladders: set KRT to a KiCadRoutingTools checkout')
AWX = os.path.join(os.path.abspath(KRT), 'awx')
head = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=AWX, capture_output=True, text=True).stdout.strip()
if not head.startswith(COMMIT):
    print(f'baseline_ladders: WARNING: {KRT} is at {head[:8]}, the paper\'s results are of {COMMIT}', file=sys.stderr)
if not os.path.isfile(os.path.join(AWX, 'tmp', 'zynq', 'zynqF.kicad_pcb')):
    print('baseline_ladders: WARNING: no awx/tmp/zynq/zynqF.kicad_pcb -- build it first (awx README, '
          '"Building the zynq article")', file=sys.stderr)

# Freerouting's score charges a clearance violation 1e6 against 5e6 an unrouted connection and stops its autorouter
# after ten passes without a gain. Under its octagonal shapes the bench's FIXED copper, legal at the board's
# clearance, reads as about a hundred violations it can never remove, which pin a small rung's score at zero, so it
# stopped early (the Zynq's K18, K26 and K33 at 6, 10 and 13 open); the paper runs it with that charge at 1e3
os.environ.setdefault('FREEROUTING__ROUTER__SCORING__CLEARANCE_VIOLATION_PENALTY', '1000')

sys.path.insert(0, AWX)
import baseline_bench as bb           # noqa: E402

bb.LADDERS = {
    'h3': dict(bench='fb_t2q_pairs.kicad_pcb', dest='DU1', src='U1', rungs=(15, 28, 35, 41, 51),
               human='fb_t2q_human.kicad_pcb', turn=None),
    'zynq': dict(bench='tmp/zynq/zynqF.kicad_pcb', dest='U2', src='U1', rungs=(18, 26, 33, 38, 43),
                 human='tmp/zynq/src/boards_set1/zynq_ad9364.kicad_pcb', turn=(1, 109.3, -103.2)),
}

if __name__ == '__main__':
    bb.main(sys.argv[1:])
