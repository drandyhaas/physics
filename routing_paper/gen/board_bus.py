#!/usr/bin/env python3
"""board_bus.py -- the Zynq DDR3 bus on the board itself: the whole route's bus step, the toolkit's production chain and
Freerouting on the same board, nets, layers and sizes, and the designer's board, all graded alike.

  KRT=/path/to/KiCadRoutingTools python3 gen/board_bus.py MODE OUTROOT [--bus BUS.kicad_pcb]
  KRT=... python3 gen/board_bus.py render OUTROOT OUT.png        the four graded boards drawn alike
  KRT=... python3 gen/board_bus.py table OUTROOT                 data/board/*.json and figs/table2.tex (Table II)

MODE:
  nets   OUTROOT/nets.txt: the bus step's run, the nets route_bus gave the whole route (BUS.bus/run/r1/nets.lines)
  krt    the production chain on PLANES (KRT/awx/tmp/zynq/planes.kicad_pcb, the awx README's recipe): bga_fanout at
         U1 and at U2 for the run's nets (no plane drops, as the bus step lays none), route_diff on its pairs, route.py
         on its singles -- every step on F.Cu and B.Cu, at the bus step's sizes, every net class set to them
  fr     Freerouting 2.4.1 with its defaults on PLANES, the run's nets only (KRT/awx/baseline_freerouting.py), its
         inner layers marked as power planes in its own copy so that it routes on F.Cu and B.Cu as the others do
  grade  every board -- OUTROOT/ours/bus.kicad_pcb (--bus), krt/route.kicad_pcb, fr/out.kicad_pcb and the designer's
         board -- graded on the run's nets under ONE project, PLANES's with every net class at the bus step's sizes:
         baseline_bench.grade, OUTROOT/MODE/result.json

The bus step itself is run by the awx README's own command ("The bus step on a real board"), writing
OUTROOT/ours/bus.kicad_pcb and its log OUTROOT/ours/bus.log.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import time

KRT = os.environ.get('KRT')
if not KRT or not os.path.isfile(os.path.join(KRT, 'awx', 'baseline_bench.py')):
    sys.exit('board_bus: set KRT to a KiCadRoutingTools checkout')
AWX = os.path.join(os.path.abspath(KRT), 'awx')
PYR = os.path.join(os.path.abspath(KRT), 'py_router')
sys.path.insert(0, AWX)
sys.path.insert(0, PYR)
PLANES = os.path.join(AWX, 'tmp', 'zynq', 'planes.kicad_pcb')
HUMAN = os.path.join(AWX, 'tmp', 'zynq', 'src', 'boards_set1', 'zynq_ad9364.kicad_pcb')
SZ = ['--track-width', '0.15', '--via-size', '0.45', '--via-drill', '0.3']      # the bus step's sizes
CLR, FAN = '0.09', '0.12'
LAYERS = ['F.Cu', 'B.Cu']
EXCLUDE = {'DDR3_RESET'}       # (grade: the bus step's lane, its pull's stub the chain's -- see grade)
SIZE_KINDS = {'TRACK-WIDTH', 'VIA-SIZE', 'VIA-DRILL-SIZE', 'ANNULAR-WIDTH', 'HOLE-SIZE', 'DRILL-SIZE'}


def nets_of(root):
    return [n for n in open(os.path.join(root, 'nets.txt')).read().split() if n]


def copy_board(src, dst):
    subprocess.run([sys.executable, os.path.join(PYR, 'copy_board.py'), src, dst], check=True,
                   stdout=subprocess.DEVNULL)


def normalize_project(pro):
    """every net class at the bus step's sizes, as baseline_bench.prep does on the benches: the board's own classes
    (wider tracks, a larger clearance) do not escape its arrays, and every router is to route and be graded alike"""
    pj = json.load(open(pro))
    for c in pj['net_settings']['classes']:
        c.update(clearance=float(CLR), track_width=0.15, via_diameter=0.45, via_drill=0.3, diff_pair_width=0.15)
    json.dump(pj, open(pro, 'w'), indent=2)


def run(log, argv):
    t = time.time()
    with open(log, 'w') as f:
        rc = subprocess.run(argv, cwd=AWX, stdout=f, stderr=subprocess.STDOUT).returncode
    return rc, round(time.time() - t, 1)


def main():
    mode, root = sys.argv[1], os.path.abspath(sys.argv[2])
    bus = sys.argv[sys.argv.index('--bus') + 1] if '--bus' in sys.argv else os.path.join(root, 'ours', 'bus.kicad_pcb')
    if mode == 'nets':
        src = os.path.splitext(bus)[0] + '.bus/run/r1/nets.lines'     # the whole route's run, as route_bus gave it
        os.makedirs(root, exist_ok=True)
        shutil.copy(src, os.path.join(root, 'nets.txt'))
        print(f'{len(nets_of(root))} nets')
        return
    import pairs as _pairs
    nets = nets_of(root)
    legs = [l for pr in _pairs.pair_names(nets, admit_all=True).values() for l in pr]
    singles = [n for n in nets if n not in legs]
    pats = lambda ns: [f'*{n}' for n in ns]
    out = os.path.join(root, mode)
    os.makedirs(out, exist_ok=True)
    I = lambda f: os.path.join(out, f)
    timing = {}
    if mode == 'krt':
        copy_board(PLANES, I('in.kicad_pcb'))
        normalize_project(I('in.kicad_pro'))
        cur = I('in.kicad_pcb')
        for part in ('U1', 'U2'):
            fo = [sys.executable, os.path.join(PYR, 'bga_fanout.py'), cur, '--output', I(f'fo_{part}.kicad_pcb'),
                  '--component', part, '--layers'] + LAYERS + ['--nets'] + pats(nets)
            if legs:
                fo += ['--diff-pairs'] + pats(legs)
            rc, timing[f'fanout_{part}'] = run(I(f'fanout_{part}.log'), fo + ['--track-width', FAN, '--clearance', CLR,
                                                                            '--plane-drop', 'off'] + SZ[2:])
            cur = I(f'fo_{part}.kicad_pcb')
        sizes = ['--clearance-ceiling', CLR, '--layers'] + LAYERS + SZ
        if legs:
            rc, timing['route_diff'] = run(I('route_diff.log'), [sys.executable, os.path.join(PYR, 'route_diff.py'), cur,
                                           '--output', I('diff.kicad_pcb'), '--nets'] + pats(legs) +
                                           ['--no-gnd-vias'] + sizes)
            cur = I('diff.kicad_pcb')
        rc, timing['route'] = run(I('route.log'), [sys.executable, os.path.join(PYR, 'route.py'), cur, '--output',
                                  I('route.kicad_pcb'), '--nets'] + pats(singles) + sizes)
    elif mode == 'fr':
        import baseline_freerouting
        copy_board(PLANES, I('in.kicad_pcb'))
        normalize_project(I('in.kicad_pro'))
        txt = open(I('in.kicad_pcb')).read()
        txt2 = re.sub(r'(\(\d+ "In\d+\.Cu") signal', r'\1 power', txt)
        assert txt2 != txt, 'board_bus: no inner signal layer to mark as a plane'
        open(I('in.kicad_pcb'), 'w').write(txt2)
        timing.update(baseline_freerouting.route(I('in.kicad_pcb'), nets, 'U1', out, fanout=True))
    elif mode == 'grade':
        import baseline_bench as bb
        from kicad_parser import parse_kicad_pcb
        import contextlib
        import io
        # the shared project: the bus step's output project (the fab floor it wrote), every class normalized. Each
        # router writes the floor it routed to, so the SIZE findings are reported apart -- each board's smallest track,
        # via and drill -- and DRC counts the electrical ones alone, at the one clearance
        pro = os.path.join(root, 'graded.kicad_pro')
        shutil.copy(os.path.splitext(bus)[0] + '.kicad_pro', pro)
        normalize_project(pro)
        # RESET is the bus step's lane between the arrays, its pull's stub left to the rest of the chain (route_bus:
        # a PULL beside the ends): graded on the other nets, as the others route the whole net
        nets = [n for n in nets if n not in EXCLUDE]
        legs = [l for l in legs if l not in EXCLUDE]
        for m, board in (('ours', bus), ('krt', os.path.join(root, 'krt', 'route.kicad_pcb')),
                         ('fr', os.path.join(root, 'fr', 'out.kicad_pcb')), ('human', HUMAN)):
            o = os.path.join(root, m)
            os.makedirs(o, exist_ok=True)
            g = os.path.join(o, 'graded.kicad_pcb')
            shutil.copy(board, g)
            shutil.copy(pro, os.path.join(o, 'in.kicad_pro'))
            r = dict(mode=m, nets=len(nets), pairs=len(legs) // 2, excluded=sorted(EXCLUDE))
            r.update(bb.grade(o, g, None, nets, own_pro=False, others=False))
            r['drc_size'] = {k: v for k, v in r['drc_kinds'].items() if k in SIZE_KINDS}
            r['drc'] = sum(v for k, v in r['drc_kinds'].items() if k not in SIZE_KINDS)
            with contextlib.redirect_stdout(io.StringIO()):
                p = parse_kicad_pcb(g)
            nm = {i: n.name.split('/')[-1] for i, n in p.nets.items()}
            segs = [x for x in p.segments if nm.get(x.net_id) in set(nets)]
            vias = [x for x in p.vias if nm.get(x.net_id) in set(nets)]
            r['min_track'] = round(min((x.width for x in segs), default=0), 4)
            r['min_via'] = round(min((x.size for x in vias), default=0), 4)
            r['min_drill'] = round(min((x.drill for x in vias), default=0), 4)
            prev = os.path.join(o, 'result.json')
            if os.path.isfile(prev):
                r['timing'] = json.load(open(prev)).get('timing', {})
            json.dump(r, open(prev, 'w'), indent=1)
            print(f"{m:6s} open {r['open']} {r['open_nets']} drc {r['drc']} vias {r['vias']} mm {r['mm']} "
                  f"min track/via/drill {r['min_track']}/{r['min_via']}/{r['min_drill']} (size findings "
                  f"{r['drc_size']})", flush=True)
        return
    elif mode == 'table':
        # the grades kept beside the paper (data/board/MODE.json) and Table II's rows from them (figs/table2.tex)
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        data = os.path.join(here, 'data', 'board')
        os.makedirs(data, exist_ok=True)
        rows = []
        for m, name in (('ours', 'whole route'), ('krt', 'KRT chain'), ('fr', 'Freerouting'), ('human', 'human')):
            src = os.path.join(root, m, 'result.json')
            if not os.path.isfile(src):
                sys.exit(f'board_bus: no grade for {m}')
            shutil.copy(src, os.path.join(data, f'{m}.json'))
            r = json.load(open(src))
            rows.append(f"{name} & {r['open']} & {r['drc']} & {r['vias']} & {r['mm']} & "
                        f"{r['min_track']:g} & {r['min_via']:g} \\\\")
        out_ = os.path.join(here, 'figs', 'table2.tex')
        open(out_, 'w').write('% generated by gen/board_bus.py table from data/board -- do not edit\n' +
                              '\n'.join(rows) + '\n')
        print(f'wrote {out_}')
        return
    elif mode == 'render':
        # the four boards drawn alike (KRT's awx/baseline_render.py, as Fig. 5), labelled from their grades
        png = os.path.abspath(sys.argv[3])
        names = (('ours', 'whole route'), ('krt', 'KRT chain'), ('fr', 'Freerouting'), ('human', 'human'))
        items = []
        for m, name in names:
            r = json.load(open(os.path.join(root, m, 'result.json')))
            # (portrait panels in one row: the label its vias and copper, its open nets and DRC in Table II)
            lab = f"{name}: {r['vias']} vias, {r['mm']} mm"
            items.append(f"{lab}={os.path.join(root, m, 'graded.kicad_pcb')}")
        graded = [n for n in nets if n not in EXCLUDE]
        subprocess.run([sys.executable, os.path.join(AWX, 'baseline_render.py'), png, ','.join(graded)] + items +
                       ['--no-others', '--cols', '4'], cwd=AWX, check=True)
        return
    else:
        sys.exit(f'board_bus: unknown mode {mode}')
    json.dump(dict(mode=mode, nets=len(nets), timing=timing), open(I('result.json'), 'w'), indent=1)
    print(mode, timing)


if __name__ == '__main__':
    main()
