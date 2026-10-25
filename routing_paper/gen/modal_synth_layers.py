#!/usr/bin/env python3
"""modal_synth_layers.py -- the paper's generated buses: KiCadRoutingTools' awx/synth_layers.py, every case on two,
three, four and six routing layers, one Modal container a case.

  KRT=/path/to/KiCadRoutingTools python3 gen/modal_synth_layers.py OUTDIR [--layers 2,3,4,6] [--only TAG,..]

Each container runs `synth_layers.py --layers N --only TAG --jobs 1` in KRT's own Modal image (awx/modal_whole.py,
run_stage: the checkout as it stands, built for Linux), and sends back the case's directory and its row of
layers.tsv. OUTDIR/TAG_LN/ holds them (row.tsv, run.log); a case already there is not run again. OUTDIR/layers.tsv
is every row, in order. KRT must be a checkout of the commit the paper names, with no local edits.
"""
import argparse
import csv
import io
import os
import subprocess
import sys
import tarfile

COMMIT = '09a3bf3b'                    # main.tex's \krtcommit

KRT = os.environ.get('KRT')
if not KRT or not os.path.isfile(os.path.join(KRT, 'awx', 'modal_whole.py')):
    sys.exit('modal_synth_layers: set KRT to a KiCadRoutingTools checkout')
AWX = os.path.join(os.path.abspath(KRT), 'awx')


def cases():
    out = subprocess.run([sys.executable, 'synth_layers.py', '--list'], cwd=AWX, capture_output=True, text=True).stdout
    return [ln.split()[0] for ln in out.splitlines() if ln.strip() and ' K=' in ln]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('outdir')
    ap.add_argument('--layers', default='2,3,4,6')
    ap.add_argument('--only')
    a = ap.parse_args()
    head = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=AWX, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=AWX, capture_output=True,
                           text=True).stdout.strip()
    if not head.startswith(COMMIT) or dirty:
        print(f'modal_synth_layers: WARNING: {KRT} is at {head[:8]}{" with local edits" if dirty else ""}, '
              f'the paper\'s results are of {COMMIT}', file=sys.stderr)
    tags = a.only.split(',') if a.only else cases()
    jobs = [(t, int(n)) for n in a.layers.split(',') for t in tags]
    out = os.path.abspath(a.outdir)
    todo = [(t, n) for t, n in jobs if not os.path.isfile(os.path.join(out, f'{t}_L{n}', 'row.tsv'))]
    print(f'{len(jobs)} runs, {len(todo)} to do', flush=True)
    sys.path.insert(0, AWX)
    import modal
    import modal_whole
    empty = io.BytesIO()
    tarfile.open(fileobj=empty, mode='w:gz').close()
    args = [(f'python3 synth_layers.py --layers {n} --only {t} --jobs 1 --outdir /tmp/sl', {}, empty.getvalue(),
             ['/tmp/sl/layers.tsv', f'/tmp/sl/{t}_L{n}']) for t, n in todo]
    with modal.enable_output(), modal_whole.app.run():
        for i, res in enumerate(modal_whole.run_stage.starmap(args, order_outputs=True)):
            t, n = todo[i]
            d = os.path.join(out, f'{t}_L{n}')
            os.makedirs(d, exist_ok=True)
            with tarfile.open(fileobj=io.BytesIO(res['tgz']), mode='r:gz') as tf:
                tf.extractall(d, filter='data')
            open(os.path.join(d, 'run.log'), 'w').write(res['log'])
            tsv = os.path.join(d, 'tmp', 'sl', 'layers.tsv')
            rows = list(csv.reader(open(tsv), delimiter='\t')) if os.path.isfile(tsv) else []
            with open(os.path.join(d, 'row.tsv'), 'w') as f:
                f.write('\n'.join('\t'.join(r) for r in rows) + '\n' if rows else f'{t}\tL{n}\tNO ROW rc={res["rc"]}\n')
            print(f'{t} L{n}: rc {res["rc"]} -- {rows[-1] if rows else "no row"}', flush=True)
    header, body = None, []
    for t, n in jobs:
        p = os.path.join(out, f'{t}_L{n}', 'row.tsv')
        if os.path.isfile(p):
            rows = [ln.rstrip('\n') for ln in open(p) if ln.strip()]
            if len(rows) > 1:
                header = header or rows[0]
            body.append(rows[-1])
    with open(os.path.join(out, 'layers.tsv'), 'w') as f:
        f.write('\n'.join(([header] if header else []) + body) + '\n')
    print(f'wrote {out}/layers.tsv', flush=True)


if __name__ == '__main__':
    main()
