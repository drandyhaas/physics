r"""Build the arXiv upload: arxiv/ and arxiv_submission.tar.gz beside main.tex.

    python3 gen/make_arxiv.py        (after a build that kept main.bbl: tectonic -X compile --keep-intermediates main.tex)

arXiv compiles the source with pdfLaTeX and does NOT run BibTeX, so main.bbl goes in; it publishes the source, so the
author-notes comment block at the top of main.tex is stripped. Only the files the paper uses are included, and Table
I's rows (figs/table1.tex) are written into main.tex where it reads them: it reads them by the \input primitive, which
arXiv's scan of the upload does not recognize, so a separate file is marked unused and offered for deletion.
"""
import os
import shutil
import tarfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "arxiv")
FILES = ["main.bbl", "figs/k51_four.png", "figs/zynq_board_four.png", "figs/lp_cartoon.tex"]
TABLE = "\\csname @@input\\endcsname figs/table1.tex"
# the local build's font switch (XeTeX only; arXiv compiles with pdfLaTeX, which never takes it) is left out
XETEX = ["% (a XeTeX build -- tectonic -- sets its text", "\\usepackage{iftex}", "\\ifXeTeX\\usepackage[OT1]{fontenc}\\fi"]

shutil.rmtree(OUT, ignore_errors=True)
os.makedirs(os.path.join(OUT, "figs"))
lines = open(os.path.join(HERE, "main.tex")).read().split("\n")
# the leading block of full-line comments (build notes and the author's open items) is not for publication
i = 0
while i < len(lines) and lines[i].startswith("%"):
    i += 1
rows = [ln for ln in open(os.path.join(HERE, "figs", "table1.tex")).read().split("\n") if ln and not ln.startswith("%")]
body = []
for ln in lines[i:]:
    if any(ln.startswith(x) for x in XETEX):
        continue
    body += rows if ln.startswith(TABLE) else [ln]
assert len(body) == len(lines) - i - 1 + len(rows) - len(XETEX), "make_arxiv: an expected line is not there exactly once"
open(os.path.join(OUT, "main.tex"), "w").write("\n".join(body))
for f in FILES:
    shutil.copy(os.path.join(HERE, f), os.path.join(OUT, f))
tgz = os.path.join(HERE, "arxiv_submission.tar.gz")
with tarfile.open(tgz, "w:gz") as t:
    for f in ["main.tex"] + FILES:
        t.add(os.path.join(OUT, f), arcname=f)
print(f"stripped {i} comment lines, Table I's {len(rows)} rows inlined; wrote {tgz} ({os.path.getsize(tgz) // 1024} KB):", ["main.tex"] + FILES)
