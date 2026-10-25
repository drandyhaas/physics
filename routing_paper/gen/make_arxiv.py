"""Build the arXiv upload: arxiv/ and arxiv_submission.tar.gz beside main.tex.

    python3 gen/make_arxiv.py        (after a build that kept main.bbl: tectonic -X compile --keep-intermediates main.tex)

arXiv compiles the source with pdfLaTeX and does NOT run BibTeX, so main.bbl goes in; it publishes the source, so the
author-notes comment block at the top of main.tex is stripped. Only the files the paper uses are included.
"""
import os
import shutil
import tarfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "arxiv")
FILES = ["main.bbl", "figs/k51_four.png", "figs/zynq_k43_four.png", "figs/lp_cartoon.tex", "figs/table1.tex"]

shutil.rmtree(OUT, ignore_errors=True)
os.makedirs(os.path.join(OUT, "figs"))
lines = open(os.path.join(HERE, "main.tex")).read().split("\n")
# the leading block of full-line comments (build notes and the author's open items) is not for publication
i = 0
while i < len(lines) and lines[i].startswith("%"):
    i += 1
open(os.path.join(OUT, "main.tex"), "w").write("\n".join(lines[i:]))
for f in FILES:
    shutil.copy(os.path.join(HERE, f), os.path.join(OUT, f))
tgz = os.path.join(HERE, "arxiv_submission.tar.gz")
with tarfile.open(tgz, "w:gz") as t:
    for f in ["main.tex"] + FILES:
        t.add(os.path.join(OUT, f), arcname=f)
print(f"stripped {i} comment lines; wrote {tgz} ({os.path.getsize(tgz) // 1024} KB):", ["main.tex"] + FILES)
