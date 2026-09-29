"""Two boards' bus copper side by side, drawn identically, for the paper (Fig. 5).

    python3 gen/fig_board_pair.py OURS.kicad_pcb HUMAN.kicad_pcb NETS.lines OUT.png "left label" "right label"

Only the run's nets are drawn (both boards), over every pad; the other nets' copper and the planes are left out.
F is red and B blue (the renderer's light palette); a differential pair's legs are drawn in PAIR; every via of
the run's nets is the SAME marker on both boards (a ring of VIA_MM), not its drawn size, so the two boards read
alike. White background, black labels, cropped to the copper. Needs the repo's py_router (KRT below).
"""
import os
import sys

from PIL import Image, ImageDraw, ImageFont

KRT = os.environ.get("KRT", "/Users/andy/Documents/KiCadRoutingTools/.claude/worktrees/bus622-take5")
sys.path.insert(0, os.path.join(KRT, "py_router"))
import render_theme                      # noqa: E402
from kicad_parser import parse_kicad_pcb  # noqa: E402
from route_render import BoardRenderer    # noqa: E402

PAIR = (230, 140, 0)          # a differential pair's legs
VIA_MM = 0.5                  # every via's marker diameter
VIA_RING = (40, 40, 40)
SIZE, SS = 1800, 3
FONT = "/System/Library/Fonts/Helvetica.ttc"

L = render_theme.LIGHT
PAPER = render_theme.Theme("paper", {**L._rgb, "ground": (255, 255, 255), "board_body": (255, 255, 255),
                                     "board_edge": (200, 200, 200), "pad": (175, 160, 120)},
                           L._mark, L.layers, 255)


def bus(pcb, names):
    ids = {i for i, n in pcb.nets.items() if n.name.split("/")[-1] in names}
    pairs = {n[:-1] for n in names if n[-1] in "PN" and (n[:-1] + ("N" if n[-1] == "P" else "P")) in names}
    pair_ids = {i for i, n in pcb.nets.items() if n.name.split("/")[-1][:-1] in pairs
                and n.name.split("/")[-1] in names}
    segs = [s for s in pcb.segments if s.net_id in ids]
    vias = [v for v in pcb.vias if v.net_id in ids]
    return segs, vias, pair_ids


def bbox(segs, vias):
    xs = [s.start_x for s in segs] + [s.end_x for s in segs] + [v.x for v in vias]
    ys = [s.start_y for s in segs] + [s.end_y for s in segs] + [v.y for v in vias]
    return min(xs), min(ys), max(xs), max(ys)


def panel(path, names, view):
    pcb = parse_kicad_pcb(path)
    segs, vias, pair_ids = bus(pcb, names)
    r = BoardRenderer(pcb, size=SIZE, supersample=SS, show_zones=False, view=view, theme=PAPER, layer_alpha=255)

    def markers(d, rr):
        rad = rr.tf.length(VIA_MM / 2)
        w = max(1, int(round(rr.tf.length(0.07))))
        for v in vias:
            x, y = rr.tf.pt(v.x, v.y)
            d.ellipse([x - rad, y - rad, x + rad, y + rad], fill=(255, 255, 255), outline=VIA_RING, width=w)

    img = r.frame(segments=[s for s in segs if s.net_id not in pair_ids], vias=[],
                  highlight_segments=[s for s in segs if s.net_id in pair_ids], highlight_color=PAIR,
                  overlays=[markers])
    # the canvas is the whole board's; crop to the view
    (x0, y0), (x1, y1) = r.tf.pt(view[0], view[1]), r.tf.pt(view[2], view[3])
    k = img.width / (r.tf.pt(r.bounds[2], r.bounds[3])[0] + r.tf.pt(r.bounds[0], r.bounds[1])[0]) if False else 1 / SS
    box = [int(min(x0, x1) * k), int(min(y0, y1) * k), int(max(x0, x1) * k) + 1, int(max(y0, y1) * k) + 1]
    return img.crop(box), len(vias)


def main():
    ours, human, nets_file, out, lab_l, lab_r = sys.argv[1:7]
    names = {n for ln in open(nets_file) for n in ln.split()}
    boxes = [bbox(*bus(parse_kicad_pcb(p), names)[:2]) for p in (ours, human)]
    m = 1.0
    view = (min(b[0] for b in boxes) - m, min(b[1] for b in boxes) - m,
            max(b[2] for b in boxes) + m, max(b[3] for b in boxes) + m)
    (a, na), (b, nb) = panel(ours, names, view), panel(human, names, view)
    print(f"via markers: ours {na}, human {nb}; view {tuple(round(v, 2) for v in view)}")
    font = ImageFont.truetype(FONT, int(a.height * 0.055))
    top, gap = int(a.height * 0.09), int(a.width * 0.03)
    W, H = a.width + gap + b.width, top + max(a.height, b.height)
    img = Image.new("RGB", (W, H), (255, 255, 255))
    img.paste(a, (0, top))
    img.paste(b, (a.width + gap, top))
    d = ImageDraw.Draw(img)
    d.text((int(a.width * 0.01), int(top * 0.15)), lab_l, fill=(0, 0, 0), font=font)
    d.text((a.width + gap + int(b.width * 0.01), int(top * 0.15)), lab_r, fill=(0, 0, 0), font=font)
    img.save(out)
    print("wrote", out, img.size)


if __name__ == "__main__":
    main()
