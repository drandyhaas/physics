"""A tiny dependency-free vector-PDF canvas built on macOS Quartz + AppKit.

Exists because this repo's figures want real vector output (no rasterised
matplotlib, no LaTeX install) with exact text metrics and Unicode maths.
Text goes through NSAttributedString, so glyph fallback is automatic.

Coordinates are PDF points, origin bottom-left, y up. `text()` takes the
BASELINE y, not the bounding-box bottom.
"""
import math
from Quartz import (
    CGPDFContextCreateWithURL, CGRectMake, CGContextBeginPage, CGContextEndPage,
    CGPDFContextClose, CGContextSaveGState, CGContextRestoreGState,
    CGContextSetRGBFillColor, CGContextSetRGBStrokeColor, CGContextSetLineWidth,
    CGContextSetLineCap, CGContextSetLineJoin, CGContextSetLineDash,
    CGContextBeginPath, CGContextMoveToPoint, CGContextAddLineToPoint,
    CGContextAddCurveToPoint, CGContextAddQuadCurveToPoint, CGContextAddArc,
    CGContextClosePath, CGContextFillPath, CGContextStrokePath, CGContextDrawPath,
    CGContextClip, CGContextAddEllipseInRect, CGContextFillEllipseInRect,
    CGContextStrokeEllipseInRect, CGContextFillRect, CGContextAddRect,
    CGContextSetAlpha, CGContextTranslateCTM, CGContextRotateCTM,
    CGContextScaleCTM, CGContextSetFlatness,
    kCGLineCapRound, kCGLineCapButt, kCGLineJoinRound, kCGPathFillStroke,
    kCGPathStroke, kCGPathFill,
)
from CoreFoundation import CFURLCreateWithFileSystemPath, kCFURLPOSIXPathStyle
from AppKit import (NSGraphicsContext, NSAttributedString, NSFont, NSColor,
                    NSFontAttributeName, NSForegroundColorAttributeName,
                    NSKernAttributeName)

# ---------------------------------------------------------------- fonts
_FALLBACKS = {
    'body':  ["Charter-Roman", "Palatino-Roman", "TimesNewRomanPSMT", "Times-Roman"],
    'bodyi': ["Charter-Italic", "Palatino-Italic", "TimesNewRomanPS-ItalicMT", "Times-Italic"],
    'bodyb': ["Charter-Bold", "Palatino-Bold", "TimesNewRomanPS-BoldMT", "Times-Bold"],
    'sans':  ["HelveticaNeue", "Helvetica"],
    'sansm': ["HelveticaNeue-Medium", "Helvetica"],
    'sansb': ["HelveticaNeue-Bold", "Helvetica-Bold"],
    'math':  ["STIXGeneral-Regular", "TimesNewRomanPSMT", "Times-Roman"],
    'mathi': ["STIXGeneral-Italic", "TimesNewRomanPS-ItalicMT", "Times-Italic"],
    'mono':  ["Menlo-Regular", "Courier"],
}
_cache = {}

def font(role, size):
    key = (role, round(size, 2))
    if key in _cache:
        return _cache[key]
    for name in _FALLBACKS[role]:
        f = NSFont.fontWithName_size_(name, size)
        if f is not None:
            _cache[key] = f
            return f
    f = NSFont.systemFontOfSize_(size)
    _cache[key] = f
    return f

def font_report():
    return {r: font(r, 10).fontName() for r in _FALLBACKS}

# ---------------------------------------------------------------- canvas
class Canvas:
    def __init__(self, path, width, height, info=None):
        url = CFURLCreateWithFileSystemPath(None, path, kCFURLPOSIXPathStyle, False)
        self.box = CGRectMake(0, 0, width, height)
        self.w, self.h = width, height
        aux = {}
        for k, v in (info or {}).items():
            aux[k] = v
        self.ctx = CGPDFContextCreateWithURL(url, self.box, aux or None)
        self._ns = None

    # -- pages ---------------------------------------------------------
    def page(self):
        CGContextBeginPage(self.ctx, self.box)
        CGContextSetFlatness(self.ctx, 0.05)
        self._ns = NSGraphicsContext.graphicsContextWithCGContext_flipped_(self.ctx, False)

    def endpage(self):
        CGContextEndPage(self.ctx)

    def close(self):
        CGPDFContextClose(self.ctx)

    # -- state ---------------------------------------------------------
    def save(self):    CGContextSaveGState(self.ctx)
    def restore(self): CGContextRestoreGState(self.ctx)
    def alpha(self, a): CGContextSetAlpha(self.ctx, a)

    def _fill(self, c):
        r, g, b = c[:3]; a = c[3] if len(c) > 3 else 1.0
        CGContextSetRGBFillColor(self.ctx, r, g, b, a)

    def _stroke(self, c, w):
        r, g, b = c[:3]; a = c[3] if len(c) > 3 else 1.0
        CGContextSetRGBStrokeColor(self.ctx, r, g, b, a)
        CGContextSetLineWidth(self.ctx, w)

    def dash(self, pattern, phase=0):
        if pattern:
            CGContextSetLineDash(self.ctx, phase, pattern, len(pattern))
        else:
            CGContextSetLineDash(self.ctx, 0, None, 0)

    def caps(self, round_=True):
        CGContextSetLineCap(self.ctx, kCGLineCapRound if round_ else kCGLineCapButt)
        CGContextSetLineJoin(self.ctx, kCGLineJoinRound)

    # -- primitives ----------------------------------------------------
    def rect(self, x, y, w, h, fill=None, stroke=None, lw=1):
        if fill:
            self._fill(fill); CGContextFillRect(self.ctx, CGRectMake(x, y, w, h))
        if stroke:
            self._stroke(stroke, lw)
            CGContextBeginPath(self.ctx); CGContextAddRect(self.ctx, CGRectMake(x, y, w, h))
            CGContextStrokePath(self.ctx)

    def rrect(self, x, y, w, h, r, fill=None, stroke=None, lw=1):
        r = min(r, w / 2, h / 2)
        CGContextBeginPath(self.ctx)
        CGContextMoveToPoint(self.ctx, x + r, y)
        CGContextAddLineToPoint(self.ctx, x + w - r, y)
        CGContextAddArc(self.ctx, x + w - r, y + r, r, -math.pi / 2, 0, 0)
        CGContextAddLineToPoint(self.ctx, x + w, y + h - r)
        CGContextAddArc(self.ctx, x + w - r, y + h - r, r, 0, math.pi / 2, 0)
        CGContextAddLineToPoint(self.ctx, x + r, y + h)
        CGContextAddArc(self.ctx, x + r, y + h - r, r, math.pi / 2, math.pi, 0)
        CGContextAddLineToPoint(self.ctx, x, y + r)
        CGContextAddArc(self.ctx, x + r, y + r, r, math.pi, 1.5 * math.pi, 0)
        CGContextClosePath(self.ctx)
        self._paint(fill, stroke, lw)

    def circle(self, cx, cy, r, fill=None, stroke=None, lw=1):
        rect = CGRectMake(cx - r, cy - r, 2 * r, 2 * r)
        if fill:
            self._fill(fill); CGContextFillEllipseInRect(self.ctx, rect)
        if stroke:
            self._stroke(stroke, lw); CGContextStrokeEllipseInRect(self.ctx, rect)

    def line(self, x1, y1, x2, y2, color, lw=1):
        self._stroke(color, lw)
        CGContextBeginPath(self.ctx)
        CGContextMoveToPoint(self.ctx, x1, y1)
        CGContextAddLineToPoint(self.ctx, x2, y2)
        CGContextStrokePath(self.ctx)

    def polyline(self, pts, stroke=None, lw=1, fill=None, close=False):
        if len(pts) < 2:
            return
        CGContextBeginPath(self.ctx)
        CGContextMoveToPoint(self.ctx, pts[0][0], pts[0][1])
        for x, y in pts[1:]:
            CGContextAddLineToPoint(self.ctx, x, y)
        if close:
            CGContextClosePath(self.ctx)
        self._paint(fill, stroke, lw)

    def bezier(self, p0, c1, c2, p1, stroke=None, lw=1):
        CGContextBeginPath(self.ctx)
        CGContextMoveToPoint(self.ctx, *p0)
        CGContextAddCurveToPoint(self.ctx, c1[0], c1[1], c2[0], c2[1], p1[0], p1[1])
        self._paint(None, stroke, lw)

    def arc(self, cx, cy, r, a0, a1, stroke=None, lw=1, fill=None, clockwise=0):
        CGContextBeginPath(self.ctx)
        CGContextAddArc(self.ctx, cx, cy, r, a0, a1, clockwise)
        self._paint(fill, stroke, lw)

    def _paint(self, fill, stroke, lw):
        if fill and stroke:
            self._fill(fill); self._stroke(stroke, lw)
            CGContextDrawPath(self.ctx, kCGPathFillStroke)
        elif fill:
            self._fill(fill); CGContextDrawPath(self.ctx, kCGPathFill)
        elif stroke:
            self._stroke(stroke, lw); CGContextDrawPath(self.ctx, kCGPathStroke)

    def clip_circle(self, cx, cy, r):
        CGContextBeginPath(self.ctx)
        CGContextAddEllipseInRect(self.ctx, CGRectMake(cx - r, cy - r, 2 * r, 2 * r))
        CGContextClip(self.ctx)

    def clip_poly(self, pts):
        CGContextBeginPath(self.ctx)
        CGContextMoveToPoint(self.ctx, pts[0][0], pts[0][1])
        for x, y in pts[1:]:
            CGContextAddLineToPoint(self.ctx, x, y)
        CGContextClosePath(self.ctx)
        CGContextClip(self.ctx)

    def arrowhead(self, x, y, ang, size, color, slim=0.42):
        """Filled triangle with its TIP at (x, y) pointing along ang."""
        b = (x - size * math.cos(ang), y - size * math.sin(ang))
        n = (-math.sin(ang), math.cos(ang))
        w = size * slim
        self.polyline([(x, y),
                       (b[0] + n[0] * w, b[1] + n[1] * w),
                       (b[0] - n[0] * w, b[1] - n[1] * w)],
                      fill=color, close=True)

    def arrow(self, x1, y1, x2, y2, color, lw=1, head=5):
        ang = math.atan2(y2 - y1, x2 - x1)
        bx, by = x2 - head * 0.85 * math.cos(ang), y2 - head * 0.85 * math.sin(ang)
        self.line(x1, y1, bx, by, color, lw)
        self.arrowhead(x2, y2, ang, head, color)

    # -- text ----------------------------------------------------------
    def _attr(self, s, role, size, color, tracking=0.0):
        rgba = tuple(color) + ((1.0,) if len(color) == 3 else ())
        at = {NSFontAttributeName: font(role, size),
              NSForegroundColorAttributeName:
                  NSColor.colorWithRed_green_blue_alpha_(*rgba)}
        if tracking:
            at[NSKernAttributeName] = tracking
        return NSAttributedString.alloc().initWithString_attributes_(s, at)

    def measure(self, s, role, size, tracking=0.0):
        return float(self._attr(s, role, size, (0, 0, 0), tracking).size()[0])

    def text(self, s, x, y, role='body', size=10, color=(0, 0, 0), align='l',
             tracking=0.0):
        """Draw `s` with its BASELINE at y; align in ('l','c','r')."""
        a = self._attr(s, role, size, color, tracking)
        w = float(a.size()[0])
        if align == 'c': x -= w / 2
        elif align == 'r': x -= w
        NSGraphicsContext.saveGraphicsState()
        NSGraphicsContext.setCurrentContext_(self._ns)
        a.drawAtPoint_((x, y + font(role, size).descender()))
        NSGraphicsContext.restoreGraphicsState()
        return w

    def para(self, lines, x, y, role='body', size=9.6, leading=13.2,
             color=(0, 0, 0), align='l', tracking=0.0):
        for i, ln in enumerate(lines):
            self.text(ln, x, y - i * leading, role, size, color, align, tracking)
        return y - (len(lines) - 1) * leading

    # -- maths ---------------------------------------------------------
    def _mruns(self, s, size, dy=0.0, scale=1.0, italic=True, out=None):
        """Parse a tiny TeX-ish markup into baseline-shifted runs.

        `_{..}` `^{..}` nest and shrink; `\\rm{..}` / `\\it{..}` force upright or
        italic. Bare ASCII letters are italic by default (maths convention);
        digits and operators stay upright.
        """
        out = [] if out is None else out
        s = s.replace('\\,', '\u2009').replace('\\;', '\u2005')
        i, buf = 0, ''
        def flush():
            nonlocal buf
            if buf:
                out.append((buf, dy, scale, italic)); buf = ''
        while i < len(s):
            c = s[i]
            if c in '^_':
                flush(); i += 1
                if i < len(s) and s[i] == '{':
                    d, j = 1, i + 1
                    while j < len(s) and d:
                        d += (s[j] == '{') - (s[j] == '}'); j += 1
                    grp, i = s[i + 1:j - 1], j
                else:
                    grp, i = s[i], i + 1
                shift = (0.44 if c == '^' else -0.21) * size * scale
                self._mruns(grp, size, dy + shift, scale * 0.74, italic, out)
            elif c == '\\' and s[i + 1:i + 4] in ('rm{', 'it{'):
                flush()
                mode = s[i + 1:i + 3]; d, j = 1, i + 4
                while j < len(s) and d:
                    d += (s[j] == '{') - (s[j] == '}'); j += 1
                self._mruns(s[i + 4:j - 1], size, dy, scale, mode == 'it', out)
                i = j
            else:
                buf += c; i += 1
        flush()
        return out

    def _msplit(self, runs):
        """Split each run so letters get italic and non-letters upright."""
        LET = set('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
                  'αβγδεζηθικλμνξπρστυφχψωΓΔΘΛΞΠΣΦΨΩ')
        out = []
        for txt, dy, sc, it in runs:
            cur, curit = '', None
            for ch in txt:
                want = it and (ch in LET)
                if curit is None or want == curit:
                    cur += ch; curit = want
                else:
                    out.append((cur, dy, sc, curit)); cur, curit = ch, want
            if cur:
                out.append((cur, dy, sc, curit))
        return out

    def math_width(self, s, size):
        runs = self._msplit(self._mruns(s, size))
        return sum(self.measure(t, 'mathi' if it else 'math', size * sc)
                   for t, _, sc, it in runs)

    def math(self, s, x, y, size=11, color=(0, 0, 0), align='l'):
        """Typeset tiny-TeX `s` with its main baseline at y."""
        runs = self._msplit(self._mruns(s, size))
        w = self.math_width(s, size)
        if align == 'c': x -= w / 2
        elif align == 'r': x -= w
        for t, dy, sc, it in runs:
            role = 'mathi' if it else 'math'
            x += self.text(t, x, y + dy, role, size * sc, color)
        return w

    # -- line breaking, with inline maths --------------------------------
    # A $...$ span inside a paragraph is one atomic token, typeset by math()
    # rather than as body text. Body faces lack the Unicode sub/superscript
    # letters, so writing "z_j" as a literal glyph font-falls-back and spaces
    # wrongly; this routes it through STIXGeneral instead.
    def _tokens(self, s):
        out, parts = [], s.split('$')
        for k, part in enumerate(parts):
            if k % 2:
                pre = parts[k - 1]
                out.append([part, True, bool(out) and (pre == '' or pre[-1:].isspace())])
            else:
                words = part.split()
                for wi, w in enumerate(words):
                    out.append([w, False, True if wi else
                                (bool(out) and part[:1].isspace())])
        return out

    def _tokw(self, tok, role, size, tracking):
        return (self.math_width(tok[0], size) if tok[1]
                else self.measure(tok[0], role, size, tracking))

    def _spacew(self, role, size, tracking):
        return (self.measure("n n", role, size, tracking)
                - 2 * self.measure("n", role, size, tracking))

    def wrap(self, s, role, size, width, tracking=0.0):
        """Greedy wrap to `width` points, measured in the real font(s).
        Returns a list of token-lines; use flow() to draw them."""
        sw = self._spacew(role, size, tracking)
        lines, cur, x = [], [], 0.0
        for tok in self._tokens(s):
            w = self._tokw(tok, role, size, tracking)
            adv = (sw if (tok[2] and cur) else 0.0) + w
            if cur and x + adv > width:
                lines.append(cur); cur, x = [tok], w
            else:
                cur.append(tok); x += adv
        if cur:
            lines.append(cur)
        return lines

    def _draw_line(self, toks, x, y, role, size, color, tracking, gap):
        for k, tok in enumerate(toks):
            if k and tok[2]:
                x += gap
            if tok[1]:
                x += self.math(tok[0], x, y, size, color)
            else:
                x += self.text(tok[0], x, y, role, size, color, 'l', tracking)
        return x

    def flow_lines(self, lines, x, y, width, role='body', size=10, leading=14.5,
                   color=(0, 0, 0), tracking=0.0, just=True, justify_last=False):
        """Draw pre-wrapped token-lines; returns the baseline y after the last."""
        sw = self._spacew(role, size, tracking)
        for i, toks in enumerate(lines):
            gap = sw
            if just and (justify_last or i < len(lines) - 1):
                n = sum(1 for k, t in enumerate(toks) if k and t[2])
                if n:
                    solid = sum(self._tokw(t, role, size, tracking) for t in toks)
                    gap = max(sw * 0.72, (width - solid) / n)
            self._draw_line(toks, x, y - i * leading, role, size, color, tracking, gap)
        return y - len(lines) * leading

    def flow(self, s, x, y, width, role='body', size=10, leading=14.5,
             color=(0, 0, 0), tracking=0.0, just=True):
        """Wrap and draw a paragraph; returns the baseline y AFTER the last line."""
        return self.flow_lines(self.wrap(s, role, size, width, tracking), x, y,
                               width, role, size, leading, color, tracking, just)

    def flow_columns(self, s, xs, y, width, role='body', size=10, leading=14.5,
                     color=(0, 0, 0), tracking=0.0):
        """Balance one paragraph across len(xs) columns of equal width, so the
        columns end level instead of one running long."""
        lines = self.wrap(s, role, size, width, tracking)
        n = len(xs)
        per = -(-len(lines) // n)
        for k, x in enumerate(xs):
            chunk = lines[k * per:(k + 1) * per]
            if chunk:
                self.flow_lines(chunk, x, y, width, role, size, leading, color,
                                tracking, just=True,
                                justify_last=(k < n - 1 and len(chunk) == per))
        return y - per * leading

    def flow_height(self, s, role, size, leading, width, tracking=0.0):
        return len(self.wrap(s, role, size, width, tracking)) * leading

    # -- general path + bitmap ----------------------------------------
    def bezpath(self, cmds, fill=None, stroke=None, lw=1):
        """cmds: ('m',x,y) ('l',x,y) ('c',x1,y1,x2,y2,x,y) ('z',)"""
        CGContextBeginPath(self.ctx)
        for cm in cmds:
            k = cm[0]
            if k == 'm':   CGContextMoveToPoint(self.ctx, cm[1], cm[2])
            elif k == 'l': CGContextAddLineToPoint(self.ctx, cm[1], cm[2])
            elif k == 'c': CGContextAddCurveToPoint(self.ctx, *cm[1:])
            elif k == 'z': CGContextClosePath(self.ctx)
        self._paint(fill, stroke, lw)

    def image_rgba(self, data, w, h, x, y, dw, dh):
        """Draw a premultiplied-RGBA bytes buffer. Row 0 of `data` is the TOP."""
        from Quartz import (CGImageCreate, CGColorSpaceCreateDeviceRGB,
                            CGDataProviderCreateWithCFData, CGContextDrawImage,
                            kCGImageAlphaPremultipliedLast,
                            kCGRenderingIntentDefault, CGRectMake as _R)
        prov = CGDataProviderCreateWithCFData(bytes(data))
        img = CGImageCreate(w, h, 8, 32, w * 4, CGColorSpaceCreateDeviceRGB(),
                            kCGImageAlphaPremultipliedLast, prov, None, True,
                            kCGRenderingIntentDefault)
        CGContextDrawImage(self.ctx, _R(x, y, dw, dh), img)
