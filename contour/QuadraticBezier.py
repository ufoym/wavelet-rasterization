import numpy as np
from collections import namedtuple
from util.solver import quadric

# -----------------------------------------------------------------------------
Point = namedtuple('Point', 'x y')
# -----------------------------------------------------------------------------

class QuadraticBezier:

    def __init__(self, x0, y0, x1, y1, x2, y2):
        self.x0, self.y0, self.x1, self.y1 = x0, y0, x1, y1
        self.x2, self.y2 = x2, y2

    def evaluate(self, t):
        return (self.x0*(1-t)**2 + 2*self.x1*(1-t)*t + self.x2*t**2,
                self.y0*(1-t)**2 + 2*self.y1*(1-t)*t + self.y2*t**2)

    def subsection(self, t0, t1):
        x0, y0 = self.evaluate(t0)
        x2, y2 = self.evaluate(t1)
        u0 = 1.0 - t0
        xm = u0 * self.x1 + t0 * self.x2
        ym = u0 * self.y1 + t0 * self.y2
        t  = (t1 - t0) / u0
        x1 = (1-t) * x0 + t * xm
        y1 = (1-t) * y0 + t * ym

        sec = QuadraticBezier( x0, y0, x1, y1, x2, y2 )
        return sec

    def clip(self, left, right, bottom, top):
        def is_t_in(t, eps = 1e-5):
            pt = self.evaluate(t)
            return left-eps<=pt[0]<=right+eps and top-eps<=pt[1]<=bottom+eps

        ax = self.x0 - 2*self.x1 + self.x2
        bx = -2*self.x0 + 2*self.x1
        _cx = self.x0
        ay = self.y0 - 2*self.y1 + self.y2
        by = -2*self.y0 + 2*self.y1
        _cy = self.y0
        ts = [0]
        ts += quadric(ax, bx, _cx-left)
        ts += quadric(ax, bx, _cx-right)
        ts += quadric(ay, by, _cy-bottom)
        ts += quadric(ay, by, _cy-top)
        ts.append(1)
        ts = [t for t in ts if 0 <= t <= 1 and is_t_in(t)]
        ts = sorted(ts)
        ts = [t for i, t in enumerate(ts) if t != ts[i-1]]
        pairs = [(ts[i-1], t) for i, t in enumerate(ts) \
                    if i > 0 and is_t_in((t + ts[i-1]) * 0.5)]
        sections = [self.subsection(a, b) for a, b in pairs]
        return sections

    def get_KL(self, eps = 1e-5):
        Kx, Ky, Lx, Ly = 0, 0, 0, 0
        for sec in self.clip(0, 1, 1, 0):
            v2 = Point(sec.x0, sec.y0)
            v1 = Point(sec.x1, sec.y1)
            v0 = Point(sec.x2, sec.y2)
            if abs(v0.x-1) < eps and abs(v1.x-1) < eps and abs(v2.x-1) < eps\
            or abs(v0.y-1) < eps and abs(v1.y-1) < eps and abs(v2.y-1) < eps:
                continue

            Kx += 1./4 * (v0.y - v2.y)
            Ky += 1./4 * (v2.x - v0.x)
            Lx += 1./24* (3 * v0.x*v0.y + 2 * v0.y*v1.x - 2 * v0.x*v1.y \
                        + v0.y*v2.x + 2 * v1.y*v2.x -(v0.x+2*v1.x+3*v2.x)*v2.y)
            Ly += 1./24* (2 * v1.y*v2.x + v0.y * (2*v1.x+v2.x) - 2 * v1.x*v2.y\
                        + 3 * v2.x*v2.y - v0.x * (3*v0.y+2*v1.y+v2.y))
        return Point(Kx, Ky), Point(Lx, Ly)

# -----------------------------------------------------------------------------

class Contour:

    def __init__(self, contour):
        self.contour = contour

    def __str__(self):
        info = ' '.join('(%2.1f, %2.1f)' % (s[0], s[1]) for s in self.contour)
        return ' :\t'.join(['QuadraticBezier', info])

    def process(self, method):
        self.contour = [method(p) for p in self.contour]

    def area(self):
        def det(a, b):  return a[0] * b[1] - a[1] * b[0]
        s = 0
        for v0, v1, v2 in self.each():
            s += 1./3 * det(v0,v1) + 1./3 * det(v1,v2) + 1./6 * det(v0,v2)
        return s

    def each(self): 
        for i in xrange(0, len(self.contour), 2):
            v2 = self.contour[i]
            v1 = self.contour[i-1]
            v0 = self.contour[i-2]
            yield v0, v1, v2

    def to_lines(self):
        tts = np.linspace(0,1, num=50)
        for section in self.each():
            section = (s[i] for s in section for i in xrange(2))
            bezier = QuadraticBezier(*section)
            for start, end in zip(tts[:-1], tts[1:]):
                sx, sy = bezier.evaluate(start)
                ex, ey = bezier.evaluate(end)
                yield sx, sy, ex, ey

    def get_KL(self, section):
        bezier = QuadraticBezier(*section)
        return bezier.get_KL()

# -----------------------------------------------------------------------------