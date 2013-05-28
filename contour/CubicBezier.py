import numpy as np
from collections import namedtuple
from util.solver import cubic

# -----------------------------------------------------------------------------
Point = namedtuple('Point', 'x y')
# -----------------------------------------------------------------------------

class CubicBezier:

    def __init__(self, x0, y0, x1, y1, x2, y2, x3, y3):
        self.x0, self.y0, self.x1, self.y1 = x0, y0, x1, y1
        self.x2, self.y2, self.x3, self.y3 = x2, y2, x3, y3

    def evaluate(self, t):
        return (self.x0*(1-t)**3 + 3*self.x1*(1-t)**2*t \
                + 3*self.x2*(1-t)*t**2 + self.x3*t**3,
                self.y0*(1-t)**3 + 3*self.y1*(1-t)**2*t \
                + 3*self.y2*(1-t)*t**2 + self.y3*t**3)

    def subsection(self, t0, t1):
        u0 = 1.0 - t0
        u1 = 1.0 - t1

        qxa = self.x0*u0*u0 + self.x1*2*t0*u0 + self.x2*t0*t0
        qxb = self.x0*u1*u1 + self.x1*2*t1*u1 + self.x2*t1*t1
        qxc = self.x1*u0*u0 + self.x2*2*t0*u0 + self.x3*t0*t0
        qxd = self.x1*u1*u1 + self.x2*2*t1*u1 + self.x3*t1*t1

        qya = self.y0*u0*u0 + self.y1*2*t0*u0 + self.y2*t0*t0
        qyb = self.y0*u1*u1 + self.y1*2*t1*u1 + self.y2*t1*t1
        qyc = self.y1*u0*u0 + self.y2*2*t0*u0 + self.y3*t0*t0
        qyd = self.y1*u1*u1 + self.y2*2*t1*u1 + self.y3*t1*t1

        sec = CubicBezier(  qxa*u0 + qxc*t0, qya*u0 + qyc*t0,
                            qxa*u1 + qxc*t1, qya*u1 + qyc*t1,
                            qxb*u0 + qxd*t0, qyb*u0 + qyd*t0,
                            qxb*u1 + qxd*t1, qyb*u1 + qyd*t1 )
        return sec

    def clip(self, left, right, bottom, top):
        def is_t_in(t, eps = 1e-5):
            pt = self.evaluate(t)
            return left-eps<=pt[0]<=right+eps and top-eps<=pt[1]<=bottom+eps

        ax = -self.x0 + 3*self.x1 - 3*self.x2 + self.x3
        bx = 3*self.x0 - 6*self.x1 + 3*self.x2
        cx, _dx = 3*self.x1 - 3*self.x0, self.x0
        ay = -self.y0 + 3*self.y1 - 3*self.y2 + self.y3
        by = 3*self.y0 - 6*self.y1 + 3*self.y2
        cy, _dy = 3*self.y1 - 3*self.y0, self.y0
        ts = [0]
        ts += cubic(ax, bx, cx, _dx-left)
        ts += cubic(ax, bx, cx, _dx-right)
        ts += cubic(ay, by, cy, _dy-bottom)
        ts += cubic(ay, by, cy, _dy-top)
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
            v3 = Point(sec.x0, sec.y0)
            v2 = Point(sec.x1, sec.y1)
            v1 = Point(sec.x2, sec.y2)
            v0 = Point(sec.x3, sec.y3)
            if abs(v0.x - 1) < eps and abs(v1.x - 1) < eps \
                and abs(v2.x - 1) < eps and abs(v3.x - 1) < eps\
            or abs(v0.y - 1) < eps and abs(v1.y - 1) < eps \
                and abs(v2.y - 1) < eps and abs(v3.y - 1) < eps:
                continue

            Kx += 1./4 * (v0.y - v3.y)
            Ky += 1./4 * (v3.x - v0.x)
            Lx += 1./80* (6 * v2.y*v3.x + 3 * v1.y*(v2.x+v3.x) \
                        + v0.y * (6*v1.x+3*v2.x+v3.x) \
                        - 6 * v2.x*v3.y - 3 * v1.x*(v2.y+v3.y) \
                        - 10 * v3.x*v3.y \
                        + v0.x * (10*v0.y-6*v1.y-3*v2.y-v3.y) )
            Ly += 1./80* (6 * v2.y*v3.x + 3 * v1.y*(v2.x+v3.x) \
                        + v0.y * (6*v1.x+3*v2.x+v3.x) \
                        - 6 * v2.x*v3.y - 3 * v1.x*(v2.y+v3.y) \
                        + 10 * v3.x*v3.y \
                        - v0.x * (10*v0.y+6*v1.y+3*v2.y+v3.y) )
        return Point(Kx, Ky), Point(Lx, Ly)

# -----------------------------------------------------------------------------

class Contour:

    def __init__(self, contour):
        self.contour = contour

    def __str__(self):
        info = ' '.join('(%2.1f, %2.1f)' % (s[0], s[1]) for s in self.contour)
        return ' :\t'.join(['CubicBezier', info])

    def process(self, method):
        self.contour = [method(p) for p in self.contour]

    def area(self):
        def det(a, b):  return a[0] * b[1] - a[1] * b[0]
        s = 0
        for v0, v1, v2, v3 in self.each():
            s += 3./10 * det(v0,v1) + 3./20 * det(v1,v2) + 3./10 * det(v2,v3) \
               + 3./20 * det(v0,v2) + 3./20 * det(v1,v3) + 1./20 * det(v0,v3)
        return s

    def each(self): 
        for i in xrange(0, len(self.contour), 3):
            v3 = self.contour[i]
            v2 = self.contour[i-1]
            v1 = self.contour[i-2]
            v0 = self.contour[i-3]
            yield v0, v1, v2, v3

    def to_lines(self):
        tts = np.linspace(0,1, num=50)
        for section in self.each():
            section = (s[i] for s in section for i in xrange(2))
            bezier = CubicBezier(*section)
            for start, end in zip(tts[:-1], tts[1:]):
                sx, sy = bezier.evaluate(start)
                ex, ey = bezier.evaluate(end)
                yield sx, sy, ex, ey

    def get_KL(self, section):
        bezier = CubicBezier(*section)
        return bezier.get_KL()

# -----------------------------------------------------------------------------