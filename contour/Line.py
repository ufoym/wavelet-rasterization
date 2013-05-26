import numpy as np
from collections import namedtuple

# -----------------------------------------------------------------------------
Point = namedtuple('Point', 'x y')
# -----------------------------------------------------------------------------

class Line:

    def __init__(self, x0, y0, x1, y1):
        self.x0, self.y0, self.x1, self.y1 = x0, y0, x1, y1

    def clip(self, left, right, bottom, top):
        t0, t1 = 0, 1
        xdelta = self.x1 - self.x0
        ydelta = self.y1 - self.y0
        for edge in xrange(4):#traverse through left, right, bottom, top edges.
            if   edge == 0:   p, q = -xdelta, -(left-self.x0) 
            elif edge == 1:   p, q =  xdelta,  (right-self.x0)
            elif edge == 2:   p, q =  ydelta,  (bottom-self.y0)
            elif edge == 3:   p, q = -ydelta, -(top-self.y0)
            if p == 0 and q < 0:    return []
            if p < 0:
                r = q / float(p)
                if r > t1:          return []
                elif r > t0:        t0 = r   # line is clipped!
            elif p > 0:
                r = q / float(p)
                if r < t0:          return []
                elif r < t1:        t1 = r   # line is clipped!
        clipped_line = Line(self.x0 + t0*xdelta, self.y0 + t0*ydelta, 
                            self.x0 + t1*xdelta, self.y0 + t1*ydelta)
        return [clipped_line]

    def get_KL(self, eps = 1e-5):
        Kx, Ky, Lx, Ly = 0, 0, 0, 0
        for sec in self.clip(0, 1, 1, 0):
            v1 = Point(sec.x0, sec.y0)
            v0 = Point(sec.x1, sec.y1)
            if abs(v0.x - 1) < eps and abs(v1.x - 1) < eps \
            or abs(v0.y - 1) < eps and abs(v1.y - 1) < eps:
                continue

            Kx += 1./4 * (v0.y - v1.y)
            Ky += 1./4 * (v1.x - v0.x)
            Lx += 1./8 * (v0.y-v1.y) * (v0.x+v1.x)
            Ly += 1./8 * (v1.x-v0.x) * (v0.y+v1.y)
        return Point(Kx, Ky), Point(Lx, Ly)

# -----------------------------------------------------------------------------

class Contour:

    def __init__(self, contour):
        self.contour = contour

    def __str__(self):
        info = ' '.join('(%2.1f, %2.1f)' % (s[0], s[1]) for s in self.contour)
        return ' :\t'.join(['Line', info])

    def process(self, method):
        self.contour = [method(p) for p in self.contour]

    def area(self):
        def det(a, b):  return a[0] * b[1] - a[1] * b[0]
        return 0.5 * sum(det(a, b) for (a, b) in \
            zip(self.contour, self.contour[1:] + [self.contour[0]]))

    def each(self): 
        for i, v1 in enumerate(self.contour):
            v0 = self.contour[i-1]
            yield v0, v1

    def to_lines(self):
        for v0, v1 in self.each():
            yield v0[0], v0[1], v1[0], v1[1]

    def get_KL(self, section):
        line = Line(*section)
        return line.get_KL()

# -----------------------------------------------------------------------------