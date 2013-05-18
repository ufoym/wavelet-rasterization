import Polygon, numpy as np, math
from collections import namedtuple

# -----------------------------------------------------------------------------
Point = namedtuple('Point', 'x y')

class Rasterizer:

    def __init__(self, poly, w, h):
        def normalize(p):
            return (p[0]/float(w), p[1]/float(h))
        box = [(0,0), (1,0), (1,1), (0,1)]
        lattice = [Point(x,y) for x in xrange(w) for y in xrange(h)]
        self._poly = Polygon.Polygon([normalize(p) for p in poly])
        self._box = [Point(*normalize(p)) for p in box]
        self._lattice = [Point(*normalize(p)) for p in lattice]
        self._w = w
        self._h = h

    def _psi(self, p, e, j, k):
        def psi_1d(p, e):
            if e == 0:   return 1 if 0 <= p < 1 else 0
            else:        return (1 if 0<=p<0.5 else -1) if 0 <= p < 1 else 0 
        return 2**j * psi_1d(2**j*p.x-k.x, e.x) * psi_1d(2**j*p.y-k.y, e.y)

    def _c(self, e, j, k):
        s = 0
        for p in self._lattice:
            psi_p = self._psi(p, e, j, k)
            if psi_p != 0:
                box = Polygon.Polygon([(b.x+p.x, b.y+p.y) for b in self._box])
                px_poly = box & self._poly
                s += px_poly.area() * psi_p
        return s

    def _g(self, p):
        s = 0
        E_pi = [Point(0,0), Point(0,1), Point(1,0), Point(1,1)]
        # ---------------------------------------------------------------------
        j = 0
        for kx in xrange(2**j):
            for ky in xrange(2**j):
                k = Point(kx, ky)
                s += self._c(E_pi[0], j,k) * self._psi(p, E_pi[0], j,k)
        # ---------------------------------------------------------------------
        max_j = int(math.ceil(math.log(max(w,h), 2)))-1
        for j in xrange(max_j+1):
            for kx in xrange(2**j):
                for ky in xrange(2**j):
                    k = Point(kx, ky)
                    for e in E_pi[1:]:
                        s += self._c(e, j, k) * self._psi(p, e, j, k)
        return s

    def get(self):
        px_arr = [self._g(p) for p in self._lattice]
        return np.reshape(px_arr, (self._w, self._h))

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    import cv2, time
    from random import randint
    w, h, z = 8, 8, 30

    while True:
        poly = [(randint(1,w-1), randint(1,w-1)) for i in xrange(3)]

        ts = time.time()
        raster = Rasterizer(poly, w, h).get()
        print '%s\ttime: %2.1fs' % (poly, time.time()-ts)

        raster = np.array(raster*255+0.5, np.uint8)
        raster = cv2.cvtColor(raster, cv2.COLOR_GRAY2BGR)
        raster = cv2.resize(raster, (w*z,h*z), interpolation=cv2.INTER_NEAREST)

        poly_vis = np.array([(p[1]*z, p[0]*z) for p in poly])
        cv2.polylines(raster, [poly_vis],True,(255,255,0),lineType=cv2.CV_AA)
        cv2.namedWindow('raster')
        cv2.imshow('raster', raster)

        cv2.waitKey(0)