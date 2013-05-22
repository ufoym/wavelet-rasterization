import Polygon, numpy as np, math
from collections import namedtuple
from util.LiangBarsky import LiangBarsky

# -----------------------------------------------------------------------------
Point = namedtuple('Point', 'x y')

class Rasterizer:

    def __init__(self, poly, w, h):
        def normalize(p):
            return (p[0]/float(w), p[1]/float(h))
        box = [(0,0), (1,0), (1,1), (0,1)]
        lattice = [Point(x,y) for x in xrange(w) for y in xrange(h)]
        self._poly = Polygon.Polygon([normalize(p) for p in poly])
        self._normed_poly = [normalize(p) for p in poly]
        self._box = [Point(*normalize(p)) for p in box]
        self._lattice = [Point(*normalize(p)) for p in lattice]
        self._w = w
        self._h = h
        self._all_c = {}
        self._get_all_c()

    def _psi(self, p, e, j, k):
        def psi_1d(p, e):
            if e == 0:   return 1 if 0 <= p < 1 else 0
            else:        return (1 if 0<=p<0.5 else -1) if 0 <= p < 1 else 0 
        return 2**j * psi_1d(2**j*p.x-k.x, e.x) * psi_1d(2**j*p.y-k.y, e.y)

    def _c(self, j, k):
        def get_KL(_start, _end, Q):
            def transform(p):
                return Point(2**(j+1)*p[0]-k.x*2-Q.x, 
                             2**(j+1)*p[1]-k.y*2-Q.y)
            start = transform(_start)
            end = transform(_end)
            try:   
                v1, v0 = LiangBarsky(0, 1, 1, 0, start, end)
            except:
                return Point(0,0), Point(0,0)
            # compute K, L
            K = Point(0.25 * (v0.y-v1.y), 
                      0.25 * (v1.x-v0.x))
            L = Point(0.125* (v0.y-v1.y) * (v0.x+v1.x), 
                      0.125* (v1.x-v0.x) * (v0.y+v1.y))
            return K, L
        Q_00, Q_01 = Point(0, 0), Point(0, 1)
        Q_10, Q_11 = Point(1, 0), Point(1, 1)
        c10, c01, c11 = 0, 0, 0
        for i, v1 in enumerate(self._normed_poly):
            v0 = self._normed_poly[i-1]
            KQ00, LQ00 = get_KL(v0, v1, Q_00)
            KQ01, LQ01 = get_KL(v0, v1, Q_01)
            KQ10, LQ10 = get_KL(v0, v1, Q_10)
            KQ11, LQ11 = get_KL(v0, v1, Q_11)
            c10 += LQ00.x + LQ01.x + KQ10.x \
                 - LQ10.x + KQ11.x - LQ11.x
            c01 += LQ00.y + LQ10.y + KQ01.y \
                 - LQ01.y + KQ11.y - LQ11.y
            c11 += LQ00.x - LQ01.x + KQ10.x \
                 - LQ10.x - KQ11.x + LQ11.x
        return c01, c10, c11

    def _c00(self):
        s = 0
        e = Point(0,0)
        j = 0
        k = Point(0,0)
        for p in self._lattice:
            psi_p = self._psi(p, e, j, k)
            if psi_p != 0:
                box = Polygon.Polygon([(b.x+p.x, b.y+p.y) for b in self._box])
                px_poly = box & self._poly
                s += px_poly.area() * psi_p
        return s

    def _get_all_c(self):
        max_wh = max(w,h)
        max_j = int(math.ceil(math.log(max_wh, 2)))-1
        for j in xrange(max_j+1):
            for kx in xrange(2**j):
                for ky in xrange(2**j):
                    k = Point(kx, ky)
                    self._all_c[(j,k)] = self._c(j, k)

    def _g(self, p):
        s = 0
        E_pi = [Point(0,0), Point(0,1), Point(1,0), Point(1,1)]
        # ---------------------------------------------------------------------
        s += self._c00() * self._psi(p, E_pi[0], 0, Point(0,0))
        # ---------------------------------------------------------------------
        max_wh = max(w,h)
        max_j = int(math.ceil(math.log(max_wh, 2)))-1
        for j in xrange(max_j+1):
            for kx in xrange(2**j):
                for ky in xrange(2**j):
                    k = Point(kx, ky)
                    cs = self._all_c[(j,k)]
                    for i, e in enumerate(E_pi[1:]):
                        psi = self._psi(p, e, j, k)
                        if psi > 0: sign = 1
                        elif psi < 0: sign = -1
                        else: sign = 0
                        s += cs[i] * sign
        return s

    def get(self):
        px_arr = [self._g(p) for p in self._lattice]
        return np.reshape(px_arr, (self._w, self._h))

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    import cv2, time
    from random import randint
    w, h, z = 8, 8, 30

    def area(poly):
        return 0.5 * sum(x0*y1 - x1*y0 \
            for ((x0, y0), (x1, y1)) in zip(poly, poly[1:] + [poly[0]]))

    while True:
        poly = [(randint(1,w-1), randint(1,w-1)) for i in xrange(3)]
        if area(poly) < 0:
            continue

        ts = time.time()
        raster = Rasterizer(poly, w, h).get()
        print '%s\ttime: %2.1fs' % (poly, time.time()-ts)

        raster = np.array(raster*255+0.5, np.uint8)
        raster = cv2.cvtColor(raster, cv2.COLOR_GRAY2BGR)
        raster = cv2.resize(raster, (w*z,h*z), interpolation=cv2.INTER_NEAREST)

        poly_vis = np.array([(int(p[1]*z), int(p[0]*z)) for p in poly])
        cv2.polylines(raster, [poly_vis],True,(255,255,0),lineType=cv2.CV_AA)
        cv2.namedWindow('raster')
        cv2.imshow('raster', raster)

        cv2.waitKey(0)