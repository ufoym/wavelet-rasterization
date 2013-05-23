import math, collections
from util.LiangBarsky import LiangBarsky

# -----------------------------------------------------------------------------
Point = collections.namedtuple('Point', 'x y')
def area(poly):
    return 0.5 * sum(x0*y1 - x1*y0 \
        for ((x0, y0), (x1, y1)) in zip(poly, poly[1:] + [poly[0]]))

class Rasterizer:

    def __init__(self, poly, w, h):
        self._w = w
        self._h = h
        self._max_j = int(math.ceil(math.log(max(w,h), 2)))-1
        self._wh = 2 ** (self._max_j+1)
        def normalize(p):
            return (p[0]/float(self._wh), p[1]/float(self._wh))
        self._poly = [normalize(p) for p in poly]
        self._area = area(self._poly)
        self._lattice = [Point(*normalize((x,y))) \
                        for x in xrange(self._wh) for y in xrange(self._wh)]
        self._valid_lattice = [Point(*normalize((x,y))) \
                        for x in xrange(h) for y in xrange(w)]
        # prepare all c
        self._all_c = {}
        for j in xrange(self._max_j+1):
            for kx in xrange(2**j):
                for ky in xrange(2**j):
                    k = Point(kx, ky)
                    self._all_c[(j,k)] = self._c(j, k)

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
                if v0.x == 1 and v1.x == 1 or v0.y == 1 and v1.y == 1:
                    return Point(0,0), Point(0,0)                    
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
        for i, v1 in enumerate(self._poly):
            v0 = self._poly[i-1]
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

    def _g(self, p):
        s = self._area
        E = [Point(0,1), Point(1,0), Point(1,1)]
        for j in xrange(self._max_j+1):
            for kx in xrange(2**j):
                for ky in xrange(2**j):
                    k = Point(kx, ky)
                    cs = self._all_c[(j,k)]
                    for i, e in enumerate(E):
                        psi = self._psi(p, e, j, k)
                        if psi > 0: sign = 1
                        elif psi < 0: sign = -1
                        else: sign = 0
                        s += cs[i] * sign
        return s

    def get(self):
        px_arr = [self._g(p) for p in self._valid_lattice]
        px_mat = [px_arr[i*self._w : (i+1)*self._w] for i in xrange(self._h)]
        return px_mat

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    import cv2, time, numpy as np
    from random import randint
    w, h, z = 13, 17, 30

    while True:
        poly = [(randint(1,h-1), randint(1,w-1)) for i in xrange(3)]
        if area(poly) < 0:
            continue

        ts = time.time()
        raster = Rasterizer(poly, w, h).get()
        print '%s\ttime: %2.1fs' % (poly, time.time()-ts)

        raster = np.array(np.asarray(raster)*255+0.5, np.uint8)
        raster = cv2.cvtColor(raster, cv2.COLOR_GRAY2BGR)
        raster = cv2.resize(raster, (w*z,h*z), interpolation=cv2.INTER_NEAREST)

        poly_vis = np.array([(int(p[1]*z), int(p[0]*z)) for p in poly])
        cv2.polylines(raster, [poly_vis],True,(255,255,0),lineType=cv2.CV_AA)
        cv2.namedWindow('raster')
        cv2.imshow('raster', raster)

        cv2.waitKey(100)