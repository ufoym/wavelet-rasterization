import math, copy
from Line import Point, Contour, get_KL

# -----------------------------------------------------------------------------

class Rasterizer:

    def __init__(self, contour, w, h):
        self._w = w
        self._h = h
        self._max_j = int(math.ceil(math.log(max(w,h), 2)))-1
        self._wh = 2 ** (self._max_j+1)
        def normalize(p):   return (p[0]/float(self._wh), p[1]/float(self._wh))
        self._contour = copy.deepcopy(contour)
        self._contour.process(normalize)
        self._area = self._contour.area()
        self._lattice = [Point(*normalize((x,y))) \
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
        def transform(section, Q):
            return (2**(j+1)*p[i]-k[i]*2-Q[i] \
                    for p in section for i in xrange(2))
        Q_00, Q_01 = Point(0, 0), Point(0, 1)
        Q_10, Q_11 = Point(1, 0), Point(1, 1)
        c10, c01, c11 = 0, 0, 0
        for section in self._contour.each():
            KQ00, LQ00 = get_KL(transform(section, Q_00))
            KQ01, LQ01 = get_KL(transform(section, Q_01))
            KQ10, LQ10 = get_KL(transform(section, Q_10))
            KQ11, LQ11 = get_KL(transform(section, Q_11))
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
        px_arr = [self._g(p) for p in self._lattice]
        px_mat = [px_arr[i*self._w : (i+1)*self._w] for i in xrange(self._h)]
        return px_mat

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    import cv2, time, numpy as np
    from random import randint
    w, h, z = 17, 13, 30

    while True:
        contour = Contour([(randint(1,h-1),randint(1,w-1)) for i in xrange(3)])
        if contour.area() < 0:
            continue

        ts = time.time()
        raster = Rasterizer(contour, w, h).get()
        print '%s\ttime: %2.1fs' % (contour, time.time()-ts)

        raster = np.array(np.asarray(raster)*255+0.5, np.uint8)
        raster = cv2.cvtColor(raster, cv2.COLOR_GRAY2BGR)
        raster = cv2.resize(raster, (w*z,h*z), interpolation=cv2.INTER_NEAREST)
        for sx, sy, ex, ey in contour.to_lines():
            cv2.line(raster, (int(sy*z),int(sx*z)), (int(ey*z),int(ex*z)), 
                (0,0,255), 2, cv2.CV_AA)
        cv2.namedWindow('raster')
        cv2.imshow('raster', raster)

        cv2.waitKey(1)