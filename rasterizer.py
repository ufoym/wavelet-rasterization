import math, copy
from collections import namedtuple

# -----------------------------------------------------------------------------
Point = namedtuple('Point', 'x y')
# -----------------------------------------------------------------------------

class Rasterizer:

    def __init__(self, contour, w, h):
        self.w = w
        self.h = h
        self.max_j = int(math.ceil(math.log(max(w,h), 2)))-1
        self.wh = 2 ** (self.max_j+1)
        def normalize(p):   return (p[0]/float(self.wh), p[1]/float(self.wh))
        self.contour = copy.deepcopy(contour)
        self.contour.process(normalize)
        self.area = self.contour.area()
        self.lattice = [Point(*normalize((x,y))) \
                        for x in xrange(h) for y in xrange(w)]
        # prepare all c
        self.all_c = {}
        for j in xrange(self.max_j+1):
            for kx in xrange(2**j):
                for ky in xrange(2**j):
                    self.all_c[(j, kx,ky)] = self.c(j, (kx,ky))

    def psi(self, p, e, j, k):
        def psi_1d(p, e):
            if e == 0:   return 1 if 0 <= p < 1 else 0
            else:        return (1 if 0<=p<0.5 else -1) if 0 <= p < 1 else 0 
        return 2**j * psi_1d(2**j*p.x-k.x, e.x) * psi_1d(2**j*p.y-k.y, e.y)

    def c(self, j, k):
        def transform(section, Q):
            return (2**(j+1)*p[i]-k[i]*2-Q[i] \
                    for p in section for i in xrange(2))
        Q_00, Q_01 = Point(0, 0), Point(0, 1)
        Q_10, Q_11 = Point(1, 0), Point(1, 1)
        c10, c01, c11 = 0, 0, 0
        for section in self.contour.each():
            KQ00, LQ00 = self.contour.get_KL(transform(section, Q_00))
            KQ01, LQ01 = self.contour.get_KL(transform(section, Q_01))
            KQ10, LQ10 = self.contour.get_KL(transform(section, Q_10))
            KQ11, LQ11 = self.contour.get_KL(transform(section, Q_11))
            c10 += LQ00.x + LQ01.x + KQ10.x \
                 - LQ10.x + KQ11.x - LQ11.x
            c01 += LQ00.y + LQ10.y + KQ01.y \
                 - LQ01.y + KQ11.y - LQ11.y
            c11 += LQ00.x - LQ01.x + KQ10.x \
                 - LQ10.x - KQ11.x + LQ11.x
        return c01, c10, c11

    def g(self, p):
        s = self.area
        E = [Point(0,1), Point(1,0), Point(1,1)]
        for j in xrange(self.max_j+1):
            for kx in xrange(2**j):
                for ky in xrange(2**j):
                    k = Point(kx, ky)
                    cs = self.all_c[(j, kx,ky)]
                    for i, e in enumerate(E):
                        psi = self.psi(p, e, j, k)
                        if psi > 0:     s += cs[i]
                        elif psi < 0:   s -= cs[i]
        return s

    def get(self):
        px_arr = [self.g(p) for p in self.lattice]
        px_mat = [px_arr[i*self.w : (i+1)*self.w] for i in xrange(self.h)]
        return px_mat

    def get_fast(self): # 100x faster than get()
        from util.getpx import get_px as get_cpp
        px_arr = get_cpp(self.area, self.max_j, self.all_c, self.lattice)
        px_mat = [px_arr[i*self.w : (i+1)*self.w] for i in xrange(self.h)]
        return px_mat

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    import cv2, numpy as np,time
    from contour import *

    ts = time.time()
    contour= Line.Contour([(8,8), (60,12), (20,28)])
    raster = Rasterizer(contour, 64, 64).get_fast()
    raster = np.array(np.asarray(raster)*255+0.5, np.uint8)
    cv2.imwrite('var/Line.png', raster)

    contour= QuadraticBezier.Contour([(8,8), (56,8), (56,56), (8,56)])
    raster = Rasterizer(contour, 64, 64).get_fast()
    raster = np.array(np.asarray(raster)*255+0.5, np.uint8)
    cv2.imwrite('var/QuadraticBezier.png', raster)

    contour= CubicBezier.Contour([(8,8),(12,8),(56,24),(56,56),(24,56),(8,24)])
    raster = Rasterizer(contour, 64, 64).get_fast()
    raster = np.array(np.asarray(raster)*255+0.5, np.uint8)
    cv2.imwrite('var/CubicBezier.png', raster)
    print time.time() - ts