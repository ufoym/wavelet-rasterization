import math, collections
from util.BezierClipping import Bezier, clip

# -----------------------------------------------------------------------------
Point = collections.namedtuple('Point', 'x y')
def area(polybez):
    def det(a, b):  return a[0] * b[1] - a[1] * b[0]
    s = 0
    for i in xrange(0, len(polybez), 3):
        v3 = polybez[i]
        v2 = polybez[i-1]
        v1 = polybez[i-2]
        v0 = polybez[i-3]
        s += 3./10 * det(v0, v1) + 3./20 * det(v1, v2) + 3./10 * det(v2, v3) \
           + 3./20 * det(v0, v2) + 3./20 * det(v1, v3) + 1./20 * det(v0, v3)
    return s

class Rasterizer:

    def __init__(self, polybez, w, h):
        self._w = w
        self._h = h
        self._max_j = int(math.ceil(math.log(max(w,h), 2)))-1
        self._wh = 2 ** (self._max_j+1)
        def normalize(p):
            return (p[0]/float(self._wh), p[1]/float(self._wh))
        self._polybez = [normalize(p) for p in polybez]
        self._area = area(self._polybez)
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
        def get_KL(_start, b0, b1, _end, Q):
            def transform(p):
                return Point(2**(j+1)*p[0]-k.x*2-Q.x, 
                             2**(j+1)*p[1]-k.y*2-Q.y)
            sx, sy = transform(_start)
            ex, ey = transform(_end)
            b0x, b0y = transform(b0)
            b1x, b1y = transform(b1)
            bezier = Bezier(sx,sy, b0x,b0y, b1x,b1y, ex,ey)
            Kx, Ky, Lx, Ly = 0, 0, 0, 0
            for sec in clip(0, 1, 1, 0, bezier):
                v3 = Point(sec.x0, sec.y0)
                v2 = Point(sec.x1, sec.y1)
                v1 = Point(sec.x2, sec.y2)
                v0 = Point(sec.x3, sec.y3)
                if v0.x == 1 and v1.x == 1 and v2.x == 1 and v3.x == 1\
                or v0.y == 1 and v1.y == 1 and v2.y == 1 and v3.y == 1:
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
        Q_00, Q_01 = Point(0, 0), Point(0, 1)
        Q_10, Q_11 = Point(1, 0), Point(1, 1)
        c10, c01, c11 = 0, 0, 0
        for i in xrange(0, len(self._polybez), 3):
            v1 = self._polybez[i]
            b1 = self._polybez[i-1]
            b0 = self._polybez[i-2]
            v0 = self._polybez[i-3]
            KQ00, LQ00 = get_KL(v0, b0, b1, v1, Q_00)
            KQ01, LQ01 = get_KL(v0, b0, b1, v1, Q_01)
            KQ10, LQ10 = get_KL(v0, b0, b1, v1, Q_10)
            KQ11, LQ11 = get_KL(v0, b0, b1, v1, Q_11)
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
    from util.BezierClipping import evaluate
    w, h, z = 13, 17, 30

    while True:
        polybez = [(randint(1,h-1), randint(1,w-1)) for i in xrange(3)]
        if area(polybez) < 0:
            continue
        polybez = [(0,0),(h,0),(0,w)]

        ts = time.time()
        raster = Rasterizer(polybez, w, h).get()
        print '%s\ttime: %2.1fs' % (polybez, time.time()-ts)

        raster = np.array(np.asarray(raster)*255+0.5, np.uint8)
        raster = cv2.cvtColor(raster, cv2.COLOR_GRAY2BGR)
        raster = cv2.resize(raster, (w*z,h*z), interpolation=cv2.INTER_NEAREST)

        tts = np.linspace(0,1, num=50)
        for i in xrange(0, len(polybez), 3):
            v3 = polybez[i]
            v2 = polybez[i-1]
            v1 = polybez[i-2]
            v0 = polybez[i-3]
            bezier = Bezier(v0[0],v0[1], v1[0],v1[1], v2[0],v2[1], v3[0],v3[1])
            for start, end in zip(tts[:-1], tts[1:]):
                sx, sy = evaluate(bezier, start)
                ex, ey = evaluate(bezier, end)
                cv2.line(raster, (int(sy*z),int(sx*z)), (int(ey*z),int(ex*z)), 
                    (0,0,255), 2, cv2.CV_AA)
        cv2.namedWindow('raster')
        cv2.imshow('raster', raster)

        cv2.waitKey(0)