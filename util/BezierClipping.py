import numpy as np
from collections import namedtuple

# -----------------------------------------------------------------------------
Bezier = namedtuple('Bezier', 'x0 y0 x1 y1 x2 y2 x3 y3')

def evaluate(bez, t):
    return (bez.x0*(1-t)**3 + 3*bez.x1*(1-t)**2*t \
            + 3*bez.x2*(1-t)*t**2 + bez.x3*t**3,
            bez.y0*(1-t)**3 + 3*bez.y1*(1-t)**2*t \
            + 3*bez.y2*(1-t)*t**2 + bez.y3*t**3)

def find_intersections(left, right, bottom, top, bez):
    def is_t_in(t, eps = 1e-5):
        pt = evaluate(bez, t)
        return left-eps <=pt[0]<= right+eps and top-eps <=pt[1]<= bottom+eps

    def find_cubic_real_root(a, b, c, d, eps = 1e-6):
        '''solve ax^3 + bx^2 + cx + d = 0'''
        results = []
        for r in np.roots([a, b, c, d]):
            if type(r) == np.complex128:
                if abs(r.imag) < eps:   results.append(r.real)
            else:                       results.append(r)
        return results

    ax, bx = -bez.x0+3*bez.x1-3*bez.x2+bez.x3, 3*bez.x0-6*bez.x1+3*bez.x2
    cx, _dx = 3*bez.x1-3*bez.x0, bez.x0
    ay, by = -bez.y0+3*bez.y1-3*bez.y2+bez.y3, 3*bez.y0-6*bez.y1+3*bez.y2
    cy, _dy = 3*bez.y1-3*bez.y0, bez.y0
    ts = [0]
    ts += find_cubic_real_root(ax, bx, cx, _dx-left)
    ts += find_cubic_real_root(ax, bx, cx, _dx-right)
    ts += find_cubic_real_root(ay, by, cy, _dy-bottom)
    ts += find_cubic_real_root(ay, by, cy, _dy-top)
    ts.append(1)
    intersections = [t for t in ts if 0 <= t <= 1 and is_t_in(t)]
    return sorted(intersections)

def subsection(bez, t0, t1):
    u0 = 1.0 - t0
    u1 = 1.0 - t1

    qxa = bez.x0*u0*u0 + bez.x1*2*t0*u0 + bez.x2*t0*t0
    qxb = bez.x0*u1*u1 + bez.x1*2*t1*u1 + bez.x2*t1*t1
    qxc = bez.x1*u0*u0 + bez.x2*2*t0*u0 + bez.x3*t0*t0
    qxd = bez.x1*u1*u1 + bez.x2*2*t1*u1 + bez.x3*t1*t1

    qya = bez.y0*u0*u0 + bez.y1*2*t0*u0 + bez.y2*t0*t0
    qyb = bez.y0*u1*u1 + bez.y1*2*t1*u1 + bez.y2*t1*t1
    qyc = bez.y1*u0*u0 + bez.y2*2*t0*u0 + bez.y3*t0*t0
    qyd = bez.y1*u1*u1 + bez.y2*2*t1*u1 + bez.y3*t1*t1

    sec = Bezier(   qxa*u0 + qxc*t0, qya*u0 + qyc*t0,
                    qxa*u1 + qxc*t1, qya*u1 + qyc*t1,
                    qxb*u0 + qxd*t0, qyb*u0 + qyd*t0,
                    qxb*u1 + qxd*t1, qyb*u1 + qyd*t1 )
    return sec

def clip(left, right, bottom, top, bez):
    intersections = find_intersections(left, right, bottom, top, bez)
    if len(intersections) == 0:
        return None
    else:
        t0, t1 = intersections[0], intersections[-1]
        return subsection(bez, t0, t1)

if __name__ == '__main__':
    import cv2, numpy as np
    from random import randint
    w, h = 500, 500
    rect = (int(w*0.25), int(h*0.25), int(w*0.75), int(h*0.75))

    while True:
        vis = np.zeros((h,w,4), np.uint8)
        bezier = Bezier(randint(0,w-1),randint(0,h-1), 
                        randint(0,w-1),randint(0,h-1), 
                        randint(0,w-1),randint(0,h-1), 
                        randint(0,w-1),randint(0,h-1))
        print bezier
        ts = np.linspace(0,1, num=50)
        for start, end in zip(ts[:-1], ts[1:]):
            sx, sy = evaluate(bezier, start)
            ex, ey = evaluate(bezier, end)
            cv2.line(vis, (int(sx), int(sy)), (int(ex), int(ey)), 
                (0,0,255), 2, cv2.CV_AA)
        cv2.rectangle(vis, rect[:2], rect[2:], (255,0,0))

        sec = clip(rect[0],rect[2],rect[3],rect[1], bezier)
        if sec:
            for start, end in zip(ts[:-1], ts[1:]):
                sx, sy = evaluate(sec, start)
                ex, ey = evaluate(sec, end)
                cv2.line(vis, (int(sx), int(sy)), (int(ex), int(ey)), 
                    (0,255,255), 3, cv2.CV_AA)
            cv2.circle(vis, bezier[:2], 2, (0,255,255), 5)
        cv2.circle(vis, bezier[:2], 2, (0,0,255), 5)

        cv2.namedWindow('raster')
        cv2.imshow('raster', vis[:,:,:3])
        cv2.waitKey(1)