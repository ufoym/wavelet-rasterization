import cv2, time, numpy as np
from random import randint, choice
from rasterizer import Rasterizer
from contour import *

if __name__ == '__main__':

    def single_test(tp, n, w, h, z, channel):
        contour = tp.Contour([(randint(1,h-1), randint(1,w-1)) \
            for i in xrange(n)])
        if contour.area() < 0:  return None

        ts = time.time()
        raster = Rasterizer(contour, w, h).get()
        print '%s\ttime: %2.1fs' % (contour, time.time()-ts)

        raster = np.array(np.asarray(raster)*255+0.5, np.uint8)
        raster = cv2.resize(raster, (w*z,h*z), interpolation=cv2.INTER_NEAREST)
        img = np.zeros((h*z,w*z,3), np.uint8)
        img[:,:,channel] = raster
        line_color = [0 if channel == ch else 255 for ch in xrange(3)]
        for sx, sy, ex, ey in contour.to_lines():
            cv2.line(img, (int(sy*z),int(sx*z)), (int(ey*z),int(ex*z)), 
                line_color, 2, cv2.CV_AA)
        return img

    w, h, z = 17, 13, 50
    params = [  (Line, 3, w, h, z, 0), 
                (QuadraticBezier, 4, w, h, z, 1), 
                (CubicBezier, 3, w, h, z, 2)    ]
    while True:
        vis_img = single_test(*choice(params))
        if vis_img is None:    continue
        cv2.namedWindow('raster')
        cv2.imshow('raster', vis_img)
        cv2.waitKey(1)