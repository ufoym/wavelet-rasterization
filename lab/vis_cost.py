import cv2, numpy as np, cairo
from collections import namedtuple

# -----------------------------------------------------------------------------
Color = namedtuple('Color', 'r g b')
Bezier = namedtuple('Bezier', 'x0 y0 x1 y1 x2 y2 x3 y3')

def rasterize(img, bezier, color, w, h, alpha = 1.0):
    surface = cairo.ImageSurface.create_for_data(img, cairo.FORMAT_ARGB32,w,h)
    cr = cairo.Context(surface)
    cr.move_to(*bezier[:2])
    cr.curve_to(*bezier[2:])
    cr.close_path()
    cr.set_source_rgba(color.r/255., color.g/255., color.b/255., alpha)
    cr.fill()

if __name__ == '__main__':
    w, h = 800, 800
    color = Color(255, 0, 0)

    # vis = np.zeros((h,w,3), np.uint8)
    # bezier = Bezier(w*0.25,h*0.25, w*0.5,h*0.25, w*0.75,h*0.5, w*0.75,h*0.75)
    # for i in xrange(4):
    #     cv2.circle(vis, (int(bezier[i*2]), int(bezier[i*2+1])), 3, (255,255,0))
    # tmp = np.zeros((h,w,4), np.uint8)
    # tmp[:,:,:3] = vis
    # vis = tmp


    org_img = np.zeros((h,w,4), np.uint8)
    bezier = Bezier(w*0.25,h*0.25, w*0.5,h*0.25, w*0.75,h*0.5, w*0.75,h*0.75)
    rasterize(org_img, bezier, color, w, h)

    sample_num = 200
    cost_img = np.zeros((sample_num,sample_num), np.float)
    X, Y = [], []
    for j, x in enumerate(np.linspace(0, w, num=sample_num)):
        for i, y in enumerate(np.linspace(0, h, num=sample_num)):
            cur_img = np.zeros((h,w,4), np.uint8)
            bezier = Bezier(w*0.25,h*0.25, x,y, w*0.75,h*0.5, w*0.75,h*0.75)
            rasterize(cur_img, bezier, color, w, h)
            cost = np.linalg.norm( cv2.absdiff(cur_img, org_img) )
            cost_img[i,j] = cost
            X.append(x)
            Y.append(y)
    print cost_img

    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, cost_img, cmap='autumn', cstride=2, rstride=2)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("cost")
    plt.show()

    # vis = vis[:,:,:3]
    # cv2.namedWindow('result')
    # cv2.imshow('result', vis)
    # cv2.waitKey(0)