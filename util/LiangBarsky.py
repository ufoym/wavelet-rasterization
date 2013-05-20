from collections import namedtuple
Point = namedtuple('Point', 'x y')

def LiangBarsky(left, right, bottom, top, start, end):
    """
    Parameters
    ----------
    left, right, bottom, top : numeric
        Input numerics that defines the x/y clipping values for the border.
    start, end : Point
        Input two points that defines the line segment to be clipped.

    Returns
    -------
    clipped_start, clipped_end : Point
        Output two points that defines the clipped line segment.
        Raise Exception if no clipped line exists.

    References
    ----------
        http://www.skytopia.com/project/articles/compsci/clipping.html

    """
    t0, t1 = 0, 1
    xdelta = end.x-start.x
    ydelta = end.y-start.y
    for edge in xrange(4):  # traverse through left, right, bottom, top edges.
        if   edge == 0:   p, q = -xdelta, -(left-start.x) 
        elif edge == 1:   p, q =  xdelta,  (right-start.x)
        elif edge == 2:   p, q =  ydelta,  (bottom-start.y)
        elif edge == 3:   p, q = -ydelta, -(top-start.y)
        r = q / float(p)
        if p == 0 and q < 0:    raise Exception('no line')
        if p < 0:
            if r > t1:          raise Exception('no line')
            elif r > t0:        t0 = r   # line is clipped!
        elif p > 0:
            if r < t0:          raise Exception('no line')
            elif r < t1:        t1 = r   # line is clipped!
    clipped_start = Point(start.x + t0*xdelta, start.y + t0*ydelta)
    clipped_end = Point(start.x + t1*xdelta, start.y + t1*ydelta)
    return clipped_start, clipped_end

if __name__ == '__main__':
    import cv2, numpy as np
    from random import randint
    w, h = 360, 360
    left, right, bottom, top = int(w*0.25),int(w*0.75),int(h*0.75),int(h*0.25)
    while True:
        start = Point(randint(1,w-1), randint(1,h-1))
        end = Point(randint(1,w-1), randint(1,h-1))

        vis = np.zeros((h,w,3), np.uint8)
        cv2.rectangle(vis, (left, top), (right, bottom), (255,0,0))
        cv2.line(vis, start, end, (255,255,0), lineType = cv2.CV_AA)
        cv2.putText(vis,'s',start, cv2.FONT_HERSHEY_PLAIN, 1.0, (0,255,0))
        cv2.putText(vis,'e',end, cv2.FONT_HERSHEY_PLAIN, 1.0, (0,255,0))

        try:
            cstart, cend = LiangBarsky(left, right, bottom, top, start, end)
            vstart = (int(cstart.x), int(cstart.y))
            vend = (int(cend.x), int(cend.y))
            cv2.putText(vis,'1',vstart, cv2.FONT_HERSHEY_PLAIN, 1.0,(0,0,255))
            cv2.putText(vis,'0',vend, cv2.FONT_HERSHEY_PLAIN, 1.0,(0,0,255))
            cv2.line(vis, vstart, vend, (0,255,255), lineType = cv2.CV_AA)
        except:
            pass
        cv2.namedWindow('result')
        cv2.imshow('result', vis)
        cv2.waitKey(100)