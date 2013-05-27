Wavelet Rasterization
=====================

Wavelet rasterization is a method for analytically calculating an anti-aliased rasterization of arbitrary polygons or shape bounded by Bezier curves. For more details, please read the following paper:

Manson, Josiah, and Scott Schaefer. **"Wavelet rasterization."** Computer Graphics Forum. Vol. 30. No. 2. Blackwell Publishing Ltd, 2011.

This is a python implementation of the algorithm. Currently it support three types of contours:
* Polygon
* Quadratic Bezier Contour
* Cubic Bezier Contour

Rasterizing these contours are very simple:

1. create a specific `Contour` object (`Line.Contour`, `QuadraticBezier.Contour`, or `CubicBezier.Contour`)
2. use this contour to construct a `Rasterizer` object
3. call the method `get()` of the `Rasterizer` object, you get an array of pixels, each of which has a value range from 0 to 1 that indicates the local transparent of the shape.

An example session could like:

    import cv2, numpy as np             # for image IO
    from contour import *
    from rasterizer import Rasterizer
    
    ## rasterize a polyline
    contour = Line.Contour([(2,2), (15,3), (5,7)])
    raster = Rasterizer(contour, 32, 32).get()
    raster = np.array(np.asarray(raster)*255+0.5, np.uint8)
    cv2.imwrite('var/Line.png', raster)

    ## rasterize a quadratic bezier contour
    contour = QuadraticBezier.Contour([(2,2), (14,2), (14,14), (2,14)])
    raster = Rasterizer(contour, 32, 32).get()
    raster = np.array(np.asarray(raster)*255+0.5, np.uint8)
    cv2.imwrite('var/QuadraticBezier.png', raster)

    ## rasterize a cubic bezier contour
    contour = CubicBezier.Contour([(2,2),(6,2),(14,6),(14,14),(6,14),(2,6)])
    raster = Rasterizer(contour, 32, 32).get()
    raster = np.array(np.asarray(raster)*255+0.5, np.uint8)
    cv2.imwrite('var/CubicBezier.png', raster)