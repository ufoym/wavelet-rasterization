Wavelet Rasterization
=====================

Wavelet rasterization is a method for analytically calculating an anti-aliased rasterization of arbitrary polygons or shape bounded by Bezier curves. For more details, please read the following paper:

Manson, Josiah, and Scott Schaefer. **"Wavelet rasterization."** Computer Graphics Forum. Vol. 30. No. 2. Blackwell Publishing Ltd, 2011.

This is a python implementation of the algorithm. Currently it supports three types of contours:
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
    contour = Line.Contour([(4,4), (30,6), (10,14)])
    raster = Rasterizer(contour, 32, 32).get()
    raster = np.array(np.asarray(raster)*255+0.5, np.uint8)
    cv2.imwrite('var/Line.png', raster)

Assuming everything is working OK, the examples should generate the following image:

![line](https://f.cloud.github.com/assets/2270240/566814/829ae7dc-c6a0-11e2-91ee-45184f5a8a1d.png)

You can also rasterize a quadratic bezier contour:

    contour = QuadraticBezier.Contour([(4,4), (28,4), (28,28), (4,28)])
    raster = Rasterizer(contour, 32, 32).get()
    raster = np.array(np.asarray(raster)*255+0.5, np.uint8)
    cv2.imwrite('var/QuadraticBezier.png', raster)

![quadraticbezier](https://f.cloud.github.com/assets/2270240/566816/8646ba6e-c6a0-11e2-9cd0-19cd058768b8.png)

or a cubic bezier contour:

    ## rasterize a cubic bezier contour
    contour = CubicBezier.Contour([(4,4),(12,4),(28,12),(28,28),(12,28),(4,12)])
    raster = Rasterizer(contour, 32, 32).get()
    raster = np.array(np.asarray(raster)*255+0.5, np.uint8)
    cv2.imwrite('var/CubicBezier.png', raster)

![cubicbezier](https://f.cloud.github.com/assets/2270240/566823/a0f4b2d0-c6a0-11e2-8b89-e1045d573714.png)
