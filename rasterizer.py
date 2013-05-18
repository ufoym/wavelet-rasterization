import Polygon, cv2, numpy as np, time
from collections import namedtuple

w, h, z = 8, 8, 30
poly = ((1,2), (w-1,1), (w-1,h-1))

Point = namedtuple('Point', 'x y')
lattice = [Point(x,y) for x in xrange(w) for y in xrange(h)]

def normalize_point(p):
	return (p[0]/float(w), p[1]/float(h))
normed_poly = Polygon.Polygon([normalize_point(p) for p in poly])
normed_box = [Point(*normalize_point(p)) for p in [(0,0), (1,0), (1,1), (0,1)]]
normed_lattice = [Point(*normalize_point(p)) for p in lattice]

# ------------------------------------------------------------------------------

def psi_2d(p, e, j, k):
	def psi_1d(p, e):
		if e == 0:	return 1 if 0 <= p < 1 else 0
		else:		return (1 if 0 <= p < 0.5 else -1) if 0 <= p < 1 else 0 
	return 2**j * psi_1d(2**j*p.x-k.x, e.x) * psi_1d(2**j*p.y-k.y, e.y)

def c(e, j, k):
	s = 0
	for p in normed_lattice:
		psi_p = psi_2d(p, e, j, k)
		if psi_p != 0:
			box = Polygon.Polygon([(b.x+p.x, b.y+p.y) for b in normed_box])
			px_poly = box & normed_poly
			s += px_poly.area() * psi_p
	return s

def g(p, d):
	s = 0
	# --------------------------------------------------------------------------
	j = 0
	for kx in xrange(2**j):
		for ky in xrange(2**j):
			k = Point(kx, ky)
			s += c(Point(0,0), j, k) * psi_2d(p, Point(0,0), j, k)
	# --------------------------------------------------------------------------
	for j in xrange(d):
		for kx in xrange(2**j):
			for ky in xrange(2**j):
				k = Point(kx, ky)
				for e in [Point(0,1), Point(1,0), Point(1,1)]:
					s += c(e, j, k) * psi_2d(p, e, j, k)
	return s

# ------------------------------------------------------------------------------
ts = time.time()
raster = [g(p, 4) for p in normed_lattice]
print 'time: %2.1fs' % (time.time()-ts)

raster = np.array(np.asarray(raster)*30, np.uint8)
raster = np.reshape(raster, (w,h))
raster = cv2.cvtColor(raster, cv2.COLOR_GRAY2BGR)
raster = cv2.resize(raster, (w*z,h*z), interpolation = cv2.INTER_NEAREST)

poly_vis = np.array([(p[1]*z, p[0]*z) for p in poly])
cv2.polylines(raster, [poly_vis],True,(255,255,0),lineType=cv2.CV_AA)
cv2.namedWindow('raster')
cv2.imshow('raster', raster)

cv2.waitKey(0)