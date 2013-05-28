from ctypes import *
import os
os.environ['PATH'] = os.path.dirname(__file__) + ';' + os.environ['PATH']
api = CDLL('solver.dll')

def linear(a, b):
    c1, c0 = c_double(a), c_double(b) 
    x = [c_double()]
    num = api.solveLinear(c1, c0, byref(x[0]))
    return [x[i].value for i in xrange(num)]

def quadric(a, b, c):
    c2, c1, c0 = c_double(a), c_double(b), c_double(c) 
    x = [c_double(), c_double()]
    num = api.solveQuadric(c2, c1, c0, byref(x[0]), byref(x[1]))
    return [x[i].value for i in xrange(num)]

def cubic(a, b, c, d):
    c3, c2, c1, c0 = c_double(a), c_double(b), c_double(c), c_double(d) 
    x = [c_double(), c_double(), c_double()]
    num = api.solveCubic(c3, c2, c1, c0, byref(x[0]), byref(x[1]), byref(x[2]))
    return [x[i].value for i in xrange(num)]

if __name__ == '__main__':
    from random import randint
    import time

    print cubic(2,-4,-22, 24) # [4, -3, 1]
    print cubic(3,-10,14, 27) # [-1]
    print cubic(1, 6, 12,  8) # [-2]

    eps = 1e-9
    ts = time.time()
    for i in xrange(1000000):
        a,b,c,d = randint(0,100),randint(0,100),randint(0,100),randint(0,100)
        for x in linear(a, b):
            if abs(a*x + b) > eps: print 'l',
        for x in quadric(a, b, c):
            if abs(a*x**2 + b*x + c) > eps: print 'q',
        for x in cubic(a, b, c, d):
            if abs(a*x**3 + b*x**2 + c*x + d) > eps: print 'c',
    print time.time() - ts