#include <math.h>
#define EQN_EPS 1e-9
#define DLLEXPORT extern "C" __declspec(dllexport)

/********************************************************
*                                                       *
* This function determines if a double is small enough  *
* to be zero. The purpose of the subroutine is to try   *
* to overcome precision problems in math routines.      *
*                                                       *
********************************************************/

static int isZero(double x)
{
    return x > -EQN_EPS && x < EQN_EPS;
}


DLLEXPORT int solveLinear(double c1, double c0, 
                          double & s0)
{
    if (isZero(c1))
        return 0;
    s0 = - c0 / c1;
    return 1;
}



/********************************************************
*                                                       *
* This function determines the roots of a quadric       *
* equation.                                             *
* It takes as parameters a pointer to the three         *
* coefficient of the quadric equation (the c[2] is the  *
* coefficient of x2 and so on) and a pointer to the     *
* two element array in which the roots are to be        *
* placed.                                               *
* It outputs the number of roots found.                 *
*                                                       *
********************************************************/

DLLEXPORT int solveQuadric(double c2, double c1, double c0, 
                           double & s0, double & s1)
{
    double p, q, D;

    // make sure we have a d2 equation

    if (isZero(c2))
        return solveLinear(c1, c0, s0);


    // normal for: x^2 + px + q
    p = c1 / (2.0 * c2);
    q = c0 / c2;
    D = p * p - q;

    if (isZero(D))
    {
        // one double root
        s0 = s1 = -p;
        return 1;
    }

    if (D < 0.0)
        // no real root
            return 0;

    else
    {
        // two real roots
        double sqrt_D = sqrt(D);
        s0 = sqrt_D - p;
        s1 = -sqrt_D - p;
        return 2;
    }
}



/********************************************************
*                                                       *
* This function determines the roots of a cubic         *
* equation.                                             *
* It takes as parameters a pointer to the four          *
* coefficient of the cubic equation (the c[3] is the    *
* coefficient of x3 and so on) and a pointer to the     *
* three element array in which the roots are to be      *
* placed.                                               *
* It outputs the number of roots found                  *
*                                                       *
********************************************************/

DLLEXPORT int solveCubic(double c3, double c2, double c1, double c0, 
                         double & s0, double & s1, double & s2)
{
    int i, num;
    double  sub,
        A, B, C,
        sq_A, p, q,
        cb_p, D;

    if (isZero(c3))
        return solveQuadric(c2, c1, c0, s0, s1);

    // normalize the equation:x ^ 3 + Ax ^ 2 + Bx  + C = 0
    A = c2 / c3;
    B = c1 / c3;
    C = c0 / c3;

    // substitute x = y - A / 3 to eliminate the quadric term: x^3 + px + q = 0

    sq_A = A * A;
    p = 1.0/3.0 * (-1.0/3.0 * sq_A + B);
    q = 1.0/2.0 * (2.0/27.0 * A *sq_A - 1.0/3.0 * A * B + C);

    // use Cardano's formula

    cb_p = p * p * p;
    D = q * q + cb_p;

    if (isZero(D))
    {
        if (isZero(q))
        {
            // one triple solution
            s0 = 0.0;
            num = 1;
        }
        else
        {
            // one single and one double solution
            double u = cbrt(-q);
            s0 = 2.0 * u;
            s1 = - u;
            num = 2;
        }
    }
    else
        if (D < 0.0)
        {
            // casus irreductibilis: three real solutions
            double phi = 1.0/3.0 * acos(-q / sqrt(-cb_p));
            double t = 2.0 * sqrt(-p);
            s0 = t * cos(phi);
            s1 = -t * cos(phi + M_PI / 3.0);
            s2 = -t * cos(phi - M_PI / 3.0);
            num = 3;
        }
        else
        {
            // one real solution
            double sqrt_D = sqrt(D);
            double u = cbrt(sqrt_D + fabs(q));
            if (q > 0.0)
                s0 = - u + p / u ;
            else
                s0 = u - p / u;
            num = 1;
        }

        // resubstitute
        sub = 1.0 / 3.0 * A;
        s0 -= sub;
        s1 -= sub;
        s2 -= sub;
        return num;
}
