#include "octree.h"
#include "insert_bez2.h"
#include "global.h"
#include "boundedarray.h"
#include "array.h"
#include <algorithm>

using namespace std;

//Array<Bez2> svg_curves;
//
//void saveSVG()
//{
//	double scale = 1024;
//
//	char fn[1024];
//	sprintf(fn, "output.svg");
//	FILE *f = fopen(fn, "wb");
//
//	fprintf(f, "<?xml version=\"1.0\" standalone=\"no\"?> <!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
//	fprintf(f, "<svg width=\"1024\" height=\"1024\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n");
//
//	// draw contours
//	for (int i = 0; i < svg_curves.s; i++)
//	{
//		fprintf(f, "<path d=\"");
//		fprintf(f, "M %f %f Q %f %f %f %f\" style=\"fill:none;stroke:black;stroke-width:1\"/>\n", 
//			svg_curves[i][0][0] * scale, scale - svg_curves[i][0][1] * scale, 
//			svg_curves[i][1][0] * scale, scale - svg_curves[i][1][1] * scale, 
//			svg_curves[i][2][0] * scale, scale - svg_curves[i][2][1] * scale);
//	}
//
//	fprintf(f, "</svg>\n");
//	fclose(f);
//}


void cut_bez2(Bez2 &c, Bez2 &a, Bez2 &b, float t)
{
	vect2f a1, b1, a2;
	float t1 = 1-t;

	a[0] = c[0];
	b[2] = c[2];

	a1 = c[0]*t1 + c[1]*t;
	b1 = c[1]*t1 + c[2]*t;
	a[1] = a1;
	b[1] = b1;

	a2 = a1*t1 + b1*t;
	a[2] = a2;
	b[0] = a[2];
}


void find_roots_bez2_analytic(Bez2 &curve, BoundedArray<float, 2> &tvals, int dir)
{
	// subtract 1 from control points, because solving for == 1, but using equation for == 0.
	float a = curve[0][dir] - 1;
	float b = curve[1][dir] - 1;
	float c = curve[2][dir] - 1;

	double disc = b*b - a*c; // discriminant

	if (disc < 0)
	{
		// discriminant < 0 -> no intersection -> all on same side.
		// curve passes through end points, so test one end point for which side.
		// this will be handled later on
	}
	else
	{

		double denom = a-2*b+c; // denominator

		if (fabs(denom) < 1e-9)
		{
			// the equation is nearly linear, so solve a linear equation
			double v = c-a;
			
			if (fabs(v) < 1e-9)
			{
				// even the linear case is degenerate, so just choose a side that contains an end point and hope for the best.
				// i.e. add no t value
			}
			else
			{
				float t = -a/v;
				if (t > 0 && t < 1)
					tvals.push_back(t);
			}
		}
		else
		{
			//double denom_inv = 1.0 / denom;
			double root = sqrt(disc);
			double num = a - b;

			// add t values
			float t;
			t = (num-root)/denom;
			if (t > 0 && t < 1)
				tvals.push_back(t);
			t = (num+root)/denom;
			if (t > 0 && t < 1)
				tvals.push_back(t);

			// ensure sorted order
			if (tvals.s > 1 && tvals[0] > tvals[1])
				swap(tvals[0], tvals[1]);
		}
	}
}

void cut_bez2_half(vect3f &c, vect3f &a, vect3f &b)
{
	a[0] = c[0];
	b[2] = c[2];

	a[1] = (c[0] + c[1])*.5;
	b[1] = (c[1] + c[2])*.5;

	a[2] = (a[1] + b[1])*.5;
	b[0] = a[2];
}

bool check_bez2_sign(vect3f &curve)
{
	if ((curve[0] < 1) != (curve[1] < 1))
		return true;
	else if ((curve[1] < 1) != (curve[2] < 1))
		return true;
	else
		return false;
}

void find_roots_bez2_bisection(vect3f &curve, BoundedArray<float, 2> &tvals, int depth, float t0, float t1)
{
	vect3f c1, c2;
	cut_bez2_half(curve, c1, c2);

	float t = (t0+t1)*.5;

	if (depth < 20)
	{
		if (check_bez2_sign(c1))
			find_roots_bez2_bisection(c1, tvals, depth+1, t0, t);
		if (check_bez2_sign(c2))
			find_roots_bez2_bisection(c2, tvals, depth+1, t, t1);
	}
	else if (check_bez2_sign(curve))
		tvals.push_back(t);
}

void find_roots_bez2_bisection(Bez2 &c, BoundedArray<float, 2> &tvals, int dir)
{
	vect3f curve;
	curve.set(c[0][dir], c[1][dir], c[2][dir]);
	find_roots_bez2_bisection(curve, tvals, 0, 0, 1);
}

float eval_bez2(vect3f &c, double t)
{
	float t1 = 1-t;
	float a1 = c[0]*t1 + c[1]*t;
	float b1 = c[1]*t1 + c[2]*t;
	return a1*t1 + b1*t;
}

#if 0
void find_roots_bez2_bisection_line(vect3f &curve, BoundedArray<float, 2> &tvals, float t0, float t1)
{
	if ((curve[0] < 1) != (curve[2] < 1))
	{
		// has one root
		if (fabs(curve[0] - 2*curve[1] + curve[2]) < 1e-6)
		{
			// linear with one root, so find linear intersection
			float t = t0 + (1 - curve[0]) * (t1 - t0) / (curve[2] - curve[0]);
			tvals.push_back(t);
		}
		else
		{
			// subdivide on the single root
			vect3f c1, c2;
			cut_bez2_half(curve, c1, c2);

			float t = (t0+t1)*.5;

			if (check_bez2_sign(c1))
				find_roots_bez2_bisection_line(c1, tvals, t0, t);
			else
				find_roots_bez2_bisection_line(c2, tvals, t, t1);
		}
	}
	else if ((curve[0] < 1) != (curve[1] < 1))
	{
		// up to two roots, so subdivide
		vect3f c1, c2;
		cut_bez2_half(curve, c1, c2);

		float t = (t0+t1)*.5;

		if (check_bez2_sign(c1))
			find_roots_bez2_bisection_line(c1, tvals, t0, t);
		if (check_bez2_sign(c2))
			find_roots_bez2_bisection_line(c2, tvals, t, t1);
	}
}
#else
void find_roots_bez2_bisection_line(vect3f &curve, BoundedArray<float, 2> &tvals, float t0, float t1)
{
	if ((curve[0] < 1) != (curve[2] < 1))
	{
		// has one root
		float t;
		float v0 = curve[0];
		float v1 = curve[1];

		for (int i = 0; ; i++)
		{
			// subdivide on the single root
			t = (t0+t1)*.5;
			float v = eval_bez2(curve, t);

			if (fabs(v0 - 2*v + v1) < 1e-8 || i > 20)
			{
				if ((v0 < 0) != (v < 0))
				{
					v1 = v;
					t1 = t;
				}
				else
				{
					v0 = v;
					t0 = t;
				}

				// linear with one root, so find linear intersection
				t = t0 + (1 - v0) * (t1 - t0) / (v1 - v0);
				tvals.push_back(t);
				break;
			}

			if ((v0 < 1) != (v < 1))
			{
				v1 = v;
				t1 = t;
			}
			else
			{
				v0 = v;
				t0 = t;
			}
		}
	}
	else if ((curve[0] < 1) != (curve[1] < 1))
	{
		// up to two roots, so subdivide
		vect3f c1, c2;
		cut_bez2_half(curve, c1, c2);

		float t = (t0+t1)*.5;

		if (check_bez2_sign(c1))
			find_roots_bez2_bisection_line(c1, tvals, t0, t);
		if (check_bez2_sign(c2))
			find_roots_bez2_bisection_line(c2, tvals, t, t1);
	}
}
#endif

void find_roots_bez2_bisection_line(Bez2 &c, BoundedArray<float, 2> &tvals, int dir)
{
	vect3f curve;
	curve.set(c[0][dir], c[1][dir], c[2][dir]);
	find_roots_bez2_bisection_line(curve, tvals, 0, 1);
}

void split_bez2(Bez2 &curve, BoundedArray<Bez2, 2> &segs1, BoundedArray<Bez2, 2> &segs2, int dir)
{
	segs1.clear();
	segs2.clear();

	//---------------------------------------------------
	// Calculate t values
	//---------------------------------------------------

	BoundedArray<float, 2> tvals;
	tvals.init();

	//find_roots_bez2_analytic(curve, tvals, dir);
	find_roots_bez2_bisection(curve, tvals, dir);
	//find_roots_bez2_bisection_line(curve, tvals, dir);

	//---------------------------------------------------
	// Split based on t values
	//---------------------------------------------------

	if (tvals.s == 0)
	{
		if (curve[0][dir] < 1)
			segs1.push_back(curve);
		else
			segs2.push_back(curve);
	}
	else if (tvals.s == 1)
	{
		segs1.s = segs2.s = 1;
		if (curve[0][dir] < 1)
			cut_bez2(curve, segs1[0], segs2[0], tvals[0]);
		else
			cut_bez2(curve, segs2[0], segs1[0], tvals[0]);
	}
	else
	{
		if (curve[0][dir] < 1)
		{
			segs1.s = 2;
			segs2.s = 1;
			Bez2 tmp;
			cut_bez2(curve, tmp, segs1[1], tvals[1]);
			cut_bez2(tmp, segs1[0], segs2[0], tvals[0]/tvals[1]);
		}
		else
		{
			segs1.s = 1;
			segs2.s = 2;
			Bez2 tmp;
			cut_bez2(curve, tmp, segs2[1], tvals[1]);
			cut_bez2(tmp, segs2[0], segs1[0], tvals[0]/tvals[1]);
		}
	}

	//---------------------------------------------------
	// Subtract 1 from second set of curve segments
	//---------------------------------------------------

	for (int i = 0; i < segs2.s; i++)
	{
		segs2[i][0][dir]--;
		segs2[i][1][dir]--;
		segs2[i][2][dir]--;
	}
}

void calc_coeffs_bez2(float *coeffs, float* v, int i, int j)
{
	// simplify the mapping from Mathematica to C by using common names
	const float p00 = v[0];
	const float p01 = v[1];
	const float p10 = v[2];
	const float p11 = v[3];
	const float p20 = v[4];
	const float p21 = v[5];

	// calc the coeffs
	float c10;
	if (i == 0)
		c10 = (3*p01- 2*p11)*p00 + 2*p01*p10  + (p01 + 2*p11)*p20 - (p00 + 2*p10 + 3*p20)*p21;
	else
		c10 = 2*p11*(p00 - p20) - p01*(-6 + 3*p00 + 2*p10 + p20) + (-6 + p00 + 2*p10 + 3*p20)*p21;

	float c01;
	if (j == 0)
		c01 = 2*p11*p20 + p01*(2*p10 + p20) +(- 2*p10 + 3*p20)*p21 - p00*(3*p01 + 2*p11 + p21);
	else
		c01 = -(p01*(2*p10 + p20)) + p20*(6 - 2*p11 - 3*p21) + 2*p10*p21 + p00*(-6 + 3*p01 + 2*p11 + p21);

	float c11;
	if (i == 0)
	{
		if (j == 0)
			c11 = 3*p00*p01 + 2*p01*p10 - 2*p00*p11 + (p01 + 2*p11)*p20 - (p00 + 2*p10 + 3*p20)*p21;
		else
			c11 = -2*p11*p20 - p01*(2*p10 + p20) + (2*p10 + 3*p20)*p21 + p00*(-3*p01 + 2*p11 + p21);
	}
	else
	{
		if (j == 0)
			c11 = 2*p11*(p00 - p20) - p01*(-6 + 3*p00 + 2*p10 + p20) + (-6 + p00 + 2*p10 + 3*p20)*p21;
		else
			c11 = 2*p11*(-p00 + p20) + p01*(-6 + 3*p00 + 2*p10 + p20) - (-6 + p00 + 2*p10 + 3*p20)*p21;
	}

	// store the coeffs (*-1/(4*6))
	coeffs[0] += c10 * -.04166666666666;
	coeffs[1] += c01 * -.04166666666666;
	coeffs[2] += c11 * -.04166666666666;
}

void insert_bez2(Node *n, Bez2 &p, int depth);
void insert_bez2_leaf(NodeData *n, Bez2 &p);
void insert_bez2_split(Node *n, Bez2 &p, int depth);
void insert_bez2_split_leaf(NodeData *n, Bez2 &p);

__forceinline void insert_bez2_split_leaf(NodeData *n, Bez2 &p)
{
	// scale by a factor of 2
	p[0] *= 2;
	p[1] *= 2;
	p[2] *= 2;

	// split line
	BoundedArray<Bez2, 2> s[4];

	// do full split
	BoundedArray<Bez2, 2> x[2];

	split_bez2(p, x[0], x[1], 0);
	for (int i = 0; i < 2; i++)
	{
		for (int c_i = 0; c_i < x[i].s; c_i++)
		{
			split_bez2(x[i][c_i], s[i], s[i+2], 1);

			for (int j = 0; j < 2; j++)
			{
				int ij = i + j*2;
				for (int c_j = 0; c_j < s[ij].s; c_j++)
				{
					// add to coefficients
					vect2f &v0 = s[ij][c_j].v[0];
					calc_coeffs_bez2(n->coeffs, v0.v, i, j);
					n->is_boundary[ij] = true;
				}
			}
		}
	}
}

void insert_bez2_split(Node *n, Bez2 &p, int depth)
{
	// scale by a factor of 2
	p[0] *= 2;
	p[1] *= 2;
	p[2] *= 2;

	// split line
	BoundedArray<Bez2, 2> s[4];

	// do full split
	BoundedArray<Bez2, 2> x[2];

	split_bez2(p, x[0], x[1], 0);
	for (int i = 0; i < 2; i++)
	{
		for (int c_i = 0; c_i < x[i].s; c_i++)
		{
			split_bez2(x[i][c_i], s[i], s[i+2], 1);

			for (int j = 0; j < 2; j++)
			{
				int ij = i + j*2;
				for (int c_j = 0; c_j < s[ij].s; c_j++)
				{
					// add to coefficients
					vect2f &v0 = s[ij][c_j].v[0];
					calc_coeffs_bez2(n->data.coeffs, v0.v, i, j);
					n->data.is_boundary[ij] = true;

					// recur
					if (depth < g.max_depth - 1)
					{
						// get a valid child pointer (possibly allocating a new node)
						Node *child;
						if (n->children[ij].invalid())
						{
							n->children[ij] = g.node_allocator.alloc();
							child = g.node_allocator.get_ref(n->children[ij]);
							child->init();
						}
						else
						{
							child = g.node_allocator.get_ref(n->children[ij]);
						}

						insert_bez2_split(child, s[ij][c_j], depth + 1);
					}
					else if (depth == g.max_depth - 1)
					{
						// get a valid child pointer (possibly allocating a new node)
						NodeData *child;
						if (n->children[ij].invalid())
						{
							n->children[ij] = g.leaf_allocator.alloc();
							child = g.leaf_allocator.get_ref(n->children[ij]);
							child->init();
						}
						else
						{
							child = g.leaf_allocator.get_ref(n->children[ij]);
						}

						insert_bez2_split_leaf(child, s[ij][c_j]);
					}
				}
			}
		}
	}
}

void insert_bez2_root(Bez2 &p)
{
	//svg_curves.clear();

	// calculate root coefficient
	g.coeff_const += -((2*p[1][1]*p[2][0] + p[0][1]*(2*p[1][0] + p[2][0]) - 2*p[1][0]*p[2][1] - p[0][0]*(2*p[1][1] + p[2][1]))/3.) * .5;
	
	// add line to tree
	insert_bez2_split(g.root, p, 0);

	//saveSVG();
}