#include "octree.h"
#include "insert_line.h"
#include "global.h"

void split_line(Line &p, Line &a, Line &b, bool &ina, bool &inb, int dir)
{
	vect2f *vp = &p[0];
	const float xp = vp->v[dir];
	const char inp = xp < 1;

	vect2f &v = p[1];
	const float &x = v.v[dir];
	const char in = x < 1;

	int sel = inp*2 + in;

	switch (sel)
	{
	case 3:
		//if (inp && in)
		{
			ina = true;
			inb = false;
			a = p;
		}
		break;
	case 2:
		//else if (inp && !in)
		{
			const float div = x - xp;
			//assert(div > 0); // since the sign is different there must be at least an infinitessimal step difference

			float t = (1 - xp)/div;
			// check slows performance by about 5%
			if (t < 0 || t > 1) // only numerical instability can make t be out of range, so assume half because the line is nearly tangent
				t = .5;

			vect2f m;
			int axis1 = 1 - dir;
			float *pp1 = vp->v + axis1;
			m[axis1] = (v[axis1] - *pp1)*t + *pp1;
			m[dir] = 1; // force to be exactly on boundary

			ina = inb = true;
			a[0] = p[0];
			a[1] = m;

			b[0] = m;
			b[1] = p[1];
		}
		break;
	case 1:
		//else if (!inp && in)
		{
			const float div = x - xp;
			//assert(div < 0); // since the sign is different there must be at least an infinitessimal step difference

			float t = (1 - xp)/div;
			// check slows performance by about 5%
			if (t < 0 || t > 1) // only numerical instability can make t be out of range, so assume half because the line is nearly tangent
				t = .5;

			vect2f m;
			int axis1 = 1 - dir;
			float *pp1 = vp->v + axis1;
			m[axis1] = (v[axis1] - *pp1)*t + *pp1;
			m[dir] = 1; // force to be exactly on boundary

			ina = inb = true;
			b[0] = p[0];
			b[1] = m;

			a[0] = m;
			a[1] = p[1];
		}
		break;
	case 0:
		//else
		{
			ina = false;
			inb = true;
			b = p;
		}
	}

	b[0].v[dir]--;
	b[1].v[dir]--;
}

//#define USE_SSE
#ifdef USE_SSE
// I'm so sad that this is slower than the non parallel code
float __declspec(align(16)) x_neg[4] = {1, -1, 1, 0};
float __declspec(align(16)) x_c[16] = 
{
	0.000000, 0.000000, 0.000000, 0,
	0.250000, 0.000000, 0.250000, 0,
	0.000000, 0.250000, 0.000000, 0,
	0.250000, 0.250000, -0.250000, 0
};
float __declspec(align(16)) x_l[16] = 
{
	0.125000, 0.125000, 0.125000, 0,
	-0.125000, 0.125000, -0.125000, 0,
	0.125000, -0.125000, -0.125000, 0,
	-0.125000, -0.125000, 0.125000, 0
};
float *x_cp = x_c;
float *x_lp = x_l;

void calc_coeffs(float *coeffs, float* v, int i, int j)
{
	__asm {
		mov		eax, v
		mov		ebx, coeffs

		movaps	xmm0, [eax] // v (v0, v1)
		movaps xmm1, xmm0
		shufps	xmm0, xmm0, 0x04 // v0
		shufps	xmm1, xmm1, 0xae // v1

		movaps	xmm7, [ebx]

		mov		eax, j
		add		eax, eax
		add		eax, i
		shl		eax, 4
		mov		ecx, x_cp
		add		ecx, eax
		mov		edx, x_lp
		add		edx, eax

		movaps	xmm4, [x_neg]
		movaps	xmm5, [ecx]
		movaps	xmm6, [edx]

		// xmm2 = lin
		movaps	xmm2, xmm1
		addps	xmm2, xmm0

		// xmm3 = norm
		movaps	xmm3, xmm1
		subps	xmm3, xmm0
		shufps	xmm3, xmm3, 0xd1
		mulps	xmm3, xmm4

		// calc coeffs
		mulps	xmm6, xmm2
		addps	xmm5, xmm6
		mulps	xmm5, xmm3
		addps	xmm5, xmm7

		// store result
		movaps	[ebx], xmm5
	}
}
#else

// specialized code is about 5-10% faster with intel compiler, slower with MS compiler
__forceinline void calc_coeffs00(float *coeffs, float* v)
{
	float norm[2] = {(v[3] - v[1])*.125, (v[0] - v[2])*.125};
	float lin[2] = {v[0] + v[2], v[1] + v[3]};

	coeffs[1 - 1] += (lin[0])*norm[0]; // 1,0
	coeffs[3 - 1] += (lin[0])*norm[0]; // 1,1
	coeffs[2 - 1] += (lin[1])*norm[1]; // 0,1
}

__forceinline void calc_coeffs01(float *coeffs, float* v)
{
	float norm[2] = {(v[3] - v[1])*.125, (v[0] - v[2])*.125};
	float lin[2] = {v[0] + v[2], v[1] + v[3]};

	coeffs[1 - 1] += (lin[0])*norm[0]; // 1,0
	coeffs[3 - 1] -= (lin[0])*norm[0]; // 1,1
	coeffs[2 - 1] += (2 - lin[1])*norm[1]; // 0,1
}

__forceinline void calc_coeffs10(float *coeffs, float* v)
{
	float norm[2] = {(v[3] - v[1])*.125, (v[0] - v[2])*.125};
	float lin[2] = {v[0] + v[2], v[1] + v[3]};

	coeffs[1 - 1] += (2 - lin[0])*norm[0]; // 1,0
	coeffs[3 - 1] += (2 - lin[0])*norm[0]; // 1,1
	coeffs[2 - 1] += (lin[1])*norm[1]; // 0,1
}

__forceinline void calc_coeffs11(float *coeffs, float* v)
{
	float norm[2] = {(v[3] - v[1])*.125, (v[0] - v[2])*.125};
	float lin[2] = {v[0] + v[2], v[1] + v[3]};

	coeffs[1 - 1] += (2 - lin[0])*norm[0]; // 1,0
	coeffs[3 - 1] -= (2 - lin[0])*norm[0]; // 1,1
	coeffs[2 - 1] += (2 - lin[1])*norm[1]; // 0,1
}

__forceinline void calc_coeffs_line(float *coeffs, float* v, int i, int j)
{
	if ( i == 0 )
	{
		if ( j == 0 )
		{
			calc_coeffs00(coeffs,v);
		}
		else
		{
			calc_coeffs01(coeffs,v);
		}
	}
	else
	{
		if ( j == 0 )
		{
			calc_coeffs10(coeffs,v);
		}
		else
		{
			calc_coeffs11(coeffs,v);
		}
	}
}
/*/


const float C1[2] = {0,.25};
const float L1[2] = {1*.125,-1*.125};

const float C2[2][2] = {{0,0},{1*.25,-1*.25}};
const float L2[2][2] = {{1*.125,-1*.125},{-1*.125,1*.125}};
__forceinline void calc_coeffs(float *coeffs, float* v, int i, int j)
{
	float norm[2] = {v[3] - v[1], v[0] - v[2]};
	float lin[2] = {(v[0] + v[2]), (v[1] + v[3])};

	coeffs[1 - 1] += (C1[i] + L1[i]*lin[0])*norm[0]; // 1,0
	coeffs[2 - 1] += (C1[j] + L1[j]*lin[1])*norm[1]; // 0,1

	coeffs[3 - 1] += (C2[i][j] + L2[i][j]*lin[0])*norm[0]; // 1,1
}*/
#endif

void insert_line(Node *n, Line &p, int depth);
void insert_line_leaf(NodeData *n, Line &p);
void insert_line_split(Node *n, Line &p, int depth);
void insert_line_split_leaf(NodeData *n, Line &p);

__forceinline void insert_line_split_leaf(NodeData *n, Line &p)
{
	// scale by a factor of 2 (integer trick, add one to exponent)
#if 0
	{
		int *iter = (int*)p.v;
		int *end = iter + 2 * 2;
		for (; iter != end; iter++)
			*iter += 0x800000;
	}
#else
	p[0] *= 2;
	p[1] *= 2;
#endif

	// split line
	__declspec(align(16)) Line s[4];
	bool sin[4];

	// do full split
	__declspec(align(16)) Line x[2];
	bool xin[2];

	split_line(p, x[0], x[1], xin[0], xin[1], 0);
	for (int i = 0; i < 2; i++)
	{
		if (xin[i] == false)
			continue;
		split_line(x[i], s[i], s[i+2], sin[i], sin[i+2], 1);

		for (int j = 0; j < 2; j++)
		{
			int ij = i + j*2;
			if (sin[ij] == false)
				continue;

			// add to coefficients
			vect2f &v0 = s[ij].v[0];
			calc_coeffs_line(n->coeffs, v0.v, i, j);
			n->is_boundary[ij] = true;
		}
	}
}

void insert_line_split(Node *n, Line &p, int depth)
{
	// scale by a factor of 2 (integer trick, add one to exponent)
#if 0
	{
		int *iter = (int*)p.v;
		int *end = iter + 2 * 2;
		for (; iter != end; iter++)
			*iter += 0x800000;
	}
#else
	p[0] *= 2;
	p[1] *= 2;
#endif

	// split line
	__declspec(align(16)) Line s[4];
	bool sin[4];

	// do full split
	__declspec(align(16)) Line x[2];
	bool xin[2];

	split_line(p, x[0], x[1], xin[0], xin[1], 0);
	for (int i = 0; i < 2; i++)
	{
		if (xin[i] == false)
			continue;
		split_line(x[i], s[i], s[i+2], sin[i], sin[i+2], 1);

		for (int j = 0; j < 2; j++)
		{
			int ij = i + j*2;
			if (sin[ij] == false)
				continue;

			// add to coefficients
			vect2f &v0 = s[ij].v[0];
			calc_coeffs_line(n->data.coeffs, v0.v, i, j);
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

				insert_line(child, s[ij], depth + 1);
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

				insert_line_leaf(child, s[ij]);
			}
		}
	}
}

void insert_line_leaf(NodeData *n, Line &p)
{
	// check if all in one cell
	int i = p.v[0].v[0] < .5 ? 0 : 1;
	int j = p.v[0].v[1] < .5 ? 0 : 1;

	{
		int ni = p.v[1].v[0] < .5 ? 0 : 1;
		int nj = p.v[1].v[1] < .5 ? 0 : 1;
		if (i != ni || j != nj)
		{
			insert_line_split_leaf(n, p);
			return;
		}
	}

	// put tri in [0,1]
	for (int q = 0; q < 2; q++)
	{
		p.v[q].v[0] = p.v[q].v[0]*2 - i;
		p.v[q].v[1] = p.v[q].v[1]*2 - j;
	}

	// add to coefficients
	vect2f &v0 = p.v[0];
	calc_coeffs_line(n->coeffs, v0.v, i, j);
	n->is_boundary[i + j*2] = true;
}

void insert_line(Node *n, Line &p, int depth)
{
	for (;;)
	{
		// check if all in one cell
		int i = p.v[0].v[0] < .5 ? 0 : 1;
		int j = p.v[0].v[1] < .5 ? 0 : 1;
		int ij = (j<<1) | i; // for whatever reason this trick only seems to give gains if used here
		
		int ni = p.v[1].v[0] < .5 ? 0 : 1;
		int nj = p.v[1].v[1] < .5 ? 0 : 1;
		int nij = (nj<<1) | ni;
		if (ij != nij)
		{
			insert_line_split(n, p, depth);
			return;
		}

		// put tri in [0,1]
		for (int q = 0; q < 2; q++)
		{
			p.v[q].v[0] = p.v[q].v[0]*2 - i;
			p.v[q].v[1] = p.v[q].v[1]*2 - j;
		}

		// add to coefficients
		vect2f &v0 = p.v[0];
		calc_coeffs_line(n->data.coeffs, v0.v, i, j);
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

			n = child;
			depth++;
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

			insert_line_leaf(child, p);
			return;
		}
	}
}

void insert_line_root(Line &p)
{
	// calculate root coefficient
	g.coeff_const += (p[0][0]*p[1][1] - p[0][1]*p[1][0]) * .5;
	
	// add line to tree
	insert_line(g.root, p, 0);
}