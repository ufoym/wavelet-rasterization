#include "write_grid.h"
#include "global.h"
#include "timer.h"
#include "vect.h"

#include <windows.h>

//__forceinline unsigned char get_char_val(float val)
//{
//	return (unsigned char)(val * 255);
//}

void write_value_to_grid(float val, vect2i off, int res, float *grid)
{
	// don't draw black pixels
	/*if (val < (1.0 / 512.0))
		return; */

	// get the pixel val
	//unsigned char charval = get_char_val(val);

	float *py = grid + off[1] * g.grid_res + off[0];
	float *endy = py + g.grid_res * res;
	for (; py != endy; py += g.grid_res)
	{
		float *px = py;
		float *endx = px + res;
		for (; px != endx; px++)
			*px = val;
	}
}

void write_node_to_grid_leaf(NodeData *n, float val, vect2i offset, int res, float *grid)
{
	float cvals[4];

	// calc values
	for (int k = 0; k < 4; k++)
	{
		cvals[k] = val;
		for (int i = 0; i < 3; i++)
		{
			const int ii = (i + 1) & k;
			const int x = ii & 1;
			const int y = (ii >> 1) & 1;
			const int sign = 1 - 2*(x ^ y);
			cvals[k] += sign * n->coeffs[i];
		}
	}

	// process all children
	for (int k = 0; k < 4; k++)
	{
		const int x = k & 1;
		const int y = (k >> 1) & 1;

		// leaf, so write value to disk
		if (n->is_boundary[k])
		{
			float &px = grid[(offset[1] + y) * g.grid_res + (offset[0] + x)];
			px = cvals[k];
		}
		else if (cvals[k] > .5)
		{
			float &px = grid[(offset[1] + y) * g.grid_res + (offset[0] + x)];
			px = 1;
		}
	}
}

void write_node_to_grid(Node *n, float val, vect2i offset, int res, float *grid)
{
	float cvals[4];

	// calc values
	for (int k = 0; k < 4; k++)
	{
		cvals[k] = val;
		for (int i = 0; i < 3; i++)
		{
			const int ii = (i + 1) & k;
			const int x = ii & 1;
			const int y = (ii >> 1) & 1;
			const int sign = 1 - 2*(x ^ y);
			cvals[k] += sign * n->data.coeffs[i];
		}
	}

	// process all children
	int res2 = res / 2;

	for (int k = 0; k < 4; k++)
	{
		const int x = k & 1;
		const int y = (k >> 1) & 1;

		vect2i off;
		off[0] = offset[0] + x*res2;
		off[1] = offset[1] + y*res2;

		const float cval = cvals[k];

		if (n->children[k].invalid())
		{
			// leaf, so write values to disk
			if (n->data.is_boundary[k])
				write_value_to_grid(cval, off, res2, grid);
			else if (cval > .5)
				write_value_to_grid(1, off, res2, grid); // this should be optimized
		}
		else if (res > 4)
		{
			// internal node, so recur
			Node *cn = g.node_allocator.get_ref(n->children[k]);
			write_node_to_grid(cn, cval, off, res2, grid);
		}
		else
		{
			// recur on leaf node
			NodeData *cn = g.leaf_allocator.get_ref(n->children[k]);
			write_node_to_grid_leaf(cn, cval, off, res2, grid);
		}
	}
}

void write_grid(float *grid)
{
	//double time_start = get_time();

	vect2i offset; offset.set(0);
	write_node_to_grid(g.root, g.coeff_const, offset, g.grid_res, grid);

	//double time_save = get_time();
	//printf("time to convert coefficients to grid = %fs\n", time_save - time_start);
}
