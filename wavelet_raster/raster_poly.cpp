#include "raster_poly.h"
#include <windows.h>
#include "timer.h"
#include "vect.h"
#include "global.h"
#include <stdlib.h>
#include "array.h"
#include "insert_line.h"
#include "insert_bez2.h"
#include "write_grid.h"


void calc_coeffs_test()
{
	// add lines into tree
	//double time_start = get_time();
	g.coeff_const = 0;

	Line l;
	Bez2 q;
	
	l[0].set(0, 0);
	l[1].set(1, 0);
	insert_line_root(l);
	
	l[0].set(0, 1);
	l[1].set(0, 0);
	insert_line_root(l);

	if (0)
	{
	l[0].set(1, 0);
	l[1].set(0, 1);
	insert_line_root(l);
	}
	else
	{
	q[0].set(1, 0);
	q[1].set(.75, .75);
	q[2].set(0, 1);
	insert_bez2_root(q);
	}

	//double time_end = get_time();
	//printf("time to create tree and calculate coefficients = %fs\n", time_end - time_start);
}

void calc_coeffs(Array<Line> &lines, Array<Bez2> &bez2s)
{
	// add lines into tree
	//double time_start = get_time();
	g.coeff_const = 0;

	for (int i = 0; i < lines.s; i++)
	{
		Line l = lines[i];
		insert_line_root(l);
	}

	for (int i = 0; i < bez2s.s; i++)
	{
		Bez2 b = bez2s[i];
		insert_bez2_root(b);
	}

	//double time_end = get_time();
	//printf("time to create tree and calculate coefficients = %fs\n", time_end - time_start);
}


void raster_poly(Array<Line> &lines, Array<Bez2> &bez2s, float *grid)
{
	//double time_start = get_time();

	g.root = g.node_allocator.get_ref(g.node_allocator.alloc());
	g.root->init();

	//calc_coeffs_test();
	calc_coeffs(lines, bez2s);
	write_grid(grid);

	g.node_allocator.clear();
	g.leaf_allocator.clear();
	
	//double time_end = get_time();
	//printf("time to raster poly = %fs\n", time_end - time_start);
}
