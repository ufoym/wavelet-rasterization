#include <stdio.h>
#include "global.h"
#include "raster_poly.h"
#include "write_grid.h"
#include "timer.h"
#include "vect.h"
#include "get_font.h"
#include "save.h"
#include <algorithm>

using namespace std;

int iters = 1;

void raster_star_noisy(int depth, int arms)
{
	srand(10202);

	Array<Line> lines;
	Array<Bez2> bez2s;

	// create poly
	const int sides = arms*2;
	const double pi = 3.1415926535;
	
	float noise = 1e-6;
	for (int i = 0; i < sides; i++)
	{
		float t = i * (1.0 / sides) * 2 * pi + pi/3;
		float tp = (i + 1) * (1.0 / sides) * 2 * pi + pi/3;

		float r = .3;
			
		Line line;
		line[0].set(cos(t) * r + .5 + rand()*noise, sin(t) * r + .5 + rand()*noise);
		line[1].set(cos(tp) * r + .5 + rand()*noise, sin(tp) * r + .5 + rand()*noise);
		lines.push(line);
	}
	
	// save obj for other rasterizers
	save_obj(lines);

	// create grid
	g.max_depth = depth - 1;
	g.grid_res = 1 << (g.max_depth + 1);

	Array2D<float> grid;
	grid.resize(g.grid_res, g.grid_res);

	for (int k = 0; k < g.grid_res*g.grid_res; k++)
		grid.data[k] = 0;

	// rasterize
	double t_start = get_time();
	for (int i = 0; i < iters; i++)
		raster_poly(lines, bez2s, grid.data);
	double t_end = get_time();
	printf("time wavelet = %f\n", (t_end - t_start) / iters);

	// save timing
	FILE *f = fopen("timings.txt", "a");
	fprintf(f, "time wavelet = %f\n", (t_end - t_start) / iters);
	fclose(f);

	// save image
	Array2D<vect3ub> ourgrid;
	ourgrid.resize(g.grid_res, g.grid_res);

	for (int k = 0; k < g.grid_res*g.grid_res; k++)
	{
		if (grid.data[k] > 1)
			grid.data[k] = 1;
		else if (grid.data[k] < 0)
			grid.data[k] = 0;

		ourgrid.data[k] = (1-grid.data[k])*255;
	}

	save_png("output.png", ourgrid);
}

void raster_star(int depth, int arms, float R1 = .19, float R2 = .5)
{
	Array<Line> lines;
	Array<Bez2> bez2s;

	// create poly
	const int sides = arms*2;
	const double pi = 3.1415926535;
	
	for (int i = 0; i < sides; i++)
	{
		float t = i * (1.0 / sides) * 2 * pi + pi/2;
		float tp = (i + 1) * (1.0 / sides) * 2 * pi + pi/2;

		float rad1 = R1, rad2 = R2;
		if ((i % 2) == 0)
			swap(rad1, rad2);
			
		Line line;
		line[0].set(cos(t) * rad1 + .5, sin(t) * rad1 + .5);
		line[1].set(cos(tp) * rad2 + .5, sin(tp) * rad2 + .5);
		lines.push(line);
	}

	// save obj for other rasterizers
	save_obj(lines);

	// create grid
	g.max_depth = depth - 1;
	g.grid_res = 1 << (g.max_depth + 1);

	Array2D<float> grid;
	grid.resize(g.grid_res, g.grid_res);

	for (int k = 0; k < g.grid_res*g.grid_res; k++)
		grid.data[k] = 0;

	// rasterize
	double t_start = get_time();
	for (int i = 0; i < iters; i++)
		raster_poly(lines, bez2s, grid.data);
	double t_end = get_time();
	printf("time wavelet = %f\n", (t_end - t_start) / iters);

	// save timing
	FILE *f = fopen("timings.txt", "a");
	fprintf(f, "time wavelet = %f\n", (t_end - t_start) / iters);
	fclose(f);

	// save image
	Array2D<vect3ub> ourgrid;
	ourgrid.resize(g.grid_res, g.grid_res);

	for (int k = 0; k < g.grid_res*g.grid_res; k++)
		ourgrid.data[k] = (1-grid.data[k])*255;

	save_png("output.png", ourgrid);
}


void raster_spiral(int depth, int arms, int segs)
{
	Array<Line> lines;
	Array<Bez2> bez2s;

	g.max_depth = depth - 1;
	g.grid_res = 1 << (g.max_depth + 1);

	// create poly
	const int sides = arms*2;
	const double pi = 3.1415926535;
	
	for (int i = 0; i < sides; i++)
	{
		float t1 = i * (1.0 / sides) * 2 * pi + pi/2;
		float t2 = (i + 1) * (1.0 / sides) * 2 * pi + pi/2;

		float r1 = .1;
		float r2 = .5;
		
		float toff = 2*pi;
		if ((i % 2) == 0)
		{
			t1 += toff;
			swap(r1, r2);
		}
		else
		{
			t2 += toff;
		}

		for (int j = 0; j < segs; j++)
		{
			float s1 = j * (1.0 / segs);
			float s2 = (j + 1) * (1.0 / segs);
			
			float tt1 = t1 * (1-s1) + t2 * s1;
			float tt2 = t1 * (1-s2) + t2 * s2;

			float rr1 = r1 * (1-s1) + r2 * s1;
			float rr2 = r1 * (1-s2) + r2 * s2;

			Line line;
			line[0].set(cos(tt1) * rr1 + .5, sin(tt1) * rr1 + .5);
			line[1].set(cos(tt2) * rr2 + .5, sin(tt2) * rr2 + .5);
			lines.push(line);
		}
	}

	// save obj for other rasterizers
	save_obj(lines);

	// create grid
	Array2D<float> grid;
	grid.resize(g.grid_res, g.grid_res);

	for (int k = 0; k < g.grid_res*g.grid_res; k++)
		grid.data[k] = 0;

	// rasterize
	double t_start = get_time();
	for (int i = 0; i < iters; i++)
		raster_poly(lines, bez2s, grid.data);
	double t_end = get_time();
	printf("time wavelet = %f\n", (t_end - t_start) / iters);

	// save timing
	FILE *f = fopen("timings.txt", "a");
	fprintf(f, "time wavelet = %f\n", (t_end - t_start) / iters);
	fclose(f);

	// save image
	Array2D<vect3ub> ourgrid;
	ourgrid.resize(g.grid_res, g.grid_res);

	for (int k = 0; k < g.grid_res*g.grid_res; k++)
		ourgrid.data[k] = (1-grid.data[k])*255;

	save_png("output.png", ourgrid);
}


void raster_spiral_out_of_core(int depth, int arms, int segs)
{
	g.max_depth = depth - 1;
	g.grid_res = 1 << (g.max_depth + 1);

	g.root = g.node_allocator.get_ref(g.node_allocator.alloc());
	g.root->init();
	g.coeff_const = 0;

	// create poly
	const int sides = arms*2;
	const double pi = 3.1415926535;
	
	for (int i = 0; i < sides; i++)
	{
		float t1 = i * (1.0 / sides) * 2 * pi + pi/2;
		float t2 = (i + 1) * (1.0 / sides) * 2 * pi + pi/2;

		float r1 = .1;
		float r2 = .5;
		
		float toff = 2*pi;
		if ((i % 2) == 0)
		{
			t1 += toff;
			swap(r1, r2);
		}
		else
		{
			t2 += toff;
		}

		for (int j = 0; j < segs; j++)
		{
			float s1 = j * (1.0 / segs);
			float s2 = (j + 1) * (1.0 / segs);
			
			float tt1 = t1 * (1-s1) + t2 * s1;
			float tt2 = t1 * (1-s2) + t2 * s2;

			float rr1 = r1 * (1-s1) + r2 * s1;
			float rr2 = r1 * (1-s2) + r2 * s2;

			Line line;
			line[0].set(cos(tt1) * rr1 + .5, sin(tt1) * rr1 + .5);
			line[1].set(cos(tt2) * rr2 + .5, sin(tt2) * rr2 + .5);
			insert_line_root(line);
		}
	}

	// create grid
	Array2D<float> grid;
	grid.resize(g.grid_res, g.grid_res);

	for (int k = 0; k < g.grid_res*g.grid_res; k++)
		grid.data[k] = 0;

	// rasterize
	//raster_poly(lines, bez2s, grid.data);
	write_grid(grid.data);

	g.node_allocator.clear();
	g.leaf_allocator.clear();

	// save image
	Array2D<vect3ub> ourgrid;
	ourgrid.resize(g.grid_res, g.grid_res);

	for (int k = 0; k < g.grid_res*g.grid_res; k++)
		ourgrid.data[k] = (1-grid.data[k])*255;

	save_png("output.png", ourgrid);
}

int main(int argc, char **argv)
{
	// defaults
	int px = 512;
	int poly_res = 500000;
	string method = "c";
	string font = "vladimir";

	// read input
	if (argc > 1)
	{
		px = atoi(argv[1]);

		poly_res = atoi(argv[2]);

		method = argv[3];

		if (method == "f")
		{
			font = argv[4];
		}
	}
	
	int depth = log((double)px)/log(2.0) + .5;

	// raster
	if (method == "n")
		raster_star_noisy(depth, poly_res);
	else if (method == "st")
		raster_star(depth, poly_res);
	else if (method == "c")
		raster_star(depth, poly_res, .45, .45);
	else if (method == "sp")
		raster_spiral(depth, poly_res, 10000);
	else if (method == "f")
	{
		vect2d times;
		times = 0;

		for (char c = 33; c < 127; c++)
			get_font(font, c, px, &times);

		printf("time to raster with FT   = %f\n", times[0]);
		printf("time to raster with Ours = %f\n", times[1]);
	}
}
