#include "save.h"
#include "FreeImagePlus.h"
#include <fstream>

using namespace std;

void save_obj(Array<Line> &lines)
{
#if 0
	FILE *f = fopen("output.obj", "wb");

	for (int i = 0; i < lines.s; i++)
	{
		fprintf(f, "v %16.15f %16.15f 0\n", lines[i][0][0], lines[i][0][1]);
		fprintf(f, "v %16.15f %16.15f 0\n", lines[i][1][0], lines[i][1][1]);

		fprintf(f, "f %d %d\n", i*2 + 1, i*2 + 2);
	}

	fclose(f);
#else
	ofstream f("output.poly", ios::binary);

	__int32 s = lines.s;
	f.write((char*)&s, 4);
	f.write((char*)lines.data, 16*s);

	f.close();
#endif
}

void save_png(string fn, Array2D<vect3ub> &grid)
{
	fipImage img(FIT_BITMAP, grid.size[0], grid.size[1], 24);

	for (int j = 0; j < grid.size[1]; j++)
	{
		vect3ub *imgline = (vect3ub*)img.getScanLine(j);
		vect3ub *gridline = (vect3ub*)(grid.data + grid.size[0]*j);

		for (int i = 0; i < grid.size[0]; i++)
		{
			vect3ub &ipx = imgline[i];
			vect3ub &gpx = gridline[i];
			ipx[0] = gpx[2];
			ipx[1] = gpx[1];
			ipx[2] = gpx[0];
		}
	}

	img.save(fn.c_str());
}