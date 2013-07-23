#include <string>
#include <iostream>
#include <algorithm>
#include <direct.h>
#include  <freetype/ftglyph.h>
#include  <freetype/freetype.h>
#include  <freetype/ftoutln.h>
#include  <freetype/ftbbox.h>
#include "get_font.h"
#include "save.h"
#include "global.h"
#include "raster_poly.h"
#include "timer.h"

extern int iters;

//******************* check error code ********************
void Check(FT_Error ErrCode, const char* OKMsg, const char* ErrMsg)
{
	if(ErrCode != 0)
	{
		std::cout << ErrMsg << ": " << ErrCode << "\n";
		std::cout << "program halted\n";
		exit(1);
	}
}

//******************** get outline ************************
int GetOutLine(FT_Glyph glyph, FT_OutlineGlyph* Outg)
{
	int Err = 0;

	switch ( glyph->format )
	{
	case FT_GLYPH_FORMAT_BITMAP:
		Err = 1;
		break;

	case FT_GLYPH_FORMAT_OUTLINE:
		*Outg = (FT_OutlineGlyph)glyph;
		break;

	default:
		;
	}
	return Err;
}

Array<Line> lines;
Array<Bez2> bez2s;
vect2f cursor;

int move_to( const FT_Vector*  to, void* user)
{
	cursor.set(to->x/64., to->y/64.);
	return 0;
}

int line_to( const FT_Vector*  to, void* user)
{
	Line l;
	l[0] = cursor;
	cursor.set(to->x/64., to->y/64.);
	l[1] = cursor;

	std::swap(l[0], l[1]);
	lines.push(l);
	return 0;
}

int conic_to( const FT_Vector*  ctrl,  const FT_Vector*  to, void* user)
{
	Bez2 l;
	l[0] = cursor;
	l[1].set(ctrl->x/64., ctrl->y/64.);
	cursor.set(to->x/64., to->y/64.);
	l[2] = cursor;

	std::swap(l[0], l[2]);
	bez2s.push(l);
	return 0;
}

int cubic_to( const FT_Vector*  ctrl1,  const FT_Vector*  ctrl2,  const FT_Vector*  to, void* user)
{
	printf("cubic_to is unsupported\n");
	cursor.set(to->x/64., to->y/64.);
	return 0;
}

//******************** M A I N ************************
void get_font(string FontName_base, char letter, int ptsize, vect2d *times)
{
	_mkdir("font_rasters");
	std::string FontName = "fonts/" + FontName_base + ".ttf";
	char buf[256];
	string letter_name;
	letter_name += letter;
	if (letter >= 'A' && letter <= 'Z')
		letter_name += letter;

	//-----------------------------------------
	// Load contour of glyph
	//-----------------------------------------

	lines.clear();
	bez2s.clear();

	FT_Face face;
	FT_Library    library;
	FT_Error error;

	error = FT_Init_FreeType( &library );
	Check(error, "", "error initializing FT lib");

	error = FT_New_Face( library, FontName.c_str(), 0, &face );   
	Check(error, "",  "error loading font");

	FT_F26Dot6 font_size = ptsize*64;
	error = FT_Set_Char_Size( face, font_size, font_size, 72, 72 );
	Check(error, "", "error setting char size");

	FT_UInt   glyph_index = FT_Get_Char_Index( face, letter ); 
	FT_Int32  load_flags = FT_LOAD_NO_HINTING | FT_LOAD_NO_BITMAP;   
	error = FT_Load_Glyph( face,  glyph_index, load_flags );      
	Check(error, "", "error loading glyph");

	FT_Glyph glyph; 
	error = FT_Get_Glyph( face->glyph, &glyph );  
	Check(error, "", "error getting glyph");

	FT_OutlineGlyph Outg;
	error = GetOutLine(glyph, &Outg);
	Check(error,"", "error getting outline"); 

	// use my own callcacks to walk over the contour
	FT_Outline_Funcs func_interface;
	func_interface.shift = 0;
	func_interface.delta = 0;
	func_interface.move_to = move_to;
	func_interface.line_to = line_to;
	func_interface.conic_to = conic_to;
	func_interface.cubic_to = cubic_to;
	FT_Outline_Decompose(&Outg->outline, &func_interface, 0);

	//-----------------------------------------
	// Rasterize using FreeType
	//-----------------------------------------

	// rasterize using FreeType so that a comparison to my code can be made
	double t_ft_start = get_time();
	FT_Render_Glyph(face->glyph, FT_RENDER_MODE_NORMAL);
	double t_ft_stop = get_time();

	// save timing
	FILE *fa = fopen("timings.txt", "a");
	fprintf(fa, "time FreeType = %f\n", t_ft_stop - t_ft_start);
	fclose(fa);

	// create grid / calculate resolution
	g.max_depth = 0;
	g.grid_res = 2;

	int maxsize = max(face->glyph->bitmap.width, face->glyph->bitmap.rows);
	while (g.grid_res < maxsize)
	{
		g.grid_res *= 2;
		g.max_depth++;
	}

	vect2i offset;
	offset.set((g.grid_res - face->glyph->bitmap.width) / 2, (g.grid_res - face->glyph->bitmap.rows) / 2);

	// copy bitmap into a buffer
	Array2D<vect3ub> ftgrid; 
	ftgrid.resize(g.grid_res, g.grid_res);
	
	for (int k = 0; k < g.grid_res*g.grid_res; k++)
		ftgrid.data[k].set(0);

	for (int j = 0; j < face->glyph->bitmap.rows; j++)
	{
		unsigned char *row = face->glyph->bitmap.buffer + (face->glyph->bitmap.rows-j-1)*face->glyph->bitmap.pitch;

		for (int i = 0; i < face->glyph->bitmap.width; i++)
		{
			ftgrid(i + offset[0], j + offset[1]) = row[i];
		}
	}

	sprintf(buf, "font_rasters/%s_%04d_%s_1_ft.png", FontName_base.c_str(), ptsize, letter_name.c_str());
	save_png(buf, ftgrid);

	//-----------------------------------------
	// Rasterize using our method
	//-----------------------------------------

	// get bbox
	FT_BBox bbox;
	FT_Outline_Get_BBox(&Outg->outline, &bbox);
	//printf("bbox = (%f, %f), (%f, %f)\n", bbox.xMin/64., bbox.yMin/64., bbox.xMax/64., bbox.yMax/64.);

	// fit in box
	vect2f ext; ext.set(bbox.xMax/64. - bbox.xMin/64., bbox.yMax/64. - bbox.yMin/64.);
	float maxext = std::max(ext[0], ext[1]) * 1.1;

	for (int i = 0; i < lines.s; i++)
	{
		//printf("line\n");
		for (int j = 0; j < 2; j++)
		{
			lines[i][j][0] = (lines[i][j][0] - floor(bbox.xMin/64.) + offset[0]) / g.grid_res;
			lines[i][j][1] = (lines[i][j][1] - floor(bbox.yMin/64.) + offset[1]) / g.grid_res;
			//printf("%f %f\n", lines[i][j][0], lines[i][j][1]);
		}
	}
	for (int i = 0; i < bez2s.s; i++)
	{
		//printf("bez2\n");
		for (int j = 0; j < 3; j++)
		{
			bez2s[i][j][0] = (bez2s[i][j][0] - floor(bbox.xMin/64.) + offset[0]) / g.grid_res;
			bez2s[i][j][1] = (bez2s[i][j][1] - floor(bbox.yMin/64.) + offset[1]) / g.grid_res;
			//printf("%f %f\n", bez2s[i][j][0], bez2s[i][j][1]);
		}
	}

	//vect3ub *grid = new vect3ub [g.grid_res*g.grid_res];
	Array2D<float> grid;
	Array2D<vect3ub> ourgrid;
	grid.resize(g.grid_res, g.grid_res);
	ourgrid.resize(g.grid_res, g.grid_res);

	// rasterize
	for (int k = 0; k < g.grid_res*g.grid_res; k++)
		grid.data[k] = 0;
	double t_ours_start = get_time();
	raster_poly(lines, bez2s, grid.data);
	double t_ours_stop = get_time();

	// save timing
	FILE *f = fopen("timings.txt", "a");
	fprintf(f, "time wavelet = %f\n", (t_ours_stop - t_ours_start));
	fclose(f);

	// copy into color image and save
	for (int k = 0; k < g.grid_res*g.grid_res; k++)
	{
		if (grid.data[k] >= 1)
			ourgrid.data[k] = 255;
		else if (grid.data[k] <= 0)
			ourgrid.data[k] = 0;
		else
			ourgrid.data[k] = grid.data[k]*255;
	}
	sprintf(buf, "font_rasters/%s_%04d_%s_2_ours.png", FontName_base.c_str(), ptsize, letter_name.c_str());
	save_png(buf, ourgrid);

	//-----------------------------------------
	// Calculate difference between renders
	//-----------------------------------------
	Array2D<vect3ub> grid_dif;
	grid_dif.resize(g.grid_res, g.grid_res);
	
	for (int i = 0; i < grid_dif.data_size; i++)
	{
		int dif = (grid.data[i]*255 - ftgrid.data[i][0]) * 10;
		if (dif > 255)
			dif = 255;
		else if (dif < -255)
			dif = -255;

#if 1
		grid_dif.data[i] = 255;

		if (dif < 0)
		{
			grid_dif.data[i][0] += dif;
			grid_dif.data[i][1] += dif;
		}
		else
		{
			grid_dif.data[i][1] -= dif;
			grid_dif.data[i][2] -= dif;
		}
#else
		grid_dif.data[i] = 0;

		if (dif < 0)
			grid_dif.data[i][2] -= dif;
		else
			grid_dif.data[i][0] += dif;
#endif
	}
	
	sprintf(buf, "font_rasters/%s_%04d_%s_3_dif.png", FontName_base.c_str(), ptsize, letter_name.c_str());
	save_png(buf, grid_dif);

	//printf("--== timing comparison ==--\n");
	//printf("time freetype = %f\n", t_ft_stop - t_ft_start);
	//printf("time wavelets = %f\n", t_ours_stop - t_ours_start);
	//printf("times slower = %f\n", (t_ours_stop - t_ours_start) / (t_ft_stop - t_ft_start) );

	if (times)
	{
		times->v[0] += t_ft_stop - t_ft_start;
		times->v[1] += t_ours_stop - t_ours_start;
	}
}
