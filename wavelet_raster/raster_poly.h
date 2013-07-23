#pragma once

#include <string>
#include  "insert_line.h"
#include  "insert_bez2.h"
#include  "array.h"
#include "vect.h"

using namespace std;

void calc_coeffs(Array<Line> &lines, Array<Bez2> &bez2s);
void raster_poly(Array<Line> &lines, Array<Bez2> &bez2s, float *grid);