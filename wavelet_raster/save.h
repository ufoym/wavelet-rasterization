#pragma once

#include "insert_line.h"
#include "array2d.h"
#include "array.h"
#include "vect.h"
#include <string>

using namespace std;

void save_obj(Array<Line> &lines);
void save_png(string fn, Array2D<vect3ub> &grid);