#pragma once
#include <assert.h>

//#define USE_GC

#ifdef USE_GC
#include "gc_cpp.h"
#endif

template <class T, bool POINTER_FREE = false>
struct Array2D
{
	Array2D(){data = 0; clear();}
	Array2D(const Array2D<T> &a){data = 0; copy(a);}
#ifdef USE_GC
#else
	~Array2D(){clear();}
#endif

	T *data;
	int data_size;
	int size[2];

	void operator=(const Array2D<T> &a) {copy(a);}
	void copy(const Array2D<T> &a)
	{
		resize(a.size[0], a.size[1]);
		for (int i = 0; i < data_size; i++)
			data[i] = a.data[i];
	}

	void resize(int w, int h)
	{
		clear();

		data_size = w*h;
		size[0] = w;
		size[1] = h;

#ifdef USE_GC
		if (POINTER_FREE)
			data = new (PointerFreeGC) T[data_size];
		else
			data = new (UseGC) T[data_size];
#else
		data = new T[data_size];
#endif
	}

	T &operator()(int x, int y)
	{
		return data[y*size[0] + x];
	}

	void clear()
	{
#ifdef USE_GC
#else
		delete[] data;
#endif
		data = 0;
		size[0] = size[1] = 0;
		data_size = 0;
	}
};
