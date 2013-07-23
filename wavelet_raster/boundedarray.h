#pragma once

#include "assert.h"

template <class T, int capacity = 9>
struct BoundedArray
{
	int s; // size
	T elem[capacity];

	void init(){s = 0;}
	void push_back(T &e){elem[s++] = e;	assert(s <= capacity);}
	T &operator[](int i){return elem[i]; assert(i >= 0 && i < s);}
	int size(){return s;}
	void resize(int s){this->s = s; assert(s <= capacity);}
	void clear(){s = 0;}
	void reserve(int s){}
};
