#pragma once
#include "bucket_allocator.h"
#include "boundedarray.h"

struct NodeData
{
	void init()
	{
#if 0
		for (int i = 0; i < 3; i++)
			coeffs[i] = 0;
		for (int i = 0; i < 4; i++)
			is_boundary[i] = false;
#else
		int *ptr = (int*)&coeffs[0];
		for (int i = 0; i < sizeof(NodeData)/4; i++)
			ptr[i] = 0;
#endif
	}
	float coeffs[3];
	bool is_boundary[4]; // store if there is an intersecting line
};

// octree node that uses my custom allocator. is 64 bytes long, so fits exactly in one cache line (or is at least cache aligned)
struct Node
{
	void init()
	{
#if 0
		for (int i = 0; i < 4; i++)
			children[i].table = children[i].offset = 0;
		data.init();
#else
		int *ptr = (int*)&children[0];
		for (int i = 0; i < sizeof(Node)/4; i++)
			ptr[i] = 0;
#endif
	}

	BucketPointer children[4];
	NodeData data;
};

