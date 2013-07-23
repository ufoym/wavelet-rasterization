#pragma once

#include "octree.h"

struct Global
{
	int max_depth; // max depth of the tree
	int grid_res; // resoultion of the grid of values (2^(depth+1)) (stride in y)
	Node *root; // pointer to root of tree
	float coeff_const; // the constant offset coefficient for phiXphiXphi (C_000)
	BucketAllocator<Node> node_allocator;
	BucketAllocator<NodeData> leaf_allocator;
};

extern Global g;
