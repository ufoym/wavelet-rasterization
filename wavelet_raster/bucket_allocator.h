#pragma once

#include <assert.h>
#include <malloc.h>

// pointer into custom allocator data
struct BucketPointer
{
	unsigned __int16 table, offset;
	bool invalid(){return 0 == *(int*)this;}
	//bool invalid(){return 0 == table && 0 == offset;}
};

// this allocator will allocate memory of uniform types very quickly and cannot free individual elements (up to 2^32 entries (possibly more bytes))
// allocator has a constant 64k overhead, all memory can be cleared simultaneously
template<class T, int ALIGNMENT = 64>
struct BucketAllocator
{
	BucketAllocator()
	{
		init();
	}

	int table; // which table is currently being operated on
	int count; // number of elements already allocated in current table
	T* memory[65536]; // 16 bit addressable (64k) tables of data, each of which will hold 64k entries (at end so that first entries share cache line with counters)

	T* get_ref(BucketPointer p)
	{
		return memory[p.table] + p.offset;
	}

	BucketPointer alloc()
	{
		// alloc more mem if ran out of current table
		if (count == 65536)
		{
			count = 0;
			table++;
			memory[table] = (T*)_aligned_malloc(sizeof(T)*65536, ALIGNMENT);
		}

		BucketPointer p;
		p.table = table;
		p.offset = count;
		count++;

		return p;
	}

	void clear()
	{
		for (int i = 0; i <= table; i++)
			_aligned_free(memory[i]);
		
		init();
	}

private:
	void init()
	{
		// the 0th entry is wasted to preserve null pointer semantics
		memory[0] = (T*)_aligned_malloc(sizeof(T)*65536, ALIGNMENT);
		table = 0;
		count = 1;
	}
};
