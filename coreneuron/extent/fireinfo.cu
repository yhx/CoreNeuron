
#include "fireinfo"

__device__ int commitFireInfo(FireInfo *shared_buf, volatile unsigned int size, FireInfo *global_buf, int *global_size, int offset) 
{
	__shared__ volatile unsigned int start_loc;
	if (threadIdx.x == 0) {
		start_loc = atomicAdd(global_size, (int)size);
	}
	__syncthreads();

	for (int idx=threadIdx.x; idx<size; idx+=blockDim.x) {
		global_buf[offset + start_loc + idx].n = shared_buf[idx].n;
		global_buf[offset + start_loc + idx].time = shared_buf[idx].time;
	}

	return 0;
}

