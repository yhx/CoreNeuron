

#include "extent.h"

__device__ int commit2globalTable(int *shared_buf, volatile unsigned int size, int *global_buf, int * global_size, int offset) 
{
	__shared__ volatile unsigned int start_loc;
	if (threadIdx.x == 0) {
		start_loc = atomicAdd(global_size, (int)size);
	}
	__syncthreads();

	for (int idx=threadIdx.x; idx<size; idx+=blockDim.x) {
		global_buf[offset + start_loc + idx] = shared_buf[idx];
	}

	return 0;
}


__global__  record_netcon(NetCon *n, PreSyn *pre, int *nids, int size) {

	double time_a[MAX_BLOCK_SIZE];
	NetCon 
	int block_num = gridDim.x;
	int num_per_block = (idx_size-1+block_num)/block_num;
	int num_per_block_1 = num_per_block - 1;
	int offset = size - block_num * num_per_block_1;

	int idx_size_t = 0;
	int idx_offset = 0;
	int block_idx = blockIdx.x;

	if (block_idx < offset) {
		idx_size_t = num_per_block;
		idx_offset = block_idx * num_per_block;
	} else if (block_idx < block_num) {
		idx_size_t = num_per_block_1;
		idx_offset = offset * num_per_block + (block_idx-offset) * num_per_block_1;
	} else {
		idx_size_t = 0;
	}

	int tid = threadIdx.x;
	int thread_num = blockDim.x;

	for  (int i=0; i < idx_size_t; i++) {
		PreSyn *p = pre + nids[i];
		for (int j=tid; j<p->nc_cnt_; j+=thread_num) {
			NetCon *d = netcon_in_presyn_order[p->nc_index_+j];
			if (d->active_ && d->target_) {


			}
		}
	}
}
