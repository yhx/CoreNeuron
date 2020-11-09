

#include "connection.h"

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

__global__  gen_spike(FireInfo *spike, int *s_size, PreSyn *pre, int *nids, int size, double time) {

	__shared__ FireInfo spike_a[MAX_BLOCK_SIZE];
	__shared__ volatile int count;

	if (threadIdx.x == 0) {
		count == 0;
	}
	__syncthreads();

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

	int count_t = 0;
	bool commit = false;

	for  (int i=0; i < idx_size_t; i++) {
		PreSyn *p = pre + nids[i];
		for (int j=tid; j<p->nc_cnt_; j+=thread_num) {
			NetCon *d = netcon_in_presyn_order[p->nc_index_+j];
			commit = d->active_ && d->target_;
			if (commit) {
				count_t = atomicAdd((int *)&count, 1);
				if (count_t < MAX_BLOCK_SIZE) {
					spike_a[count_t].time = time + d->delay_;
					spike_a[count_t].n =  d;
					commit = false;
				} else {
					commitFireInfo(spike_a, MAX_BLOCK_SIZE, spike, &s_size[0], 0);
				}
			}
			__syncthreads();

			if (threadIdx.x == 0 && count >= MAX_BLOCK_SIZE) {
				count == 0;
			}
			__syncthreads();
			if (commit) {
				count_t = atomicAdd((int *)&count, 1);
				if (count_t < MAX_BLOCK_SIZE) {
					spike_a[count_t].time = time + d->delay_;
					spike_a[count_t].n =  d;
					commit = false;
				}
			}
			__syncthreads();
			if (count > 0) {
				commitFireInfo(spike_a, MAX_BLOCK_SIZE, spike, &s_size[0], 0);
			}
			if (threadIdx.x == 0) {
				count == 0;
			}
			__syncthreads();
		}
	}
}

__global__  fire_spike(FireInfo *spike, int *s_size, FireInfo *spike1, int s1_size, FireInfo *spike2, int *s2_size, double time) {
	__shared__ FireInfo spike_a[MAX_BLOCK_SIZE];
	__shared__ FireInfo spike_b[MAX_BLOCK_SIZE];
	__shared__ volatile int count_a;
	__shared__ volatile int count_b;

	if (threadIdx.x == 0) {
		count_a == 0;
		count_b == 0;
	}

	bool commit_a = false, commit_b = false;
	int count_a_t = 0, count_b_t = 0;
	__syncthreads();
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	for (int id=tid; id<s_size[0]; id+=blockDim.x*gridDim.x) {
		commit_a = spike[id].time <= time;
		commit_b = ~commit_a;
		if (commit_a) {
			count_a_t = atomicAdd((int *)&count_a, 1);
			if (count_a_t < MAX_BLOCK_SIZE) {
				spike_a[count_a_t].time = time + d->delay_;
				spike_a[count_a_t].n =  d;
				commit_a = false;
			} else {
				commitFireInfo(spike_a, MAX_BLOCK_SIZE, spike1, &s1_size[0], 0);
			}
		}
		__syncthreads();

		if (threadIdx.x == 0 && count >= MAX_BLOCK_SIZE) {
			count_a == 0;
		}
		__syncthreads();
		if (commit) {
			count_a_t = atomicAdd((int *)&count, 1);
			if (count_a_t < MAX_BLOCK_SIZE) {
				spike_a[count_a_t].time = time + d->delay_;
				spike_a[count_a_t].n =  d;
				commit = false;
			}
		}
		__syncthreads();
		if (count_a > 0) {
			commitFireInfo(spike_a, MAX_BLOCK_SIZE, spike1, &s1_size[0], 0);
		}
		if (threadIdx.x == 0) {
			count_a == 0;
		}
		__syncthreads();

		if (commit_b) {
			count_b_t = atomicAdd((int *)&count_a, 1);
			if (count_b_t < MAX_BLOCK_SIZE) {
				spike_b[count_a_t].time = time + d->delay_;
				spike_b[count_a_t].n =  d;
				commit_b = false;
			} else {
				commitFireInfo(spike_b, MAX_BLOCK_SIZE, spike2, &s2_size[0], 0);
			}
		}
		__syncthreads();

		if (threadIdx.x == 0 && count >= MAX_BLOCK_SIZE) {
			count_b == 0;
		}
		__syncthreads();
		if (commit) {
			count_b_t = atomicAdd((int *)&count_b, 1);
			if (count_b_t < MAX_BLOCK_SIZE) {
				spike_b[count_a_t].time = time + d->delay_;
				spike_b[count_a_t].n =  d;
				commit = false;
			}
		}
		__syncthreads();
		if (count_b > 0) {
			commitFireInfo(spike_b, MAX_BLOCK_SIZE, spike2, &s2_size[0], 0);
		}
		if (threadIdx.x == 0) {
			count_b == 0;
		}
		__syncthreads();

	}
}


__global__  reset_size(int *s_size, FireInfo *spike1, int s1_size, FireInfo *spike2, int *s2_size, double time) {
