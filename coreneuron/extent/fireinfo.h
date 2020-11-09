
#ifndef FIREINFO_H
#define FIREINFO_H

class FireInfo {
public:
	NetCon *n;
	double time;
};

__device__ int commitFireInfo(FireInfo *shared_buf, volatile unsigned int size, FireInfo *global_buf, int *global_size, int offset);

#endif // FIREINFO_H
