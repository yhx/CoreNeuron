
#include "connceciton.h"

Connection::Connection(int n, int *start, int *end) 
{
	_n = n;
	_start = new int[_n];
	_num = new int[_n];
}

Connection::~Connection()
{
	delete[] _start;
	delete[] _num;
	_n = 0;
}

void Connection::save(char *name)
{
	FILE *f = fopen(name, "wb+");
	fwrite(&_n, sizeof(int), 1, f);
	fwrite(&_start, sizeof(int), _n, f);
	fwrite(&_num, sizeof(int), _n, f);
	fclose(f);
}

void Connection:load(char *name)
{
	if (_n > 0) {
		delete[] _start;
		delete[] _end;
	}

	FILE *f = fopen(name, "wb+");
	fread(&_n, sizeof(int), 1, f);
	_start = new int[_n];
	_num = new int[_n];
	fread(&_start, sizeof(int), _n, f);
	fread(&_num, sizeof(int), _n, f);
	fclose(f);
}
