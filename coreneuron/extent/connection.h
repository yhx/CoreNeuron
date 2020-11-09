
#ifndef CONNECTION_H
#define CONNECTION_H

class Connection {
public:
	int *_start;
	int *_num;
	int _n;

	Connection(int n, int *start, int *end);
	~Connection();
protected:
	void save(char *name);
	void load(char *name);
};

typedef Connection PreSyn;

#endif // SPIKEINFO_H
