#ifndef __NASS_H
#define __NASS_H

#include <algorithm>
#include <cstdlib>
#include <list>
#include <mutex>

#include "dataset.h"
#include "coregraph.h"
#include "looptimer.h"
#include "workload.h"

#define PARTIME 1 // incremental partitioning time
#define GEDTIME 2 // ged computation time

extern int NTHREADS;
extern unsigned MEMLIMIT;

class nass
{
public:
	static int threshold;

private:
	loop_timer* timer;
	vector<int> res_vec;
	vector<int>* tgData;//保存各个类的图id（第一维为类id，第二维为该类的图id）
	int* dataTG;//保存每个图所属的类id

	void clearStats();

private:
	void basicFilter(coregraph* x, vector<pair<int, int> >& candidates, int threshold, int groupid, unsigned sum, int begin = 0);
	void basicFilter(coregraph* x, vector<pair<int, int> >& candidates, int threshold, int begin = 0);
	int basicFilter(coregraph* x, coregraph* y);
	int mcsbasicFilter(coregraph* x, coregraph* y);
	unsigned regenCandidates(vector<pair<int, int> >& candidates, unsigned pos, int res, int distance);
	// unsigned regenCCandidates(vector<pair<int, int> >& candidates, vector<int> c_ccandiddates, int g_pos, int c_pos, int gid, int cid, int distance);
	void readTGFiles();
	void readMCSIndexFile();
	void search(coregraph* x);

public:
	nass(const char* filename, const char* thresholdGroupFile, const char* MCSIndexFile, bool load_index = true);
	~nass(); //{ delete [] nassIndex; }

public: 
	void run(Workload& workload);

/******************************************************/
/*                    INDEXING                        */
/******************************************************/
private:
	const char* indexfile;
	const char* MCSindexfile;
	const char* thresholdGroupFile;
	int index_max_threshold;
	int mcs_max_threshold;
	vector< int >* index;////nassIndex没有distance版本
	vector< pair<int, int> >* nassIndex;
	vector< pair<int, int> >* MCSIndex;
	vector<vector< int >> thresholdGroup;//第一维下标表示当前图id，第二维保存与其ged值在5之内的图id（包括自己）
	vector<vector<int>> newThresholdGroup;//保存精简后的阈值数组

	bool use_index;
	bool* sampled;
	bool use_sample;

public: 
	void buildNassIndex(bool distributed, char* address, int sampling_rate);
	bool loadNassIndex(bool verbose = true);
	unsigned getThresholdGroupSize(){return newThresholdGroup.size();}

// for measuring indexing time
private:
    TIMESTAMP ts1, ts2;

// for multi-threaded index building
private:
	mutex id_mutex;
	mutex out_mutex;

	int  num_index_worker;
	bool done_indexing;

	unsigned next_graph_id;
	int next_worker_id;

	int getNextGraphIDLocal();
	int getNextGraphID(vector<pair<int, int> >* res = NULL, int id = -1);

	int getWorkerID();
	void outputResult(vector<pair<int, int> >& res, int gid);
	void outputThresholdGroup(vector<vector<int>>& newThresholdGroup);
	void searchWorker();

	void monitor(); // status mointoring thread

// for distributed index building
private:
	char* server_address;
	int sv_sock; // server socket
	bool distributed; // distributed index building if true
	bool remote; // client if remote == true

	int connectToCoordinator(const char* address = NULL);
	int initCoordinator(struct sockaddr_in&, struct sockaddr_in&);

	int getNextGraphIDRemote(vector<pair<int, int> >* res, int id);

	void coordinateLoop();

	void resort(vector<vector< int >>& thresholdGroup);
	void recluster(int tid);
	void refine();
        void miniRefine(vector<vector< int >>& group);
};

bool gtid(const pair<int, int>& c1, const pair<int, int>& c2);
bool gtlb(const pair<int, int>& c1, const pair<int, int>& c2);
bool mSize(const pair<int, int> &x, const pair<int, int> &y);
bool groupLength(vector<int> x, vector<int> y);
bool numSize(int x, int y);

#endif //__NASS_H
