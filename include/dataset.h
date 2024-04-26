#ifndef __DATASET_H
#define __DATASET_H

#include "coregraph.h"
#include "datafile.h"

#include<algorithm>
#include <vector>
#include <map>
#include <iostream>
#include <cstdlib>

class DataSet
{
private:
	vector<coregraph*> data;
	vector<coregraph*> mcsdata;
    vector<pair<unsigned, unsigned>> datagroup;//第一个int表示该组图的节点个数（所有图节点个数相同），第二个int表示该组有多少个图
	unsigned max_vertices;

private:
	static DataSet* instance;

private:
	DataSet(){}
	~DataSet()
	{
		// delete data graphs
		for(unsigned i = 0; i < data.size(); i++) delete data[i];
		for(unsigned i = 0; i < mcsdata.size(); i++) delete mcsdata[i];
	}

// private:
// 	static bool genum(vector<pair<int, int>>* x, vector<pair<int, int>>* y);

// public:
// 	vector<vector<pair<int, int>>> resortThreasholdGroup;//第一维下标表示当前图id，第二维中pair的第一个int表示图id，第二个int表示distances(进一步细化阈值分组)

public:
	static DataSet* getInstance()
	{
		if(instance == NULL) instance = new DataSet();
		return instance;
	}

	static void finishUp() { delete instance; }

	void buildDataSet(DataFile& file, bool verbose = false);
	void buildMcsSet(DataFile& file);

	vector<coregraph*> getData() { return data; }
	vector<coregraph*> getMCSData() { return mcsdata; }
	// vector<vector<pair<int, int>>> getGroup() {return resortThreasholdGroup;}
	vector<pair<unsigned, unsigned>> getDatagroup() { return datagroup; }
	// void resort(vector<vector<pair<int, int>>>& resortThreasholdGroup);
	// bool groupLength(vector<pair<int, int>>* x, vector<pair<int, int>>* y);
	void groupData(vector<coregraph*>& data);
	unsigned numGraphs() { return data.size(); }
	unsigned numMCSs() { return mcsdata.size(); }
    unsigned vmax() { return max_vertices; }
	coregraph* graphAt(int i){ return data[i]; }
	coregraph* mcsAt(int i){ return mcsdata[i]; }
};

#endif
