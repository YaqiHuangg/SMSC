#ifndef __CORE_GRAPH_H
#define __CORE_GRAPH_H

#include <string>
#include <map>
#include <vector>
#include <limits>
#include <iostream>

using namespace std;

typedef unsigned short vertex_t;

extern map<string, int> vlabel_map;
extern map<string, int> elabel_map;

extern vector<string> vertex_label;
extern vector<string> edge_label;

#define EPSILON numeric_limits<unsigned short>::max()

class coregraph
{
private:
	int* vlabels;  // vertex labels
	int** elabels; // adjacent matrix

	unsigned num_vertices;
	unsigned num_edges;

private:
    unsigned* vfreq;
    unsigned* efreq;

    unsigned** vertex_efreq;

private:
	vector<unsigned>* vlabel_index;//保存节点标签出现的位置，即序号

private:
	/*leaves:
      第1维下标为每个顶点id，第2维保存与第1维相应顶点相连的叶子顶点的相关信息,
      第2维里的vector保存与第1维相应顶点相连的叶子顶点id，这些叶子顶点标签相同且与u相连的边标签也相同
      按照每个顶点与其相连不同标签的叶子顶点的个数升序排列
      保存的顺序先是边标签再是顶点标签*/
	vector<vector<pair<pair<unsigned int, unsigned int>, vector<int>>>> leaves;//leaves -> vector(pair(edge label, vertex label), vector(leaves list))

public:
	coregraph(unsigned n);
	~coregraph();

public:
	unsigned vsize(){ return num_vertices; }
	unsigned esize(){ return num_edges; }

public:
    int vlabel(int v){
        if(v == EPSILON) return EPSILON;
        return vlabels[v];
    }

    int elabel(int v1, int v2){
        if(v1 == EPSILON || v2 == EPSILON) return 0;
        return elabels[v1][v2];
    }

public: // accessor to vertex & edge labels
	int operator()(int v){ 
		if(v == EPSILON) return EPSILON;
		return vlabels[v];
	}
	int operator()(int v1, int v2){
		if(v1 == EPSILON || v2 == EPSILON) return 0;
		return elabels[v1][v2];
	}

public:
	void setVertexLabel(unsigned v, string label);
	void setEdgeLabel(unsigned v1, unsigned v2, string label);
	void setEdgeBond(unsigned v1, unsigned v2, int bond);
	void setLeaves(vector<vector<pair<pair<unsigned int, unsigned int>, vector<int>>>> leaf);

public:
    void labelScan();
	void leaveScan();
	void labelRescan(unsigned sz_v, unsigned sz_e);
    
    unsigned* getVertexFrequencies(){ return vfreq; }
    unsigned* getEdgeFrequencies(){ return efreq; }
    unsigned** getVertexEdgeFrequencies(){ return vertex_efreq; }
	vector<unsigned>* getVertexLabelIndex(){ return vlabel_index; }
	vector<vector<pair<pair<unsigned int, unsigned int>, vector<int>>>> getLeaves(){ return leaves;}

	vector<short>* adjacent;//保存每个顶点与其相邻的顶点

public:
	unsigned sizeFilter(coregraph* rh);
	unsigned labelFilter(coregraph* rh);

public:
	void output(ostream& out, int number = 0);
};

#endif // __CORE_GRAPH_H
