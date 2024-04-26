#ifndef __MCS_H
#define __MCS_H

#include "nass.h"
#include "dataset.h"
#include "coregraph.h"

#include <mutex>
const int short_memory_threshold = 1e5;
const int long_memory_threshold = 1e9;

struct VtxPair {//save the mapping vertex pairs
    int v;
    int w;
    VtxPair(int v, int w): v(v), w(w) {}
};

struct Bidomain {
    int l, r;        // start indices of left and right sets
    int left_len, right_len;
    bool is_adjacent;
    Bidomain(int l, int r, int left_len, int right_len, bool is_adjacent):
            l(l),
            r(r),
            left_len (left_len),
            right_len (right_len),
            is_adjacent (is_adjacent) { };
};

class mcs
{
private:
    const char* thresholdGroup_file;//保存阈值分组文件名
    const char* mcs_file;
    const char* resort_graph_file;
    // const char* class_file;
    const char* mcs_original_file;
    const char* heuristic;
    unsigned mcsData_size;
    vector<int>* mcs_original_id;//保存各类中各个图对应公共子图匹配的顶点id
    vector<int>* mcs_original_id_t;//保存各类图组中各个图对应公共子图在mcs_original_id中的下标
    /*除最后一项，偶数项（包括0）保存提取出的mcs的顶点id，奇数项保存前一项（即偶数项）作为匹配图匹配的顶点id
    最后一项只保存最终该类图提取出的mcs的顶点id
    */
    vector<int>* mcs_id;//保存各类图提取公共子图过程中产生的所有mcs的顶点id
    vector<int>* mcs_id_t;//保存各类图组对应产生的所有mcs在mcs_id中的下标
    vector<int>* tg;
    vector<coregraph*>* mcsData;//前面保存的是顶点经过排序后的图，最后一项保存的是mcs
    vector<coregraph*> resortData;
    // vector<VtxPair>* vp;

    // const char* mcs_id_file;
    // const char* mcs_original_id_file;

private:
    TIMESTAMP ts_1, ts_2;

private:
    mutex tid_mutex;
    mutex out_mutex;
    
    unsigned next_group_id;
    unsigned next_mid;// mid为各类中产生的mcs的id
    unsigned next_moid;// moid为该类中各图对应公共子图的id
    unsigned next_process_num;

public:
    mcs(const char* thresholdGroupFile, const char* resort_graph_file, const char* mcs_original_file, const char* mcsFile, const char* heuristic, unsigned mcsDataSize);
    ~mcs();

/************************************************************************/
/******************************compute***********************************/
/************************************************************************/
public:
    void buildMCSGroup();
    void setMcsDataSize(unsigned mcsDataSize);
    // vector<vector<coregraph*>> getMCSData(){ return mcsData;}
    vector<coregraph*>* getMCSData(){ return mcsData;}
    coregraph* mcsAt(int i, int j){ return mcsData[i][j]; }

private:
	int getNextGroupID();
    unsigned getNextMID();
    unsigned getNextMOID();

    unsigned getProcess();

private:
    unsigned sum(const vector<unsigned> & vec);
    vector<unsigned> calculate_degrees(coregraph* g);//计算每个图中各个顶点的度数
    coregraph* induced_subgraph(coregraph* g, vector<int> vv);//将顶点重新标号，按照度数升序/降序分配id
    // void peak_leaves(coregraph* g);
    //提取一对图的最大公共子图
    vector<VtxPair> extractMCS(coregraph* g0, coregraph* g1);
    void solve(coregraph* g0, coregraph* g1,vector<double> &V, vector<vector<double>> &Q, vector<VtxPair> & incumbent,
        vector<VtxPair> & current, vector<int> &g0_matched, vector<int> &g1_matched, vector<Bidomain> & domains,
        vector<int> & left, vector<int> & right, unsigned int matching_size_goal);
    int calc_bound(const vector<Bidomain>& domains);
    int select_bidomain(const vector<Bidomain>& domains, const vector<int> & left,const vector<double> & lgrade,
        unsigned current_matching_size);
    int selectV_index(const vector<int>& arr,const vector<double> & lgrade, int start_idx, int len);
    int selectW_index(const vector<int> &arr, const vector<double> &rgrade, int start_idx, int len, const vector<int> &wselected);
    void remove_vtx_from_array(vector<int> &arr, int start_idx, int &len, int remove_idx);
    void remove_bidomain(vector<Bidomain>& domains, int idx);
    int remove_matched_vertex(vector<int> &arr, int start, int len, const vector<int> &matched);
    int partition(vector<int>& all_vv, int start, int len, coregraph* g, int v);
    vector<Bidomain> rewardfeed(const vector<Bidomain> & d,int bd_idx,vector<VtxPair> & current, vector<int> &g0_matched, vector<int> &g1_matched,
        vector<int> & left, vector<int> & right,vector<double> & lgrade,vector<double> & rgrade,
        coregraph* g0, coregraph* g1, int v, int w,
        bool multiway = true);
     coregraph* transToGraph(vector<int>& vertex, coregraph* g);
    void OutputMCSOriginal(vector<int> mcs_original_id_t, unsigned tid);
    void OutputMCS(coregraph* mcs, unsigned tid);
    // void OutputMCS();
    void readThresholdGroupFile();
     void readMCSid();
    void readMCSOid();
    void refineMCSid();
    // void classifyGraph();
    void resortGraph();
    void computeMCS();

/************************************************************************/
/*******************************search***********************************/
/************************************************************************/
public:
    unsigned next_compute_num;

public:
    unsigned getCompute();

public:
    pair<coregraph*, vector<VtxPair>> buildMCS(coregraph* q, coregraph* g);
    // void buildMCS(coregraph* q, vector<coregraph*> mcsdata, vector<pair<int, pair<coregraph*, vector<VtxPair>>>> &qmcsdata);

};

#endif