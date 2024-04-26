#include "mcs.h"
#include "graph.h"

#include <vector>
#include <sstream>
#include <iostream>
#include <set>
#include <atomic>
#include <limits.h>
#include <string.h>
#include <numeric>
#include <thread>
#include <exception>
#define VSCORE

using namespace std;

// 原子操作
static atomic<bool> abort_due_to_timeout;

int dl;
unsigned long long nodes;//赋初始值为0，表示公共子图中的顶点数
unsigned long long bestnodes,bestcount;
unsigned long long cutbranches;
unsigned long long conflicts;
int timeout = 3600;
clock_t start;

mcs::mcs(const char* thresholdGroupFile, const char* mcsFile, const char* heuristic, unsigned mcsDataSize)
{
    this->thresholdGroup_file = thresholdGroupFile;
    this->mcs_file = mcsFile;
    this->heuristic = heuristic;
    // this->mcsData_size = mcsDataSize;
    this->mcsData = NULL;

    setMcsDataSize(mcsDataSize);
}

mcs::~mcs()
{
    // delete [] &mcsData;
    delete [] mcsData;
}

void mcs::setMcsDataSize(unsigned mcsDataSize)
{
    mcsData_size = mcsDataSize;
}

vector<unsigned> mcs::calculate_degrees(coregraph* g)
{
    vector<unsigned> degree(g->vsize(), 0);//初始化g中每个顶点度数为0
    unsigned** vertex_efreq = g->getVertexEdgeFrequencies();
    for (unsigned v=0; v<g->vsize(); v++) {
        degree[v] = vertex_efreq[v][elabel_map.size()];
    }
    return degree;
}

coregraph* mcs::induced_subgraph(coregraph* g, vector<int> vv) 
{
    coregraph* g_sorted = new coregraph(g->vsize());

    unsigned vertexSize = g->vsize();
	// unsigned edgeSize = g->esize();

    for (unsigned i=0; i<vertexSize; i++)
        g_sorted->setVertexLabel(i, to_string(g->vlabel(vv[i])));

    for (unsigned i=0; i<vertexSize; i++){
        for (unsigned j=0; j<vertexSize; j++){
            if(g->elabel(vv[i], vv[j]))
                g_sorted->setEdgeLabel(i, j, to_string(g->elabel(vv[i], vv[j])));//vv向量标志着顶点度数顺序
        }
    }
    g_sorted->labelScan();
    // g_sorted->leaveScan();

    return g_sorted;
}

int mcs::calc_bound(const vector<Bidomain>& domains)
{
    int bound = 0;
    for (const Bidomain &bd : domains) {
        bound += min(bd.left_len, bd.right_len);
    }
    return bound;
}

int mcs::selectV_index(const vector<int>& arr,const vector<double> & lgrade, int start_idx, int len)
{
    int idx = -1;//保存arr中从start_idx开始最大分数的顶点下标i，实际下标为start_idx+i
    double max_g = -1;
    int vtx, best_vtx = INT_MAX;
    for(int i = 0; i < len; i++) {
        vtx = arr[start_idx + i];
        if(lgrade[vtx] > max_g) {
            idx = i;
            best_vtx = vtx;
            max_g = lgrade[vtx];
        }
        else if(lgrade[vtx] == max_g) {
            if(vtx < best_vtx) {
                idx = i;
                best_vtx = vtx;
            }
        }
    }
    return idx;
}

int mcs::selectW_index(const vector<int> &arr, const vector<double> &rgrade, int start_idx, int len, const vector<int> &wselected) 
{
    int idx = -1;
    double max_g = -1;
    int vtx, best_vtx = INT_MAX;
    for(int i = 0; i < len; i++) {
        vtx = arr[start_idx + i];
        if(wselected[vtx] == 0) {
            if(rgrade[vtx] > max_g) {
                idx = i;
                best_vtx = vtx;
                max_g = rgrade[vtx];
            }
            else if(rgrade[vtx] == max_g) {
                if(vtx < best_vtx) {
                    idx = i;
                    best_vtx = vtx;
                }
            }
        }
    }
    return idx;//返回分数和度数最大的顶点id
}

int mcs::select_bidomain(const vector<Bidomain>& domains, const vector<int> & left,const vector<double> & lgrade,
        unsigned current_matching_size)
{
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int tie_breaker;
    unsigned int i;  int len;
    int best = -1;
    for (i=0; i<domains.size(); i++) {
        const Bidomain &bd = domains[i];
        if (current_matching_size>0 && !bd.is_adjacent) continue;
            // 先判断arguments.heuristic == min_max ？最后给len赋值
            len = string(heuristic) == "min_max" ?
                max(bd.left_len, bd.right_len) :
                bd.left_len * bd.right_len;//最多有这么多种连线
        if (len < min_size) {//优先计算可能连线少的情况(即标签频率最低的domain),当domain中可能配对数最小时，选择domain中含有最大度的点那个domain
            min_size = len;//此时进行修改将度最大的点进行优化
            min_tie_breaker = left[bd.l + selectV_index(left, lgrade, bd.l, bd.left_len)];//序号越小度越大，选择分数最高的点，min_tie_breaker保存顶点id
            best = i;
        } else if (len == min_size) {
            tie_breaker = left[bd.l + selectV_index(left,lgrade,bd.l, bd.left_len)];
            if (tie_breaker < min_tie_breaker) {//选择序号小的breaker, 因为度数更大
                min_tie_breaker = tie_breaker;
                best = i;
            }
        }
    }
    return best;//返回标签频率最小的domain下标
}

void mcs::remove_vtx_from_array(vector<int> &arr, int start_idx, int &len, int remove_idx)
{
    len--;
    /*arr中的值已经代表了排序，所以这样交换去除点，不会影响后面的排序
    顶点已按照度数重新分配id，arr中存储的值就是顶点id*/
    swap(arr[start_idx + remove_idx], arr[start_idx + len]);
}

int mcs::partition(vector<int>& all_vv, int start, int len, coregraph* g, int v) 
{
    int i=0;
    
    for (int j=0; j<len; j++) {
        if(g->elabel(v,all_vv[start+j])){
            swap(all_vv[start+i], all_vv[start+j]);
            i++;
        }
    }

    return i;
}

coregraph* mcs::transToGraph(vector<VtxPair>& solution, coregraph* g)
{
    // solution vertices in ascending order of id (descending order of degree)
    coregraph* mcs = new coregraph(solution.size());

    // for (auto& vtx_pair : solution) {
    //     vtx_pair.v = vv0[vtx_pair.v];
    //     vtx_pair.w = vv1[vtx_pair.w];
    // }

    for(unsigned i = 0; i < solution.size(); i++){
        int v = solution[i].v;//vertex id
        mcs->setVertexLabel(i, to_string(g->vlabel(v)));
    }

    for(unsigned i = 0; i < solution.size(); i++){
        int v = solution[i].v;//vertex id
        for(unsigned j=i;j<solution.size();j++){
            int w = solution[j].v;
            if(g->elabel(v,w))
                mcs->setEdgeLabel(i,j,to_string(g->elabel(v,w)));
        }
    }
    mcs->labelScan();
    peak_leaves(mcs);

    return mcs;
}

int mcs::getNextGroupID()
{
    tid_mutex.lock();
	int tid = -1;

	// if(next_group_id == mcsData.size()){
	// 	tid_mutex.unlock();
	// 	return -1;
	// }
    if(next_group_id == mcsData_size){
		tid_mutex.unlock();
		return -1;
	}

	tid = next_group_id++;

	tid_mutex.unlock();
	return tid;
}

void mcs::OutputMCS(coregraph* mcs, unsigned tid)
{
    // cout<<"successn enter OutputMCS, please wait..."<<endl;
    static unsigned num_output = 0;
	static unsigned progress = 0;

    out_mutex.lock();

    ofstream out(mcs_file, ios::app);

    out<<"mcs # "<<tid<<endl;//threshold group id
    // int length = mcsGroup.size();
    // coregraph* g = mcsData[length-1];//the mcs of the id-th threshold group
    for(unsigned i=0;i<mcs->vsize();i++)
        out<<"v "<<i<<' '<<mcs->vlabel(i)<<endl;
    for(unsigned i=0;i<mcs->vsize();i++){
        for(unsigned j=i+1;j<mcs->vsize();j++){
            if(mcs->elabel(i,j))
                out<<"e "<<i<<' '<<j<<' '<<mcs->elabel(i,j)<<endl;
        }
    }

    num_output++;
    unsigned percent = num_output*100/mcsData_size;
    if(percent > progress){
        cout<<"\rExtracting mcs: " << percent << "%";
        cout << flush;
        progress = percent;
    }

    out << flush;
    out.close();


    // for(unsigned i=0;i<mcsData_size;i++){
    //     out<<"mcs # "<<i<<endl;//threshold group id
    //     coregraph* g = mcsData[i][mcsData[i].size()-1];//the mcs of the id-th threshold group
    //     for(unsigned j=0;j<g->vsize();j++){
    //         out<<"v "<<j<<' '<<g->vlabel(j)<<endl;
    //     }
    //     for(unsigned j=0;j<g->vsize();j++){
    //         for(unsigned k=j+1;k<g->vsize();k++){
    //             if(g->elabel(j,k))
    //                 out<<"e "<<j<<' '<<k<<' '<<g->elabel(j,k)<<endl;
    //         }
    //     }
    //     out << flush;

    //     num_output++;
    //     unsigned percent = num_output*100/mcsData_size;
    //     if(percent > progress){
    //         cout<<"\rOutputting mcs: " << percent << "%";
    //         cout << flush;
    //         progress = percent;
    //     }
    // }

    // out << endl << flush;
	// out.close();

    out_mutex.unlock();
}

int mcs::remove_matched_vertex(vector<int> &arr, int start, int len, const vector<int> &matched) {
    int p = 0;
    for(int i = 0; i < len; i++) {
        if(!matched[arr[start + i]]) {//该顶点未匹配
            swap(arr[start + i], arr[start + p]);//将未匹配的顶点换到数组前面
            p++;
        }
    }
    return p;//返回未匹配的顶点个数
}

vector<Bidomain> mcs::rewardfeed(const vector<Bidomain> & d,int bd_idx,vector<VtxPair> & current, vector<int> &g0_matched, vector<int> &g1_matched,
        vector<int> & left, vector<int> & right,vector<double> & lgrade,vector<double> & rgrade,
        coregraph* g0, coregraph* g1, int v, int w,
        bool multiway)
{
    //每个domain均分成2个new_domain分别是与当前对应点相连或者不相连
    //left中存储的是g0中与g1的公共标签的顶点id
    //right中存储的是g1中与g0的公共标签的顶点id
    current.push_back(VtxPair(v, w));
    g0_matched[v] = 1;
    g1_matched[w] = 1;

    int leaves_match_size = 0, v_leaf, w_leaf;
    vector<vector<pair<pair<unsigned int, unsigned int>, vector<int>>>> leaves0 = g0->getLeaves();
    vector<vector<pair<pair<unsigned int, unsigned int>, vector<int>>>> leaves1 = g1->getLeaves();
    //匹配分别与v、w相连的叶子顶点
   for(unsigned int i = 0, j = 0; i < leaves0[v].size() && j < leaves1[w].size(); ) {
        // pair先按first比较，如果相等，再按second比较
        if(leaves0[v][i].first < leaves1[w][j].first) i++;
        else if(leaves0[v][i].first > leaves1[w][j].first) j++;
        else {
            const vector<int> &leaf0 = leaves0[v][i].second;
            const vector<int> &leaf1 = leaves1[w][j].second;
            for(unsigned int p = 0, q = 0; p < leaf0.size() && q < leaf1.size(); ) {
                if(g0_matched[leaf0[p]]) p++;
                else if(g1_matched[leaf1[q]]) q++;
                else {
                    v_leaf = leaf0[p], w_leaf = leaf1[q];//匹配新的叶顶点
                    p++, q++;
                    current.push_back(VtxPair(v_leaf, w_leaf));
                    g0_matched[v_leaf] = 1;
                    g1_matched[w_leaf] = 1;
                    leaves_match_size++;
                }
            }
            i++, j++;
        }
    }

    vector<Bidomain> new_d;
    new_d.reserve(d.size());//为new_d分配d.size()个内存,这样每个域都有两个new_domain.
    // unsigned int old_bound = current.size() + calc_bound(d)+1;//这里的domain是已经去掉选了的点，但是该点还没纳入current所以+1
    //unsigned int new_bound(0);
    int l,r,j=-1;
    int temp=0,total=0;
    int unmatched_left_len, unmatched_right_len;
    for (const Bidomain &old_bd : d) {
        j++;
        l = old_bd.l;
        r = old_bd.r;
        if(leaves_match_size > 0 && old_bd.is_adjacent == false) {//不相邻的域
            unmatched_left_len = remove_matched_vertex(left, l, old_bd.left_len, g0_matched);//返回未匹配的顶点个数，且该函数将未匹配的顶点换到了数组前面
            unmatched_right_len = remove_matched_vertex(right, r, old_bd.right_len, g1_matched);
        }
        else {
            unmatched_left_len = old_bd.left_len;
            unmatched_right_len = old_bd.right_len;
        }
        // After these two partitions, left_len and right_len are the lengths of the
        // arrays of vertices with edges from v or w (int the directed case, edges
        // either from or to v or w)
        int left_len = partition(left, l, unmatched_left_len, g0, v);//返回与v相连的顶点数, 并且将相连顶点依次交换到left数组前面, 度数大的在前
        int right_len = partition(right, r, unmatched_right_len, g1, w);//w同上
        int left_len_noedge = unmatched_left_len - left_len;//与v不相连的顶点数
        int right_len_noedge = unmatched_right_len - right_len;//w同上
//这里传递v,w选取的下标到bd_idx，j用来计算当前所分domain的下标，这样就能将bound更精确地计算
        temp=min(old_bd.left_len,old_bd.right_len)-min(left_len,right_len)-min(left_len_noedge,right_len_noedge);
        total+=temp;//最小下界

#ifdef DEBUG

        printf("adj=%d ,noadj=%d ,old=%d ,temp=%d \n",std::min(left_len,right_len),
               std::min(left_len_noedge,right_len_noedge),std::min(old_bd.left_len,old_bd.right_len),temp);
        cout<<"j="<<j<< "  idx="<<bd_idx<<endl;
        cout<<"gl="<<lgrade[v]<<" gr="<<rgrade[w]<<endl;
#endif
        if (left_len_noedge && right_len_noedge)//new_domain存在的条件是同一domain内需要存在与该对应点同时相连，或者同时不相连的点
            new_d.push_back({l+left_len, r+right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});//不相邻分区
        if (multiway && left_len && right_len) {//与顶点相连的点,且存在边标签
            auto& adjrow_v = g0->adjacent[v];//度数大的顶点排在前面
            auto& adjrow_w = g1->adjacent[w];
            auto l_begin = begin(left) + l;
            auto r_begin = begin(right) + r;
            sort(l_begin, l_begin+left_len, [&](int a, int b)
                    { return adjrow_v[a] < adjrow_v[b]; });//当有边标签时，按照标签升序排列
            sort(r_begin, r_begin+right_len, [&](int a, int b)
                    { return adjrow_w[a] < adjrow_w[b]; });
            int l_top = l + left_len;
            int r_top = r + right_len;
            while (l<l_top && r<r_top) {
                int left_label = adjrow_v[left[l]];//边标签
                int right_label = adjrow_w[right[r]];
                if (left_label < right_label) {
                    l++;
                } else if (left_label > right_label) {
                    r++;
                } else {
                    int lmin = l;
                    int rmin = r;
                    do { l++; } while (l<l_top && adjrow_v[left[l]]==left_label);
                    do { r++; } while (r<r_top && adjrow_w[right[r]]==left_label);
                    new_d.push_back({lmin, rmin, l-lmin, r-rmin, true});//相邻且边标签相同分区
                }
            }
        } else if (left_len && right_len) {//与顶点相连的点,不存在边标签
            new_d.push_back({l, r, left_len, right_len, true});//与顶点相连的顶点 标志为Ture
        }
    }
    if(total>0){
        conflicts++;
        lgrade[v]+=total;//打分
        if(lgrade[v] > short_memory_threshold){//(10)^9
            for(int i=0; i<lgrade.size(); i++)
                lgrade[i]=lgrade[i]/2;
        }

        rgrade[w]+=total;
        if(rgrade[w] > long_memory_threshold){
            for(int i=0; i<rgrade.size(); i++)
                rgrade[i]=rgrade[i]/2;
        }
    }
#ifdef DEBUG
  cout<<"new domains are "<<endl;
  for (const Bidomain &testd : new_d){
      l = testd.l;
      r = testd.r;
      for(j=0;j<testd.left_len;j++)
          cout<<left[l+j] <<" ";
      cout<<" ; ";
      for(j=0;j<testd.right_len;j++)
          cout<<right[r+j]<<" " ;
      cout<<endl;
  }
#endif
    return new_d;//两个分区，第一个是不相邻分区，第二个是相邻分区
}

void mcs::remove_bidomain(vector<Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size()-1];
    domains.pop_back();
}

void mcs::solve(coregraph* g0, coregraph* g1,vector<double> &V, vector<vector<double>> &Q, vector<VtxPair> & incumbent,
        vector<VtxPair> & current, vector<int> &g0_matched, vector<int> &g1_matched, vector<Bidomain> & domains,
        vector<int> & left, vector<int> & right, unsigned int matching_size_goal)
{
    if(timeout && double(clock() - start) / CLOCKS_PER_SEC > timeout)
    {
        // cout <<"time out" <<endl;
        return;
    }
    nodes++;

    if (current.size() > incumbent.size()) {//incumbent 现任的
        incumbent = current;
        bestcount = cutbranches+1;
        bestnodes = nodes;
        // bestfind = clock();
    }

    unsigned int bound = current.size() + calc_bound(domains);//计算相连和不相连同构数的最大可能加上当前已经同构的点数
    if (bound <=incumbent.size() || bound < matching_size_goal){//剪枝
        cutbranches++;
        return;
    }
    if (incumbent.size()==matching_size_goal)
       return;

    int bd_idx = select_bidomain(domains, left, V, current.size());//选出domain中可能情况最少(标签频率最低)的下标
    if (bd_idx == -1) {  // In the MCCS case, there may be nothing we can branch on
      //  cout<<endl;
      //  cout<<endl;
      //  cout<<"big"<<endl;
        return;}
    Bidomain &bd = domains[bd_idx];//一个domain代表一个标签的情况

    int v, w;
    int tmp_idx;

    tmp_idx = selectV_index(left, V, bd.l, bd.left_len);
    v = left[bd.l + tmp_idx];//选取v值时不再是选序号最小的点，而是选grade最大的点，当grade均相同时，则选度最大的。
    remove_vtx_from_array(left, bd.l, bd.left_len, tmp_idx);//从left中去除了度最大的点

    // Try assigning v to each vertex w in the colour class beginning at bd.r, in turn
    vector<int> wselected(g1->vsize(), 0);//标记g1中的顶点是否被选择
    bd.right_len--;//下面循环中用的<=
    for (int i=0; i<=bd.right_len; i++) {//选取与v匹配的w
        tmp_idx = selectW_index(right, Q[v], bd.r, bd.right_len + 1, wselected);//同selectV
        w = right[bd.r + tmp_idx];
        wselected[w]=1;
        swap(right[bd.r + tmp_idx], right[bd.r + bd.right_len]);
        right[bd.r + bd.right_len] = w;
#ifdef DEBUG
     cout<<"v= "<<v<<" w= "<<w<<endl;
     unsigned int m;
     for(m=0;m<left.size();m++)
         cout<<left[m]<<" ";
     cout<<endl;
     for(m=0;m<right.size();m++)
         cout<<right[m]<<" ";
     cout<<endl;
#endif
        unsigned int cur_len = current.size();
        auto new_domains = rewardfeed(domains,bd_idx,current, g0_matched, g1_matched, left, right, V, Q[v], g0, g1, v, w);//返回与v、w不相邻和相邻的分区
        dl++;
        solve(g0, g1, V, Q, incumbent, current, g0_matched, g1_matched, new_domains, left, right, matching_size_goal);
        while(current.size() > cur_len) {
            VtxPair pr = current.back();
            current.pop_back();
            g0_matched[pr.v] = 0;
            g1_matched[pr.w] = 0;
        }
    }
    bd.right_len++;
    if (bd.left_len == 0)
        remove_bidomain(domains, bd_idx);
    solve(g0, g1, V, Q, incumbent, current, g0_matched, g1_matched, domains, left, right, matching_size_goal);
}

vector<VtxPair> mcs::extractMCS(coregraph* g0, coregraph* g1)
{
    start=clock();
    
    //存储g0中与g1公共标签的顶点id
    vector<int> left;  // the buffer of vertex indices for the left partitions
    //存储g1中与g0公共标签的顶点id
    vector<int> right;  // the buffer of vertex indices for the right partitions

    vector<int> g0_matched(g0->vsize(), 0);//标记g0中的顶点是否已匹配，1表示已匹配
    vector<int> g1_matched(g1->vsize(), 0);//标记g1中的顶点是否已匹配，1表示已匹配

    vector<double> V(g0->vsize(), 0);//保存g0中每个顶点的分数
    vector<vector<double>> Q(g0->vsize(), vector<double> (g1->vsize(), 0));//保存在选择g0中v顶点的情况下g1中各顶点的分数

    //保存g0和g1各个公共标签的顶点起始位置和顶点数量
    auto domains = vector<Bidomain> {};

    // set：不会出现重复的内容
    std::set<int> left_labels;
    std::set<int> right_labels;
    for(unsigned i=0;i<g0->vsize();i++)
        left_labels.insert(g0->vlabel(i));
    for(unsigned i=0;i<g1->vsize();i++)
        right_labels.insert(g1->vlabel(i));

    std::set<int> labels;  // labels that appear in both graphs
    //求交集，保存在labels中
    std::set_intersection(std::begin(left_labels),
                          std::end(left_labels),
                          std::begin(right_labels),
                          std::end(right_labels),
                          std::inserter(labels, std::begin(labels)));//它将元素插入到labels的std::begin(labels)位置中。
    
    // Create a bidomain for each label that appears in both graphs
    for (int label : labels) {
        int start_l = left.size();
        int start_r = right.size();

        for (unsigned int i=0; i<g0->vsize(); i++)
            if (g0->vlabel(i)==label)
                left.push_back(i);
        for (unsigned int i=0; i<g1->vsize(); i++)
            if (g1->vlabel(i)==label)
                right.push_back(i);

        int left_len = left.size() - start_l;
        int right_len = right.size() - start_r;
        domains.push_back({start_l, start_r, left_len, right_len, false});
    }

    vector<VtxPair> incumbent;

    for (unsigned k=0; k<g0->vsize(); k++) {
        unsigned int goal = g0->vsize() - k;//公共子图的目标尺寸
        auto left_copy = left;
        auto right_copy = right;
        auto domains_copy = domains;
        vector<VtxPair> current;//保存当前已同构的顶点数
        conflicts=0;
        solve(g0, g1, V, Q, incumbent, current, g0_matched, g1_matched, domains_copy, left_copy, right_copy, goal);//匹配顶点，划分分区
        if (incumbent.size() == goal || abort_due_to_timeout) break;
    }

    return incumbent;//保存匹配的顶点对
}

unsigned mcs::sum(const vector<unsigned> & vec) {
    /*在一个范围内累积值。使用运算符+()累计[begin(vec), end(vec)]范围内的值。初始值为0。这些值按顺序处理。*/
    return accumulate(std::begin(vec), std::end(vec), 0);
}

void mcs::peak_leaves(coregraph* g)
{
    vector<int> deg(g->vsize(), 0);//下标为顶点id，保存每个顶点的度数

    unsigned vertexSize = g->vsize();

    g->leaveScan();
    vector<vector<pair<pair<unsigned int, unsigned int>, vector<int>>>> leaves = g->getLeaves();

    for (unsigned i=0; i<vertexSize; i++){
        for (unsigned j=0; j<vertexSize; j++){
            if(i!=j && g->elabel(i, j))
                deg[i]++;
        }
    }

    for(int u = 0; u < vertexSize; u++){
        for(int v = 0; v < vertexSize; v++)
            if(g->elabel(u, v) && u != v && deg[v] == 1){//v为叶子顶点
                pair<unsigned int, unsigned int> labels(g->elabel(u, v), g->vlabel(v));//保存边标签及叶子顶点标签
                int pos = -1;
                for(int k = 0; ; k++){
                    if(k == int(leaves[u].size())){//g.leaves[u].size()表示与u相连的叶子顶点的个数
                        leaves[u].push_back(make_pair(labels, vector<int> ()));//vector<int>()表示vector<int>的原始对象的写法，没有名称
                    }
                    if(leaves[u][k].first == labels){
                        pos = k;
                        break;
                    }
                }
                leaves[u][pos].second.push_back(v);//保存与u相连的叶子顶点id，这些叶子顶点标签相同且与u相连的边标签也相同
            }
        sort(leaves[u].begin(), leaves[u].end());//按照每个顶点与其相连不同标签的叶子顶点的个数升序排列
    }

    g->setLeaves(leaves);
}

void mcs::readThresholdGroupFile()
{
    if(mcsData) delete [] mcsData;
    
    DataSet* dataset = DataSet::getInstance();
    // mcsData.reserve(dataset->numGraphs());
    mcsData = new vector<coregraph*>[mcsData_size];
    
    ifstream in(thresholdGroup_file);

    if(!in.is_open()){
		cerr << "Error in opening thresholdGroup file" << endl;
		return;
	}

    string line;

    while(getline(in, line)){
        stringstream ss(line);
        int tid,gid;//tid为阈值组id
        // vector<coregraph*> graphs;//save graphs of the threshold group
        // vector<vector<unsigned>> degree;//保存当前组所有图各个顶点的度数
        // vector<vector<int>> grade;//保存每个顶点的得分
        // vector<vector<int>> vv;/*当图为全连接且有自环时，按照顶点度数升序排列并存储在vv中，否则以降序排列存储顶点id*/
        // vector<coregraph*> g_sorted;//将当前组中的每个图按照度数升序/降序分配id

        ss >> tid;
        // unsigned n=0;
        // mcsData.push_back(vector<coregraph*>());
        while(ss >> gid){
            coregraph* g = dataset->graphAt(gid);
            // graphs.push_back(g);

            vector<unsigned> g_deg = calculate_degrees(g);

            vector<int> vv0(g->vsize());
            iota(begin(vv0), end(vv0), 0);
            bool g_dense = sum(g_deg) > (g->vsize())*((g->vsize())-1);
            stable_sort(begin(vv0),end(vv0), [&](int a, int b) {
                return g_dense ? (g_deg[a]<g_deg[b]) : (g_deg[a]>g_deg[b]);
            });
            // vv.insert(vv.begin()+n,vv0);

            coregraph* g_sorted = induced_subgraph(g, vv0);
            peak_leaves(g_sorted);
            // g_sorted.push_back(g0_sorted);
            // g_sorted.insert(g_sorted.begin()+n, g0_sorted);
            mcsData[tid].push_back(g_sorted);
            // n++;
        }

        // mcsData.insert(mcsData.begin()+tid,g_sorted);
    }
    in.close();
}

void mcs::computeMCS()
{
    // static int progress = 0;

    int tid = getNextGroupID();

    while(tid != -1){
        // int percent = tid*100/(mcsData_size-1);
        // if(percent > progress){
        //     cout<<"\rExtracting mcs: " << percent << "%";
        //     cout << flush;
        //     progress = percent;
        // }

        nodes = 0;bestnodes = 0;bestcount = 0;cutbranches = 0;dl=0;
        if(mcsData[tid].size()>1){
            // vector<int> lgrade(mcsData[tid][0]->vsize(),0);//保存每个顶点的得分
            // vector<int> rgrade(mcsData[tid][1]->vsize(),0);

            vector<VtxPair> solution = extractMCS(mcsData[tid][0], mcsData[tid][1]);//同构子图
            coregraph* mcs = transToGraph(solution, mcsData[tid][0]);
            
            // solution0.clear();grade0.clear();grade1.clear();
            for(unsigned i=2;i<mcsData[tid].size();i++){
                // coregraph* g2_sorted = mcs0;
                // coregraph* g_sorted = mcsData[tid][i];
                vector<int> lgrade(mcs->vsize(),0);
                vector<int> rgrade(mcsData[tid][i]->vsize(),0);
                
                vector<VtxPair> solution = extractMCS(mcs, mcsData[tid][i]);
                mcs = transToGraph(solution, mcs);
            }
            mcsData[tid].push_back(mcs);
            // if(tid == 36717){
            //     cout<<"36717th mcs vertex label:"<<endl;
            //     for(unsigned i=0;i<mcs->vsize();i++)
            //         cout<<i<<' '<<mcs->vlabel(i)<<endl;
            // }
        }
        OutputMCS(mcsData[tid][mcsData[tid].size()-1], tid);
        // if(i==0)
        //     cout<<"success computeMCS"<<endl;
        // cout<<"success compute the mcs of the "<<tid<<"th group"<<endl;
        tid = getNextGroupID();
    }
    // cout<<"success computeMCS"<<endl;
}

void mcs::buildMCSGroup()
{
    GETTIME(&ts_1);

    cout<<"Reading mcs from "<<thresholdGroup_file<<endl<<flush;
    readThresholdGroupFile();
    cout<<"success input mcsData"<<endl;

    next_group_id = 0;
    cout<<"\rExtracting mcs: 0%" << flush;
    // computeMCS();
    std::thread* worker[NTHREADS];
    for(int i = 0; i < NTHREADS; i++)
		worker[i] = new thread(&mcs::computeMCS, this);

    for(int i = 0; i < NTHREADS; i++){
		worker[i]->join();
		delete worker[i];
	}

    cout<<endl;
    cout<<"success compute MCS"<<endl;

    GETTIME(&ts_2);

    cout << "Computing MCS Time: " << ts_2.tv_sec - ts_1.tv_sec << " seconds" << endl;
}