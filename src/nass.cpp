#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits.h>

#include "nass.h"
#include "nassged.h"

using namespace std;

//statistics
long long cands_final;
long long push_count;

int nass::threshold = 0;

nass::nass(const char* filename, const char* thresholdGroupFile, const char* MCSIndexFile, bool load_index)
{
	this->indexfile = filename;
	this->MCSindexfile = MCSIndexFile;
	this->thresholdGroupFile = thresholdGroupFile;

	this->use_index = false;
	this->index = NULL;
	this->tgData = NULL;
	this->dataTG = NULL;
	this->nassIndex = NULL;
	this->MCSIndex = NULL;
	this->sampled = NULL;
	this->index_max_threshold = 0;

	if(load_index && indexfile && loadNassIndex())
		use_index = true;

	clearStats();
}

nass::~nass()
{
	delete [] index;
	delete [] tgData;
	delete [] nassIndex;
	delete [] MCSIndex;
	delete [] sampled;
	if(dataTG != NULL) delete [] dataTG;
}


void nass::clearStats()
{
	cands_final = 0;
	push_count = 0;
	res_vec.clear();
}

void nass::run(Workload& workload)
{
	readTGFiles();
	readMCSIndexFile();
	clearStats();

	timer = new loop_timer("Nass", "", workload.size(), NULL, 1); 
	timer->start();

	NassGED::initialize();
	for(unsigned i = 0; i < workload.size(); i++){
		//unsigned res_begin = res_vec.size(); // TO ACCESS iTH QUERY RESULTS

		// DataSet* dataset = DataSet::getInstance();
		// coregraph* g = dataset->graphAt(10490);
		// search(g);
		search(workload[i]);
		// cout<<"\nres = ";
		// for(int j = 0; j <res_vec.size(); j++)
		// 	cout<<res_vec[j]<<", ";
		// cout<<endl;

		/******** ACCESS QUERY RESULTS AS FOLLOWS *********
		cout << "Result for the " << i << "th query" << endl;
		DataSet* dataset = DataSet::getInstance();
		for(unsigned k = res_begin; k < res_vec.size(); k++){
			coregraph* g = dataset->graphAt(res_vec[k]);
			g->output(cout, k - res_begin);	
		}
		****************************************************/

		timer->next();
	}
	NassGED::finalize();

	timer->stop();
	long long elapsed = timer->end(true);
	cout << "===================== Results =====================" << endl;
	cout << "GED threshold: " << threshold << endl;
	cout << "# queries: " << workload.size() << endl;
	cout << "# total candidates: " << cands_final << endl;
	cout << "# total mappings: " << push_count << endl;
	cout << "# total results: " << res_vec.size() << endl;
	cout << "Processing time: ";
    fprintf(stdout, "%lld.%03lld seconds\n", elapsed/1000000, (elapsed % 1000000)/1000); 
	cout << "===================================================" << endl;

	delete timer;
}

bool gtid(const pair<int, int>& c1, const pair<int, int>& c2)
{
	return c1.first < c2.first;
}

bool gtlb(const pair<int, int>& c1, const pair<int, int>& c2)
{
	return c1.second < c2.second;
}

bool groupLength(vector<int> x, vector<int> y)
{
	return x.size()<y.size();
}

bool numSize(int x, int y)
{
	return x<y;
}

// unsigned nass::regenCCandidates(vector<pair<int, int> >& candidates, vector<int> c_ccandiddates, int g_pos, int c_pos, int gid, int cid, int distance)
// {

// }

unsigned nass::regenCandidates(vector<pair<int, int> >& candidates, unsigned pos, int res, int distance)
{
	// cout<<"hello\n";
	int max_threshold = index_max_threshold;
	if(!sampled[res]) max_threshold -= 1;

	res_vec.push_back(res);
	// cout<<"1: "<<res<<", ";
	if(distance + threshold > max_threshold) return pos;

	std::sort(candidates.begin() + pos, candidates.end(), gtid);

	vector<pair<int,int> > tmp;
	unsigned ptr1 = pos, ptr2 = 0;

	while(ptr1 < candidates.size() && ptr2 < nassIndex[res].size()){
		int gdist = nassIndex[res][ptr2].second;
		bool inexact = (gdist < 0);
		if(gdist == INEXACT_ZERO) gdist = 0;
		if(inexact) gdist = gdist * -1;

		if(gdist > threshold + distance) ptr2++;
		else if(nassIndex[res][ptr2].first < candidates[ptr1].first) ptr2++; 
		else if(nassIndex[res][ptr2].first > candidates[ptr1].first) ptr1++;
		else{
			if(gdist + distance <= threshold && !inexact){
				res_vec.push_back(nassIndex[res][ptr2].first);
				// cout<<"2: "<<nassIndex[res][ptr2].first<<", ";
			}
			else tmp.push_back(candidates[ptr1]);

			ptr1++; ptr2++;
		}
	}
	std::swap(candidates, tmp);

	std::sort(candidates.begin(), candidates.end(), gtlb);

	return 0;
}

void nass::basicFilter(coregraph* x, vector<pair<int, int> >& candidates, int threshold, int groupid, unsigned sum, int begin)
{
	DataSet* dataset = DataSet::getInstance();
	vector<coregraph*> data = dataset->getData();
	vector<pair<unsigned, unsigned>> datagroup = dataset->getDatagroup();
	ofstream outfile("data/candidates_group.txt",ios::app);
	// bool p=true;

	// for(unsigned i = begin; i < dataset->numGraphs(); i++){
	for(unsigned i = begin; i < data.size(); i++){
		//The initial sum means the position of in dataGroup[groupid]
		// coregraph* y = dataset->graphAt(i);
		// if(p){
		// 	outfile << endl;
		// 	outfile<<"NO."<<begin<<", candidates="<<endl;
		// 	outfile<<"{";
		// 	p=false;
		// }
		// outfile<<"i="<<i<<endl;

		coregraph* y = data[i];
		if(sum==datagroup[groupid].second){
			sum=1;//y enters another group
			groupid++;
		}
		else sum++;
		// outfile<<"y.vsize="<<y->vsize()<<",y.esize="<<y->esize()<<endl;
		int vdiff = (int)x->vsize() - (int)y->vsize();
		if(vdiff < 0) vdiff = -1*vdiff;
		if(vdiff > threshold){ 
			// outfile<<"break"<<endl;
			break;
		}
		// unsigned total_y=y->vsize()+y->esize();
		// if(begin<2) cout<<y->vsize()<<' ';
		if(x == y) continue; // skip the same graph
		if(x->sizeFilter(y) > unsigned(threshold)){
			// outfile<<"continue , sum="<<datagroup[groupid].second<<"and position="<<sum<<endl;
			i+=datagroup[groupid].second-sum;
			// if(i+1<data.size()){
			// 	unsigned m=i+1;
			// 	outfile<<"next:"<<m<<"y.vsize="<<data[m]->vsize()<<","<<m<<"y.esize="<<data[m]->esize()<<endl;
			// }
			groupid++;
			sum=0;
			continue;
		} 

		int lb = x->labelFilter(y); //gx.unmappedErrors(gy);
		if(lb <= threshold){
			candidates.push_back(pair<int, int>(i, lb));
			unsigned total = y->vsize()+y->esize();
			outfile << "id=" << i << "," << total << ' ';
		}
	}
	if(candidates.size()>0)
		outfile<<endl;
	outfile.close();

	vector<coregraph*>().swap(data);
	vector<pair<unsigned, unsigned>>().swap(datagroup);
}

void nass::basicFilter(coregraph* x, vector<pair<int, int> >& candidates, int threshold, int begin){
	DataSet* dataset = DataSet::getInstance();

	for(unsigned i = begin; i < dataset->numGraphs(); i++){
		coregraph* y = dataset->graphAt(i);
		if(x == y) continue; // skip the same graph
		if(x->sizeFilter(y) > unsigned(threshold)){
			continue;
		} 

		int lb = x->labelFilter(y); //gx.unmappedErrors(gy);
		if(lb <= threshold){
			candidates.push_back(pair<int, int>(i, lb));
		}
	}
}

int nass::basicFilter(coregraph* x, coregraph* y)
{
	if(x->sizeFilter(y) > unsigned(threshold))
		return x->sizeFilter(y);
	return x->labelFilter(y);
}

int nass::mcsbasicFilter(coregraph* x, coregraph* y)
{
	return x->labelFilter(y) - x->sizeFilter(y);
}

void nass::readTGFiles()
{
	DataSet* dataset = DataSet::getInstance();
	int* visited = new int[dataset->numGraphs()];
	memset(visited, 0, sizeof(int)*dataset->numGraphs());
	if(tgData) delete [] tgData;
	if(dataTG) delete [] dataTG;

    // cout << "\nReading threshold group file from " << tgfile << endl << flush;
	cout << "\nReading threshold group file from " << thresholdGroupFile << flush;

    DataSet* mcsset = DataSet::getInstance();
	tgData = new vector<int>[mcsset->numMCSs()];
	dataTG = new int[dataset->numGraphs()];

	ifstream in(thresholdGroupFile);

    if(!in.is_open()){
		cerr << "Error in opening threshold group file" << endl;
		exit(1);
	}

    string line;

    while(getline(in, line)){
        stringstream ss(line);
        int tid, gid;
        ss >> tid;//阈值组

        while(ss >> gid){
            if(visited[gid])
				continue;
			visited[gid] = 1;
			tgData[tid].push_back(gid);
			dataTG[gid] = tid;
        }
		// tmptg[tid] = tgData[tid];
    }

    in.close();

	cout << " .. Done" << endl;

	// resort(tmptg);

	// ofstream out(thresholdGroupFile, ios::app);

	// for(unsigned i = 0;i < tmptg.size();i++){
	// 	out << i;//保存阈值分组id
	// 	out << "\t";
	// 	for(unsigned j = 0; j < tmptg[i].size();j++){
	// 		out << tmptg[i][j] << ' ';//保存属于同一阈值组的图id
	// 	}
	// 	out << endl << flush;
	// }
	// out.close();
}

void nass::readMCSIndexFile()
{
	DataSet* mcsset = DataSet::getInstance();
	
	if(MCSIndex != NULL) delete [] MCSIndex;
	MCSIndex = new vector< pair<int,int> >[mcsset->numMCSs()];

	bool sanity_check[mcsset->numMCSs()];//表示当前图的索引表是否已被访问过
	for(unsigned i = 0; i < mcsset->numMCSs(); i++)
		sanity_check[i] = false;

	cout << "\nReading index from " << MCSindexfile << endl << flush;

	ifstream in(MCSindexfile);

	if(!in.is_open()){
		cerr << "Error in opening index file" << endl;
	}
		
	string line;

	use_sample = false;
	{
		getline(in, line);
		stringstream ss(line);
		char c; ss >> c;
		if(c != 't'){
			cerr << "Can't read the maximum threshold" << endl << endl;
			cerr << line << endl;
		}
		ss >> mcs_max_threshold;
		ss >> use_sample;
	}

	unsigned long long total_entries = 0;
	while(getline(in, line)){
		stringstream ss(line);
		// thresholdGroup.push_back(vector< int >());
		int cid, gid, dist;
		ss >> cid;

		sanity_check[cid] = true;

		while(ss >> gid){
			ss >> dist;	
			MCSIndex[gid].push_back(pair<int, int>(cid, dist));
			MCSIndex[cid].push_back(pair<int, int>(gid, dist));

			total_entries++;
		}
	}

	in.close();

	// cout << "============= Statistics of the MCS index =============" << endl;
	// cout << "Maximum threshold of the index: " << mcs_max_threshold << endl << flush;
	// cout << "Total number of indexed entries: " << total_entries << endl;
	// cout << "===================================================" << endl << endl;
}

bool mSize(const pair<int, int> &x, const pair<int, int> &y)
{
	return x.second < y.second;
}

#include <unistd.h>
void nass::search(coregraph* x)
{
	DataSet* dataset = DataSet::getInstance();
	DataSet* mcsset = DataSet::getInstance();

	int* visited = new int[dataset->numGraphs()];//图是否被访问
	memset(visited, 0, sizeof(int)*dataset->numGraphs());
	int* cvisited = new int[mcsset->numMCSs()];//类是否被访问
	memset(cvisited, 0, sizeof(int)*mcsset->numMCSs());
	// vector<int> c_cand(mcsset->numMCSs());
	// cged[i].first = cid, cged[i].second = 保存每一类MCS与q之间的ged减去MCS与q的size之间的差值
	vector<pair<int, int>> cged;
	// 保存候选图的id,lb
    vector<pair<int, int> > candidates;
	// 保存候选类
	vector<int> c_candidates;
	cged.reserve(mcsset->numMCSs());
	candidates.reserve(dataset->numGraphs());
	c_candidates.reserve(mcsset->numMCSs());

	// cout<<"mcsset->size = "<<mcsset->numMCSs()<<endl;

	// coregraph* q = dataset->graphAt(20560);
// 	for(int cid = 0; cid < mcsset->numMCSs(); cid++){
// 		if(cid == 20982 || cid == 21248 || cid == 21319 || cid == 21355 || cid == 21783 || cid == 22024 || cid == 22077 || cid == 22144 || cid == 22159
// 		 || cid == 22249 || cid == 22368 || cid == 22394 || cid == 22668 || cid == 22901 || cid == 24058 || cid == 24691 || cid == 24755 || cid == 25054
// 		  || cid == 25114 || cid == 25185 || cid == 25372 || cid == 25507 || cid == 25516 || cid == 25538 || cid == 25745 || cid == 25765 || cid == 25897
// 		   || cid == 26155 || cid == 26722 || cid == 26946 || cid == 27044 || cid == 27240 || cid == 27401 || cid == 28017 || cid == 28060 || cid == 28282
// 		    || cid == 28427 || cid == 28640 || cid == 28674 || cid == 29002 || cid == 29182 || cid == 29264 || cid == 29398 || cid == 29412 || cid == 29653
// 			 || cid == 29800 || cid == 29820 || cid == 29888 || cid == 30489 || cid == 30503){
// 	// coregraph* g = dataset->graphAt(31571);
// 	coregraph* q = dataset->graphAt(20560);
// 	for(int j = 0; j < tgData[cid].size(); j++){
// 	int gid = tgData[cid][j];
// 	coregraph* g = dataset->graphAt(gid);
// 	// coregraph* g = mcsset->mcsAt(20982);
// 	int lb = basicFilter(q, g);
// 	cout<<"gid = "<<gid<<", lb = "<<lb<<endl;
// // 	// if q->vsize < g->vsize, swap(q, g)
// 	NassGED* ged = new NassGED(q, g, threshold, timer);
// 	int distance = threshold + 1;
// 	distance = ged->computeGED();
// 	cout<<"distance = "<<distance<<", real_ged = "<<ged->real_ged<<endl;}}
// 	}
// return;
	// basicFilter(x, candidates, threshold);
	// sort(candidates.begin(), candidates.end(), gtlb);
	// int i = 0;
	// while(i < candidates.size()){
	// 	coregraph* y = dataset->graphAt(candidates[i].first);

	// 	NassGED* ged = new NassGED(x, y, threshold, timer);
	// 	int distance = threshold + 1;
	// 	distance = ged->computeGED();
	// 	delete ged;
		
	// 	if(distance <= threshold){
	// 		// cout<<"gid = "<<candidates[i].first<<", ";
	// 		if(use_index) i = regenCandidates(candidates, i+1, candidates[i].first, distance);
	// 		else res_vec.push_back(candidates[i++].first);
	// 	}
	// 	else i++;
	// }
	// return;

	// long long tmp = 0;
	// long long filt = 0;

	for(int i = 0; i < mcsset->numMCSs(); i++){
		if(tgData[i].size() == 1){
			if(cvisited[i])
				continue;
			cvisited[i] = 1;
			int gid = tgData[i][0];
			// if(visited[gid])
			// 	continue;
			// visited[gid] = 1;
			coregraph* y = dataset->graphAt(gid);
			if(x == y) continue;
			int lb = basicFilter(x, y);
			if(lb > threshold)
				continue;
			visited[gid] = 1;
			candidates.push_back(make_pair(gid, lb));
		}
		else{
			if(cvisited[i])
				continue;
			cvisited[i] = 1;
			long avesize = 0;
			long num = 0;
			for(int j = 0; j < tgData[i].size(); j++){
				int gid = tgData[i][j];
				coregraph* y = dataset->graphAt(gid);
				int lb = basicFilter(x, y);
				if(lb <= threshold)
					num++;
				avesize = avesize + y->vsize() + y->esize();
			}
			if(num == 0)
				continue;
			// avesize = avesize / tgData[i].size();
			// if(abs(avesize - x->esize() - x->vsize()) > threshold + 5)
			// 	continue;


			coregraph* mcs = mcsset->mcsAt(i);
			long diff1 = abs(avesize - mcs->esize() - mcs->vsize());
			// int lb = mcsbasicFilter(x, mcs);
			// int diff = lb;
			// lb = lb - abs((long)(x->vsize() + x->esize()) - (tance == 0){
					// 	for(int j = 0; j < tgData[i].size(); j++){
					// 		int gid = tgData[i][j];
					// 		coregraph* y = dataset->graphAt(gid);
					// 		if(x == y)
					// 			continue;
					// 		visited[gid] = 1;
					// 		res_vec.push_back(gid);
					// 	}
					// }long)(mcs->vsize() + mcs->esize()));
			// if(lb > threshold)
			// 	continue;

			// if(mcs->vsize() == 0){
			// 	cged.push_back(make_pair(i, lb));
			// 	// for(int j = 0; j < tgData[i].size(); j++){
			// 	// 	int gid = tgData[i][j];
			// 	// 	if(visited[gid])
			// 	// 		continue;
			// 	// 	visited[gid] = 1;
			// 	// }
			// 	continue;
			// }

			NassGED* ged = new NassGED(x, mcs, threshold, timer);
			int distance = threshold + 1;
			distance = ged->computeMCSGED();
			// tmp++;
			int diff = abs(((long)x->vsize() - (long)mcs->vsize()) + ((long)x->esize() - (long)mcs->esize()));
			// if(i == 30792)
			// 	cout<<"cid = "<<i<<", diff = "<<diff<<", distance  = "<<distance<<", real_ged = "<<ged->real_ged<<endl;
			
			// bool flag = false;
			if(distance < threshold + 1){
				// if(i == 22368)
					// cout<<"cid = "<<i<<", diff1 = "<<diff1<<", diff2 = "<<diff2<<", distance = "<<distance<<endl;
					// flag = true;tance == 0){
					// 	for(int j = 0; j < tgData[i].size(); j++){
					// 		int gid = tgData[i][j];
					// 		coregraph* y = dataset->graphAt(gid);
					// 		if(x == y)
					// 			continue;
					// 		visited[gid] = 1;
					// 		res_vec.push_back(gid);
					// 	}
					// }
				// if(distance < diff)tance == 0){
					// 	for(int j = 0; j < tgData[i].size(); j++){
					// 		int gid = tgData[i][j];
					// 		coregraph* y = dataset->graphAt(gid);
					// 		if(x == y)
					// 			continue;
					// 		visited[gid] = 1;
					// 		res_vec.push_back(gid);
					// 	}
					// }
					// cout<<"cid = "<<i<<", diff = "<<diff<<", distance = "<<distance<<endl;
				// if(distance < diff){
				// 	cged.push_back(make_pair(i, lb));
				// 	for(int j = 0; j < tgData[i].size(); j++){
				// 		int gid = tgData[i][j];
				// 		if(visited[gid])
				// 			continue;
				// 		visited[gid] = 1;
				// 	}
				// 	continue;
				// }
				// else
					// distance = distance - diff;
					// if(distance == 0){
						for(int j = 0; j < tgData[i].size(); j++){
							int gid = tgData[i][j];
							coregraph* y = dataset->graphAt(gid);
							if(x == y)
								continue;
							if(visited[gid])
								continue;
							visited[gid] = 1;
							// int diff1 = y->vsize()-mcs->vsize() + y->esize()-mcs->esize();
							// if(distance + diff1 <= threshold){
							// 	res_vec.push_back(gid);
							// 	continue;
							// }
							// else if(distance -diff1 > threshold)
							// 	continue;
							NassGED* tmp_ged = new NassGED(x, y, threshold, timer);
							int dist = tmp_ged->computeGED();
							if(dist <= threshold){
								res_vec.push_back(gid);
								for(int k = 0; k < nassIndex[gid].size(); k++){
									int id = nassIndex[gid][k].first;
									if(visited[id])
										continue;
									visited[id] = 1;
									int gdist = nassIndex[gid][k].second;
									if(dist + gdist <= threshold)
										res_vec.push_back(id);
									else if(gdist <= dist + threshold){
										coregraph* g = dataset->graphAt(id);
										int tmp_lb = basicFilter(x, g);
										candidates.emplace_back(id, tmp_lb);
									}
								}
								if(dist + threshold <= index_max_threshold){
									memset(cvisited, 1, mcsset->numMCSs());
									memset(visited, 1, dataset->numGraphs());
								}
							}
							delete tmp_ged;
						}
						distance = distance - diff;
					// }
					// 	for(int j = 0; j < tgData[i].size(); j++){
					// 		int gid = tgData[i][j];
					// 		coregraph* y = dataset->graphAt(gid);
					// 		if(x == y)
					// 			continue;
					// 		visited[gid] = 1;
					// 		res_vec.push_back(gid);
					// 	}
				// distance = ged->real_ged - diff;
			}
			else{
				// if(i == 22368) 
				// 	cout<<"cid = "<<i<<", diff1 = "<<diff1<<", diff2 = "<<diff2<<", distance  = "<<distance<<", real_ged = "<<ged->real_ged<<endl;
				if(ged->real_ged == INT_MAX || ged->real_ged < diff || ged->real_ged < distance)
					continue;
					// cout<<"cid = "<<i<<", diff1 = "<<diff1<<", diff2 = "<<diff2<<", distance  = "<<distance<<", real_ged = "<<ged->real_ged<<endl;
			// 	if(ged->real_ged < diff){
			// 		cged.push_back(make_pair(i, lb));
			// 		for(int j = 0; j < tgData[i].size(); j++){
			// 			int gid = tgData[i][j];
			// 			if(visited[gid])
			// 				continue;
			// 			visited[gid] = 1;
			// 		}
			// 		continue;
			// 	}
			// 	else
				for(int j = 0; j < tgData[i].size(); j++){
					int gid = tgData[i][j];
					coregraph* y = dataset->graphAt(gid);
					if(x == y)
						continue;
					if(visited[gid])
						continue;
					int diff1 = y->vsize()-mcs->vsize() + y->esize()-mcs->esize();
					if(distance + diff1 <= threshold){
						res_vec.push_back(gid);
						visited[gid] = 1;
						continue;
					}
					else if(distance -diff1 > threshold){
						visited[gid] = 1;
						continue;
					}
				}
				distance = ged->real_ged - diff;
			}
			delete ged;
			if(distance <= threshold){
				vector<pair<int, int>> vec;
				cged.push_back(make_pair(i, distance));
				// for(int j = 0; j < tgData[i].size(); j++){(flag == 0 || flag == 1) && 
				// 	int gid = tgData[i][j];
				// 	if(visited[gid])
				// 		continue;
				// 	visited[gid] = 1;(flag == 0 || flag == 1) && 
				// }
				if(threshold + distance > mcs_max_threshold){
				// 	if(i == 28634 || i == 22350 || i == 23706 || i == 21918 || i == 20993 || i == 25412 || i == 28641 || i == 21438
				// || i == 27876 || i == 21537 || i == 23572 || i == 20849 || i == 22920 || i == 22468 || i == 21959)
				// // 		cout<<"1: cid = "<<i<<", diff = "<<diff<<", distance  = "<<distance<<", real_ged = "<<ged->real_ged<<endl;
                //     cout<<"cid = ";
					for(int j = 0; j < MCSIndex[i].size(); j++){
                        int cid = MCSIndex[i][j].first;
                        if(cvisited[cid])
                            continue;
                        int dist = MCSIndex[i][j].second;
                        if(dist + distance <= threshold){
							// if(flag)
								// cout<<"cand-cid = "<<cid<<", ";
                            cvisited[cid] = 1;
                            cged.push_back(make_pair(cid, dist + distance));
				// 			if(cid == 28634 || cid == 22350 || cid == 23706 || cid == 21918 || cid == 20993 || cid == 25412 || cid == 28641 || cid == 21438
				// || cid == 27876 || cid == 21537 || cid == 23572 || cid == 20849 || cid == 22920 || cid == 22468 || cid == 21959)
								// cout<<"1: cid = "<<cid<<", ";
                            // for(int k = 0; k < tgData[cid].size(); k++){
                            //     int gid = tgData[cid][k];
                            //     visited[gid] = 1;
                            // }
                        }
                        else if(dist > distance + threshold){
							// if(flag)
								// cout<<"filter-cid = "<<cid<<", ";
                            cvisited[cid] = 1;
				// 			if(cid == 28634 || cid == 22350 || cid == 23706 || cid == 21918 || cid == 20993 || cid == 25412 || cid == 28641 || cid == 21438
				// || cid == 27876 || cid == 21537 || cid == 23572 || cid == 20849 || cid == 22920 || cid == 22468 || cid == 21959)
								// cout<<"2: cid = "<<cid<<", ";
						}
                    }
				// 	if(i == 28634 || i == 22350 || i == 23706 || i == 21918 || i == 20993 || i == 25412 || i == 28641 || i == 21438
				// || i == 27876 || i == 21537 || i == 23572 || i == 20849 || i == 22920 || i == 22468 || i == 21959)
				// 		cout<<endl;
                }
                else{
				// 	if(i == 28634 || i == 22350 || i == 23706 || i == 21918 || i == 20993 || i == 25412 || i == 28641 || i == 21438
				// || i == 27876 || i == 21537 || i == 23572 || i == 20849 || i == 22920 || i == 22468 || i == 21959)
				// 		cout<<"2: cid = "<<i<<", diff = "<<diff<<", distance  = "<<distance<<", real_ged = "<<ged->real_ged<<endl;
                    vector<int> vec;
					for(int j = 0; j < MCSIndex[i].size(); j++){
                        int cid = MCSIndex[i][j].first;
                        if(cvisited[cid])
                            continue;
                        int dist = MCSIndex[i][j].second;
                            
                        if(dist + distance <= threshold){
							// if(flag)
								// cout<<"cand-cid = "<<cid<<", ";
                            cvisited[cid] = 1;
                			cged.push_back(make_pair(cid, dist + distance));
				// 			if(cid == 28634 || cid == 22350 || cid == 23706 || cid == 21918 || cid == 20993 || cid == 25412 || cid == 28641 || cid == 21438
				// || cid == 27876 || cid == 21537 || cid == 23572 || cid == 20849 || cid == 22920 || cid == 22468 || cid == 21959)
				// 				cout<<"3: cid = "<<cid<<", ";
                            // for(int k = 0; k < tgData[cid].size(); k++){
                            //     int gid = tgData[cid][k];
                            //     visited[gid] = 1;
                           	// }
                        }
                        else if(dist > distance + threshold){
							// if(flag)
								// cout<<"filter-cid = "<<cid<<", ";
                            cvisited[cid] = 1;
				// 			if(cid == 28634 || cid == 22350 || cid == 23706 || cid == 21918 || cid == 20993 || cid == 25412 || cid == 28641 || cid == 21438
				// || cid == 27876 || cid == 21537 || cid == 23572 || cid == 20849 || cid == 22920 || cid == 22468 || cid == 21959)
								// cout<<"4: cid = "<<cid<<", ";
						}
                        else{
                            vec.push_back(cid);
                        }
                    }
                    if(vec.size() > 0){
                        // cout<<"clear\n";
                        memset(cvisited, 1, sizeof(int)*mcsset->numMCSs());
                        for(int i = 0; i < vec.size(); i++){
                            int cid = vec[i];
                            cvisited[cid] = 0;
                        }
                    }
                }
			}
		}
	}
	// cout<<"\ntmp = "<<tmp<<endl;
	// return;
	sort(cged.begin(), cged.end(), mSize);
	// cout<<"\ncid = \n";
	// cout<<"\ngid = \n";
	for(int i = 0; i < cged.size(); i++){
		int cid = cged[i].first;
		// int lb = cged[i].second;
		// cout<<cid<<", ";
		for(int j = 0; j < tgData[cid].size(); j++){
			int gid = tgData[cid][j];
			if(visited[gid])
				continue;
			visited[gid] = 1;
			// cout<<gid<<", ";
			// if(visited[gid] == 1)
			// 	visited[gid] = 0;
			// else
			// 	continue;
			coregraph* g = dataset->graphAt(gid);
			if(x == g)
				continue;
			int lb = basicFilter(x, g);
			if(lb > threshold)
				continue;
			candidates.emplace_back(gid, lb);
		}
	}
	// cout<<"\n";
	sort(candidates.begin(), candidates.end(), gtlb);
	// cout<<"candidates = ";
	// for(int i = 0; i < candidates.size(); i++)
	// 	cout<<candidates[i].first<<", ";
	// cout<<endl;

	for(int i = 0; i < candidates.size();){
		coregraph* y = dataset->graphAt(candidates[i].first);

		NassGED* ged = new NassGED(x, y, threshold, timer);
		int distance = threshold + 1;
		distance = ged->computeGED();
		delete ged;
		
		if(distance <= threshold){
			if(use_index) i = regenCandidates(candidates, i+1, candidates[i].first, distance);
			else res_vec.push_back(candidates[i++].first);
		}
		else i++;
	}
}
