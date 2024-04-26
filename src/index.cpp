#include "nass.h"
#include "nassged.h"
#include <cstring>
#include <sstream>
#include <thread>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <netdb.h>

using namespace std;

#define IDONLY -1
#define GETSTARTED -2
#define TERMINATE -4

#define PORTNO 3333

int NTHREADS = 8; // default number of index-building threads: 8
unsigned MEMLIMIT = 1000; // default memory consumption for index building: 1G

static unsigned memory_limit = 0;
static unsigned minimum_memory = 0;


/*********** UTILITY FUNCTIONS FOR NETWORKING ***********/
void error(const char *msg)
{
	cerr << msg << endl;
    exit(0);
}

void nread(int fd, char* datum, size_t sz)
{
	int n = read(fd, datum, sz);
	if(n < 0) error("network read error");
}

void nwrite(int fd, char* datum, size_t sz)
{
	int n = write(fd, datum, sz);
	if(n < 0) error("network write error");
}
/******************************************************/

void nass::resort(vector<vector< int >>& thresholdGroup)
{
	//按数组长度升序排列
    sort(thresholdGroup.begin(),thresholdGroup.end(), groupLength);
	
    // 删除重复项
	vector<vector<int>>::iterator iter;
	for(unsigned i=0;i<thresholdGroup.size();){
		auto iter = thresholdGroup.begin() + i;
		if(find(thresholdGroup.begin(), iter, thresholdGroup[i]) != iter)
			thresholdGroup.erase(iter);
		else
			i++;
        // for(unsigned j=i+1;j<thresholdGroup.size();j++){
        //     if(thresholdGroup[i]==thresholdGroup[j]){
        //         iter=thresholdGroup.begin()+j;
        //         thresholdGroup.erase(iter);
        //         j--;
        //     }
        // }
    }

	// vector<vector< int >> tmp;
	// tmp.reserve(thresholdGroup.size());
	// for(unsigned i=0;i<thresholdGroup.size();){
	// 	size_t length1=thresholdGroup[i].size();
    //     bool flag=false;//标记i是否为其他数组子集
	// 	// vector<int> pos;//标记包含i的数组的位置
    //     for(unsigned j=i+1;j<thresholdGroup.size();j++){
    //         size_t m=0;
    //         size_t n=0;
    //         size_t length2=thresholdGroup[j].size();
    //         vector<int> temp(length2,0);//记录j中i子集的位置
	// 		// auto iter = thresholdGroup.begin() + j;
    //         if(length1 == length2){
    //             if(thresholdGroup[i] == thresholdGroup[j]){
	// 				thresholdGroup.erase(thresholdGroup.begin() + j);
	// 				j--;
	// 			}
	// 		}
    //         else{//j的长度>=i的长度
    //             while(m<length1&&n<length2){
    //                 if(thresholdGroup[j][n]<thresholdGroup[i][m])
    //                     n++;
    //                 else if(thresholdGroup[i][m]==thresholdGroup[j][n]){
    //                     temp[n]=1;
    //                     // count--;
    //                     m++;
    //                     n++;
    //                 }
    //                 else if(thresholdGroup[j][n]>thresholdGroup[i][m])
    //                     break;
    //             }
    //             if(m<length1)//i不是j的子集
    //                 continue;
    //             else if(m==length1){//i是j的子集
    //                 if(length1==1){
    //                     flag=true;
    //                     break;
    //                 }
    //                 n=0;
    //                 flag=true;
    //                 while(true){
	// 					// pos.push_back(j);
    //                     if(n<temp.size()&&temp[n]==1){
    //                         //数组中不存在重复项
    //                         thresholdGroup[j].erase(begin(thresholdGroup[j])+n);
    //                         // 此时thresholdGroup数组size已减小，所以temp的size也应减小
    //                         temp.erase(begin(temp)+n);
    //                     }
    //                     else if(n<temp.size())
    //                         n++;
    //                     else
    //                         break;
    //                 }
    //             }
    //         }
    //     }

    //     if(length1<2&&flag){//i只有1个元素，且被包含于其他数组，可删除不在此考虑
	// 		thresholdGroup.erase(thresholdGroup.begin()+i);//删除空数组
    //     }
    //     else{
	// 		tmp.push_back(thresholdGroup[i]);
	// 		thresholdGroup.erase(thresholdGroup.begin()+i);
    //     }
	// }
	// swap(thresholdGroup, tmp);
}

// bool iferase;
// void nass::recluster(int tid)
// {
// 	iferase = false;
// 	if(thresholdGroup[tid].size()==2){
//         if(count(index[thresholdGroup[tid][1]].begin(), index[thresholdGroup[tid][1]].end(), thresholdGroup[tid][0])==0){
//             vector<int> temp1;
//             vector<int> temp2;
//             temp1.push_back(thresholdGroup[tid][0]);
//             temp2.push_back(thresholdGroup[tid][1]);
//             thresholdGroup.erase(thresholdGroup.begin()+tid);
//             thresholdGroup.insert(thresholdGroup.begin()+tid,temp1);
//             thresholdGroup.insert(thresholdGroup.begin()+tid,temp2);
//         }
//         else{
//             newThresholdGroup.push_back(thresholdGroup[tid]);
//             thresholdGroup.erase(thresholdGroup.begin()+tid);
//             iferase = true;
//         }
//     }
//     else{
//         vector<int> temp;
// 	    temp.push_back(thresholdGroup[tid][0]);
// 	    thresholdGroup[tid].erase(thresholdGroup[tid].begin());
// 	    bool flag = true;//记录gid是否与temp中每个图的ged都在5之内
// 	    for(unsigned i = 0; i < thresholdGroup[tid].size(); i++){
// 		    int gid = thresholdGroup[tid][i];
// 		    for(unsigned j = 0; j< temp.size(); j++)
// 			    if(count(index[temp[j]].begin(), index[temp[j]].end(), gid)==0){
// 				    flag=false;
// 				    break;
// 			    }
// 		    if(flag){
// 			    temp.push_back(gid);
// 			    thresholdGroup[tid].erase(thresholdGroup[tid].begin()+i);
//                 i--;//因为gid已从group中移除
// 		    }
// 		    else
// 			    flag = true;
// 	    }
//         newThresholdGroup.push_back(temp);
//     }
// }
void nass::recluster(int tid)
{
	if(thresholdGroup[tid].size()==1)
        return;
    else if(thresholdGroup[tid].size()==2){
		int id1 = thresholdGroup[tid][1];
		int id2 = thresholdGroup[tid][0];
        if(count(index[id1].begin(), index[id1].end(), id2)==0){
            vector<int> temp1;
            vector<int> temp2;
            temp1.push_back(id2);
            temp2.push_back(id1);
            thresholdGroup.erase(thresholdGroup.begin()+tid);
            thresholdGroup.insert(thresholdGroup.begin()+tid,temp1);
            thresholdGroup.insert(thresholdGroup.begin()+tid,temp2);
        }
    }
    else{
        vector<int> temp;
	    temp.push_back(thresholdGroup[tid][0]);
	    thresholdGroup[tid].erase(thresholdGroup[tid].begin());
	    bool flag = true;//记录gid是否与temp中每个图的ged都在5之内
	    for(unsigned i = 0; i < thresholdGroup[tid].size(); i++){
		    int gid = thresholdGroup[tid][i];
		    for(unsigned j = 0; j< temp.size(); j++){
				int id = temp[j];
			    if(count(index[id].begin(), index[id].end(), gid)==0){
				    flag=false;
				    break;
			    }
			}
		    if(flag){
			    temp.push_back(gid);
			    thresholdGroup[tid].erase(thresholdGroup[tid].begin()+i);
                i--;//因为gid已从group中移除
		    }
		    else
			    flag = true;
	    }
        // newThresholdGroup.push_back(temp);
		if(thresholdGroup[tid].size()==0)
            thresholdGroup.erase(thresholdGroup.begin()+tid);
		thresholdGroup.insert(thresholdGroup.begin() + tid, temp);
    }
}

void nass::refine()
{
	//细化阈值分组

	for(unsigned i=0;i<thresholdGroup.size(); i++){
		recluster(i);
	}

	resort(thresholdGroup);

	static unsigned progress = 0;
	unsigned total_size = thresholdGroup.size();
	for(unsigned i=0, sum=0;i<thresholdGroup.size();sum++){
		unsigned percent = sum*100/(total_size-1);
		if(percent >= progress){
			// cout << "\rRefineing: " << percent << "%";
			cout <<"\r refined "<<sum;
			cout << flush;
			progress = percent;
		}

        size_t length1=thresholdGroup[i].size();
        bool flag=false;//标记i是否为其他数组子集
		// vector<int> pos;//标记包含i的数组的位置
        for(unsigned j=i+1;j<thresholdGroup.size();j++){
            size_t m=0;
            size_t n=0;
            size_t length2=thresholdGroup[j].size();
            vector<int> temp(length2,0);//记录j中i子集的位置
			// auto iter = thresholdGroup.begin() + j;
            if(length1 == length2){
                if(thresholdGroup[i] == thresholdGroup[j]){
					thresholdGroup.erase(thresholdGroup.begin() + j);
					j--;
				}
			}
            else{//j的长度>=i的长度
                while(m<length1&&n<length2){
                    if(thresholdGroup[j][n]<thresholdGroup[i][m])
                        n++;
                    else if(thresholdGroup[i][m]==thresholdGroup[j][n]){
                        temp[n]=1;
                        // count--;
                        m++;
                        n++;
                    }
                    else if(thresholdGroup[j][n]>thresholdGroup[i][m])
                        break;
                }
                if(m<length1)//i不是j的子集
                    continue;
                else if(m==length1){//i是j的子集
                    if(length1==1){
                        flag=true;
                        break;
                    }
                    n=0;
                    flag=true;
                    while(true){
						// pos.push_back(j);
                        if(n<temp.size()&&temp[n]==1){
                            //数组中不存在重复项
                            thresholdGroup[j].erase(begin(thresholdGroup[j])+n);
                            // 此时thresholdGroup数组size已减小，所以temp的size也应减小
                            temp.erase(begin(temp)+n);
                        }
                        else if(n<temp.size())
                            n++;
                        else
                            break;
                    }
                }
            }
        }

        if(length1<2&&flag){//i只有1个元素，且被包含于其他数组，可删除不在此考虑
			thresholdGroup.erase(thresholdGroup.begin()+i);//删除空数组
        }
        else{
			newThresholdGroup.push_back(thresholdGroup[i]);
			thresholdGroup.erase(thresholdGroup.begin()+i);
			// bool ifresort = false;
			// recluster(i);
            // if(pos.size()>0){
			// 	for(unsigned p = pos.size()-1; p >= 0; p--)
            //         if(pos[p]-1>0&&thresholdGroup[pos[p]-1].size()<thresholdGroup[pos[p]-2].size()){
            //             ifresort = true;
            //             break;
            //         }
            // }
			// if(ifresort)
            //     resort(thresholdGroup);
        }
	}
}

void nass::miniRefine(vector<vector< int >>& group){
	for(unsigned i=0;i<group.size();){
        long length1=group[i].size();
        bool flag=false;//标记i是否为其他数组子集
        for(unsigned j=i+1;j<group.size();j++){
            long m=0;
            long n=0;
            long length2=group[j].size();
            if(length1==length2)
                continue;//不存在重复数组
            else{//j的长度>=i的长度
                while(m<length1&&n<length2){
                    if(group[j][n]<group[i][m])
                        n++;
                    else if(group[i][m]==group[j][n]){
                        m++;
                        n++;
                    }
                    else if(group[j][n]>group[i][m])
                        break;
                }
                if(m<length1)//i不是j的子集
                    continue;
                else if(m==length1){//i是j的子集
                    flag=true;
                    break;
                }
            }
        }
        if(flag)
            group.erase(group.begin()+i);
        else
            i++;
    }
}

void nass::outputThresholdGroup(vector<vector<int>>& newThresholdGroup)
{
	ofstream out(thresholdGroupFile, ios::app);

	for(unsigned i=0;i<newThresholdGroup.size();i++){
		out << i;//保存阈值分组id
		out << "\t";
		for(unsigned j=0;j<newThresholdGroup[i].size();j++){
			out << newThresholdGroup[i][j] << ' ';//保存属于同一阈值组的图id
		}
		out << endl << flush;
	}
	out.close();
}

bool nass::loadNassIndex(bool verbose)
{
	if(nassIndex) delete [] nassIndex;
	if(index) delete [] index;

	DataSet* dataset = DataSet::getInstance();
	index = new vector<int>[dataset->numGraphs()];
	nassIndex = new vector< pair<int,int> >[dataset->numGraphs()];
	sampled = new bool[dataset->numGraphs()];
	thresholdGroup.reserve(dataset->numGraphs());

	bool sanity_check[dataset->numGraphs()];//表示当前图的索引表是否已被访问过
	for(unsigned i = 0; i < dataset->numGraphs(); i++)
		sanity_check[i] = false;

	if(verbose)
		cout << "\nReading index from " << indexfile << endl << flush;

	ifstream in(indexfile);

	if(!in.is_open()){
		cerr << "Error in opening index file" << endl;
		return false;
	}
		
	string line;

	use_sample = false;
	// first read the max threshold
	{
		getline(in, line);
		stringstream ss(line);
		char c; ss >> c;
		if(c != 't'){
			cerr << "Can't read the maximum threshold" << endl << endl;
			cerr << line << endl;
			return false;
		}
		ss >> index_max_threshold;
		ss >> use_sample;
	}

	unsigned long long total_entries = 0;
	unsigned long long inexact_entries = 0;
	while(getline(in, line)){
		stringstream ss(line);
		// thresholdGroup.push_back(vector< int >());
		int cid, gid, dist;
		ss >> cid;

		sanity_check[cid] = true;

		int sampling = 1;
		if(use_sample) ss >> sampling;
		sampled[cid] = (sampling == 1) ? true : false;

		while(ss >> gid){
			ss >> dist;	
			nassIndex[gid].push_back(pair<int, int>(cid, dist));
			nassIndex[cid].push_back(pair<int, int>(gid, dist));

			total_entries++;
			if(dist < 0) inexact_entries++;
		}
	}

	in.close();

	// sanity check!! and sorting ids
	for(unsigned id = 0; id < dataset->numGraphs(); id++){
		vector<int> temp;
		temp.push_back(id);
		if(!sanity_check[id]){
			cout << "Sanity check for index failed! " << id << endl;
			return false;
		}
		std::sort(nassIndex[id].begin(), nassIndex[id].end(), gtid);

		index[id].reserve(nassIndex[id].size());
		for(unsigned i = 0; i < nassIndex[id].size(); i++){
			index[id].push_back(nassIndex[id][i].first);
			temp.push_back(nassIndex[id][i].first);
		}
		
		sort(temp.begin(), temp.end(),numSize);
		thresholdGroup.push_back(temp);
	}

	// 细化阈值分组
	resort(thresholdGroup);
	newThresholdGroup.reserve(thresholdGroup.size());
	// cout << "\rRefineing: 0%" << flush;

	refine();
	cout<<"...Done!"<<endl;

	resort(newThresholdGroup);
miniRefine(newThresholdGroup);

	vector<vector<int>> group;
	group.reserve(newThresholdGroup.size());
	int* visited = new int[dataset->numGraphs()];
	memset(visited, 0, sizeof(int)*dataset->numGraphs());
	for(unsigned i = 0; i < newThresholdGroup.size(); i++){
		vector<int> tmp;
		tmp.reserve(newThresholdGroup[i].size());
		for(unsigned j = 0; j < newThresholdGroup[i].size(); j++){
			int gid = newThresholdGroup[i][j];
			if(visited[gid]) continue;
			visited[gid] = 1;
			tmp.push_back(gid);
            }
            if(tmp.size() > 0)
			group.push_back(tmp);
        }
        resort(group);
	miniRefine(group);

	outputThresholdGroup(group);

	// if(verbose){
	// 	cout << "============= Statistics of the index =============" << endl;
	// 	cout << "Maximum threshold of the index: " << index_max_threshold << endl << flush;
	// 	cout << "Total number of indexed entries: " << total_entries << endl;
	// 	cout << "Number of inexact entries: " << inexact_entries << endl;
	// 	cout << "===================================================" << endl << endl;
	// }

	return true;
}

int nass::connectToCoordinator(const char* address)
{
	if(address == NULL) address = server_address;

    struct sockaddr_in serv_addr;
    struct hostent *server;

    int sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock < 0) error("ERROR opening socket");

    server = gethostbyname(address);
    if (server == NULL) error("ERROR, no such host\n"); 

    bzero((char *) &serv_addr, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    bcopy((char *)server->h_addr, (char *)&serv_addr.sin_addr.s_addr, server->h_length);
    serv_addr.sin_port = htons(PORTNO);

	int res = -1;
	for(int i = 0; i < 5 && res < 0; i++) // retry up to five times
		res = connect(sock, (struct sockaddr *) &serv_addr, sizeof(serv_addr));
	if(res < 0) error("ERROR connecting");

	return sock;
}

int nass::initCoordinator(struct sockaddr_in& serv_addr, struct sockaddr_in& cli_addr)
{
	int sock =  socket(AF_INET, SOCK_STREAM, 0);
	if (sock < 0) error("ERROR opening socket");

	bzero((char *) &serv_addr, sizeof(serv_addr));

	serv_addr.sin_family = AF_INET;  
	serv_addr.sin_addr.s_addr = INADDR_ANY;  
	serv_addr.sin_port = htons(PORTNO);

	int res = 1;
	setsockopt(sock, SOL_SOCKET, SO_REUSEADDR, &res, sizeof(int));

#if defined(__APPLE__)
	res = ::bind(sock, (struct sockaddr *) &serv_addr, sizeof(serv_addr));
#elif defined(__linux__)
	res = bind(sock, (struct sockaddr *) &serv_addr, sizeof(serv_addr));
#endif
	if(res < 0) error("ERROR on binding");

	listen(sock, 5);

	return sock;
}

void nass::coordinateLoop()
{
	struct sockaddr_in serv_addr, cli_addr;
	int sock = initCoordinator(serv_addr, cli_addr);

	vector<pair<int,int> > res;
	while(true){
		socklen_t clilen = sizeof(cli_addr);
		int clisock = accept(sock, (struct sockaddr *) &cli_addr, &clilen);
		if (clisock < 0){
			cerr << "ERROR on accept" << endl;
			close(clisock); break;
		}

		int recv;
		nread(clisock, (char *)&recv, sizeof(int));

		if(recv == TERMINATE){ close(clisock); break; }
		if(recv == GETSTARTED){
			DataSet* dataset = DataSet::getInstance();

			nwrite(clisock, (char *)&nass::threshold, sizeof(int));
			nwrite(clisock, (char *)&MEMLIMIT, sizeof(int));

			int one = 1, zero = 0;
			for(unsigned i = 0; i < dataset->numGraphs(); i++)
				if(sampled[i]) nwrite(clisock, (char *)&one, sizeof(int));
				else nwrite(clisock, (char *)&zero, sizeof(int));

			close(clisock);
			continue;
		}

		if(recv != IDONLY){
			int sz;
			nread(clisock, (char *)&sz, sizeof(int)); // number of results

			int gr, dist;
			for(int i = 0; i < sz; i++){
				nread(clisock, (char *)&gr, sizeof(int));
				nread(clisock, (char *)&dist, sizeof(int));
				res.push_back(pair<int, int>(gr, dist));
			}

			outputResult(res, recv);
			res.clear();
		}

		int gid = getNextGraphID();
		send(clisock, (char *)&gid, sizeof(int), 0);

		close(clisock);
	}

	close(sock);
}

int nass::getNextGraphIDRemote(vector<pair<int, int> >* res, int id)
{
	id_mutex.lock();
	int sock = connectToCoordinator();

	nwrite(sock, (char *)&id, sizeof(int));

	if(id != IDONLY){
		int n = (int)res->size();
		nwrite(sock, (char *)&n, sizeof(int));

		for(int i = 0; i < n; i++){
			nwrite(sock, (char *)&(*res)[i].first, sizeof(int));
			nwrite(sock, (char *)&(*res)[i].second, sizeof(int));
		}
	}

	int gid = -1;
	nread(sock, (char *)&gid, sizeof(int));

	close(sock);
	id_mutex.unlock();

	return gid;
}

int nass::getNextGraphIDLocal()
{
	id_mutex.lock();
	int id = -1;

	DataSet* dataset = DataSet::getInstance();
	vector<coregraph*> data = dataset->getData();
	// dataset->resort(data);
	// if(next_graph_id == dataset->numGraphs()){
	// 	id_mutex.unlock();
	// 	return -1;
	// }
	if(next_graph_id == data.size()){
		id_mutex.unlock();
		return -1;
	}

	id = next_graph_id++;

	id_mutex.unlock();
	return id;
}

int nass::getNextGraphID(vector<pair<int, int> >* res, int id)
{
	if(remote) return getNextGraphIDRemote(res, id);
	else return getNextGraphIDLocal();
}

void nass::outputResult(vector<pair<int, int> >& res, int id)
{
//	bool verbose = true;
	static unsigned num_output = 0;
	static unsigned progress = 0;

	out_mutex.lock();

	// ofstream out(indexfile, fstream::app);
	ofstream out(indexfile, ios::app);

	out << id;
	if(use_sample) out << " " << sampled[id];
	out << "\t";

	for(unsigned j = 0; j < res.size(); j++){
		out << res[j].first << " " << res[j].second << " ";//res[j].first=graphid, res[j].second=distance
	}
	out << endl << flush;
	out.close();

	num_output++;
	DataSet* dataset = DataSet::getInstance();
	unsigned percent = num_output*100/dataset->numGraphs();
	if(percent > progress){
		cout << "\rProcessing: " << percent << "%";
	//	if(verbose) cout << " (" << num_output << "/" << dataset->numGraphs() << ")";
		cout << flush;
		progress = percent;
	}

	if(num_output == dataset->numGraphs()){
		// GETTIME(&ts2);
		// cout << endl << "Indexing Time: " << ts2.tv_sec - ts1.tv_sec << " seconds" << endl;

		if(distributed && !remote){
			int sock = connectToCoordinator("localhost");
			int terminate = TERMINATE;
			nwrite(sock, (char *)&terminate, sizeof(int)); // force to terminate coordinating server
			close(sock);
		}
	}

	out_mutex.unlock();
}

int nass::getWorkerID()
{
	int id = -1;
	id_mutex.lock(); // share id_mutex with getNextID()

	id = next_worker_id++;

	id_mutex.unlock(); // since this function called only once for each thread
	return id;
}

// bool nass::groupLength(vector<pair<int, int>>& x, vector<pair<int, int>>& y)
// {
// 	return x.size() < y.size();
// }

// void nass::resort(vector<vector<pair<int, int>>>& resortThreasholdGroup)
// {
//     DataSet* dataset = DataSet::getInstance();
// 	unsigned groupNum=dataset->numGraphs();

// 	// for(unsigned i=0;i<groupNum;i++){
// 	// 	std::sort(resortThreasholdGroup[i].begin(), resortThreasholdGroup[i].end(), groupLength);
// 	// }
// 	std::sort(resortThreasholdGroup.begin(), resortThreasholdGroup.end(), groupLength);

// 	for(unsigned i=0;i<groupNum;i++){
// 		if(resortThreasholdGroup[i].size()>1){
// 			;
// 		}
// 	}
// }

void nass::searchWorker()
{
//	static int progress = 0; // progress indicator
	int wid = getWorkerID();

	vector<pair<int, int> > res;

	DataSet* dataset = DataSet::getInstance();
	int gid = getNextGraphID();
	while(gid != -1){
		coregraph* x = dataset->graphAt(gid);

		int threshold = nass::threshold - 1;
		if(sampled[gid]) threshold = nass::threshold;
		vector<pair<int, int> > candidates;
		basicFilter(x, candidates, threshold, gid+1);
		for(unsigned i = 0; i < candidates.size(); i++){
			coregraph* y = dataset->graphAt(candidates[i].first);
			// coregraph* y = data[candidates[i].first];

			if(sampled[gid] || sampled[candidates[i].first]) threshold = nass::threshold;
			else threshold = nass::threshold - 1;

			NassGED* ged = new NassGED(x, y, threshold, NULL, wid);
			int distance = threshold + 1;
			distance = ged->computeGED();
			delete ged;
		
			if(distance <= threshold)
				res.push_back(pair<int,int>(candidates[i].first, distance));
		}
		if(!remote) 
			outputResult(res, gid);
		
/*
		else{
			out_mutex.lock();
			cout << "\r";
			progress = (progress + 1) % 20;
			for(int i = 0; i < progress; i++) cout << "-";
			cout << "*";
			for(int i = progress + 1; i < 20; i++) cout << "-";
			cout << " " << flush;
			out_mutex.unlock();
		}
*/
		gid = getNextGraphID(&res, gid);
		res.clear();
	}
	// reuse id_mutex to decrease number of workers
	id_mutex.lock();
	num_index_worker -= 1;

	if(num_index_worker != 0){
		int factor = NTHREADS / num_index_worker;
		memory_limit = MEMLIMIT/factor;	
		if(memory_limit < minimum_memory) memory_limit = minimum_memory;
	}

	id_mutex.unlock();
	// vector<coregraph*>().swap(data);
	// vector<pair<unsigned, unsigned>>().swap(datagroup);
}

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>
#endif

unsigned getRSS()
{
#if defined(__APPLE__) && defined(__MACH__)
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
        (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
        return 0;

    unsigned long long mem = (unsigned long long)info.resident_size;
    mem = mem / (1024L * 1024L); // size in MB

    return (unsigned)mem;
#endif

    unsigned rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        for(unsigned i = 0; i < 23; i++) ifs >> ignore;
        ifs >> rss;
    }

    unsigned page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024;
    return (rss * page_size_kb) / 1024; // size in MB
}

void nass::monitor()
{
	unsigned margin = 0;

	while(!done_indexing){
		usleep(1000000); // sleep 1 seconds for checking memory consumption

		unsigned usage = getRSS();

		if(usage > memory_limit + margin){
			NassGED::interrupted = true; // signal workers to stop

			int wait_count = 0;
			while(wait_count < num_index_worker){ 
				// wait until all workers stop
				// busy waiting with 1 microsecond sleep
				wait_count = 0;
				for(int i = 0; i < NTHREADS; i++)
					if(NassGED::wait[i]) wait_count++;
				usleep(1);
			}

			if(usage > memory_limit) margin = memory_limit/10;

			NassGED::victim = 0;
			for(int wid = 1; wid < NTHREADS; wid++)
				if(PriorityQueue::total_size[wid] > PriorityQueue::total_size[NassGED::victim])
					NassGED::victim = wid;

			NassGED::interrupted = false;
			sleep(1); // wait another second before finding another victim
		} 
		else if(usage < memory_limit) margin = 0;
	}
}
		
void nass::buildNassIndex(bool distributed, char* address, int sampling_rate)
{
	this->distributed = distributed;
	this->use_sample = (sampling_rate == 100) ? false : true;

	remote  = false;
	server_address = address;

	// distributed setting
	std::thread* coordinator = NULL;
	if(distributed){
		if(server_address) remote = true; // client
		else coordinator = new thread(&nass::coordinateLoop, this); // server
	}

	if(use_sample && !remote){
		if(!loadNassIndex(false)){  // load an index for sampling
			cerr << "Can't load the index for sampling" << endl << flush;
			exit(1);
		}
		use_sample = true; // loadNassIndex can make use_sample false
	}

	DataSet* dataset = DataSet::getInstance();
	sampled = new bool[dataset->numGraphs()];
	for(unsigned i = 0; i < dataset->numGraphs(); i++) sampled[i] = false;

	if(remote){
		int tmp = GETSTARTED;
		int sock = connectToCoordinator();

		nwrite(sock, (char *)&tmp, sizeof(int));

		nread(sock, (char *)&tmp, sizeof(int));
		nass::threshold = tmp;
		nread(sock, (char *)&tmp, sizeof(int));
		MEMLIMIT = tmp;

		for(unsigned i = 0; i < dataset->numGraphs(); i++){
			nread(sock, (char *)&tmp, sizeof(int));
			if(tmp == 1) sampled[i] = true;
		}

		close(sock);
	}
	else if(use_sample){
		int s = 0;
		unsigned tried = 0;
		bool visited[dataset->numGraphs()];
		for(unsigned i = 0; i < dataset->numGraphs(); i++) visited[i] = false;
		for(unsigned i = 0; i < (dataset->numGraphs()*sampling_rate)/100; i++){
			do{
				s = rand() % dataset->numGraphs(); 
				if(!visited[s]) tried++;
				visited[s] = true;	
			} while(sampled[s] || (tried != dataset->numGraphs() && nassIndex[s].empty()));
			sampled[s] = true;
		}

		delete [] nassIndex;
		nassIndex = NULL;
	}
	else for(unsigned i = 0; i < dataset->numGraphs(); i++) sampled[i] = true;

	cout << "Building an index";
	if(remote) cout << " as a remote client";
	cout << endl << "threshold: " << threshold << ", # threads: " << NTHREADS << ", memory limit: " << MEMLIMIT << " MB";
	if(!remote && use_sample) cout << ", " << sampling_rate << "% sampling";
	cout << endl;

	GETTIME(&ts1);

	NassGED::initialize(); // we need # threads NassGED arrays

	if(!remote){
		ofstream out(indexfile, fstream::out);//文件以输出方式打开（内存数据输出到文件）
		out << "t\t" << nass::threshold << " " << use_sample << endl << flush;
		out.close();
	}

	next_graph_id = 0;
	next_worker_id = 0;
	std::thread* worker[NTHREADS];

	num_index_worker = NTHREADS;
	memory_limit = MEMLIMIT;
	minimum_memory = getRSS();
    minimum_memory += (MEMLIMIT - minimum_memory)/NTHREADS;

	if(!remote) cout << "\rProcessing: 0%" << flush;
	// cout<<"111111111111111\n";
	for(int i = 0; i < NTHREADS; i++)
		worker[i] = new thread(&nass::searchWorker, this);
	// cout<<"22222222222222222\n";	

	done_indexing = false;
	thread watcher(&nass::monitor, this);

	for(int i = 0; i < NTHREADS; i++){
		worker[i]->join();
		delete worker[i];
	}

	if(coordinator){
		coordinator->join();
		delete coordinator;
	}

	done_indexing = true;
	watcher.join();

	NassGED::finalize();

//	GETTIME(&ts2);
//	cout << endl << "Indexing Time: " << ts2.tv_sec - ts1.tv_sec << " seconds" << endl;

	if(!remote){ 
		loadNassIndex();// check for index
		GETTIME(&ts2);
		cout << endl << "Indexing and Refining Time: " << ts2.tv_sec - ts1.tv_sec << " seconds" << endl;
	} 
	else cout << "\rDONE           " << endl << flush;
}

