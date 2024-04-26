#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstring>
#include <cstdlib>

#include "nassged.h"
#include "inves.h"
#include "datafile.h"
#include "dataset.h"
#include "nass.h"
#include "mcs.h"
//#include "inves.h"

using namespace std;

void exit_with_usage()
{
	cerr << "usage: nass-index threshold data_file index_file [options]" << endl;
	cerr << "usage[client]: nass-index data_file server_addr" << endl;
	exit(1);
}

int main(int argc, char* argv[])
{
	if(argc < 3) exit_with_usage();

	bool distributed = false;
	for(unsigned i = 0; argv[1][i] != '\0'; i++)
		if(argv[1][i] >= '0' && argv[1][i] <= '9') continue;
		else{ distributed = true; break; }

	if(distributed && argc != 3) exit_with_usage();
	if(!distributed && argc < 8) exit_with_usage();

	char* server_address = NULL;
	char* data_file = NULL;
	char* index_file = NULL;
	char* thresholdGroup_file = NULL;
        char* resort_graph_file = NULL;
	char* mcs_original_file = NULL;
	char* mcs_file = NULL;
	char* heuristic = NULL;//min_max or min_product
	if(distributed){ data_file = argv[1]; server_address = argv[2]; }
	else{ nass::threshold = atoi(argv[1]); data_file = argv[2]; index_file = argv[3]; thresholdGroup_file = argv[4]; resort_graph_file = argv[5]; mcs_original_file = argv[6]; mcs_file = argv[7]; heuristic = argv[8];}

	int sampling_rate = 100;
	for(int i = 11; i < argc; i += 2){
		if(strcmp(argv[i], "-M") == 0)
			MEMLIMIT = atoi(argv[i+1]);
		else if(strcmp(argv[i], "-p") == 0)
			NTHREADS = atoi(argv[i+1]);
		else if(strcmp(argv[i], "-s") == 0)
			sampling_rate = atoi(argv[i+1]);
		else if(strcmp(argv[i], "--coordinator") == 0){
			distributed = true;
			i -= 1;
		}
		else{
			cerr << "Unknown option: " << argv[i] << endl;
			cerr << "usage: nass-index threshold data_file index_file thresholdGroup_file resort_graph_file mcs_original_file mcs_file heuristic[options]" << endl;
			exit(1);
		}
	}

	DataFile* file = new SyntheticFile(data_file);
	DataSet* dataset = DataSet::getInstance();
	dataset->buildDataSet(*file);
	// vector<coregraph*> data = dataset->getData();
	// cout << "hello" << endl;
	// for (unsigned i = 0; i < data.size(); i++)
	// {
	// 	cout << data[i]->vsize() << ' ';
	// }
	// cout<<endl;
	file->close();
	delete file;

	nass builder(index_file, thresholdGroup_file, false);
	builder.buildNassIndex(distributed, server_address, sampling_rate);
	unsigned mcsData_size = builder.getThresholdGroupSize();
	// unsigned mcsData_size = 13473;
	mcs extracter(thresholdGroup_file, resort_graph_file, mcs_original_file, mcs_file, heuristic, mcsData_size);
	extracter.buildMCSGroup();
	

	// all graph (wrapper) objects should be destroyed before calling finishUp()
	DataSet::finishUp();

	return 0;
}
