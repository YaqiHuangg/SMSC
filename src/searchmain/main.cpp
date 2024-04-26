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

using namespace std;

int main(int argc, char* argv[])
{
	if(argc < 7){
		cerr << "usage: nass threshold data_file tg_file mcs_file mcs_index_file [index_file1 index_file2 query_file]" << endl;
		exit(1);
	}
	// if(argc < 4){
	// 	cerr << "usage: nass threshold data_file index_file workload_file" << endl;
	// 	exit(1);
	// }

	// parse the specified threshold
	nass::threshold = 0;
	for(unsigned i = 0; argv[1][i] != '\0'; i++)
		if(argv[1][i] >= '0' && argv[1][i] <= '9')
			nass::threshold = nass::threshold*10 + (argv[1][i] - '0');
		else{
			cerr << "usage: nass threshold data_file [index_file]" << endl;
			exit(1);
		}

	const char* tgfile = NULL;
	tgfile = argv[3];

	DataFile* file = new SyntheticFile(argv[2]);
	DataSet* dataset = DataSet::getInstance();
	dataset->buildDataSet(*file);
	file->close();
	delete file;

	DataFile* mcsfile = new MCSFile(argv[4]);
	DataSet* mcsset = DataSet::getInstance();
	mcsset->buildMcsSet(*mcsfile);
	mcsfile->close();
	delete mcsfile;

	Workload workload;
	workload.generateWorkload();
	// return 0;
	// To use a custom workload from a file,
	// replace the line above with
	// workload.readWorkload(filename)
	// workload.readWorkload(argv[7]);

	const char* indexfile1 = NULL;
	// const char* indexfile2 = NULL;
	const char* mcsindexfile = NULL;
	mcsindexfile = argv[5];
	if(argc == 7) indexfile1 = argv[6];
	// if(argc == 8) indexfile2 = argv[7];
	// indexfile = argv[3];
	nass searcher(indexfile1, tgfile, mcsindexfile);

	searcher.run(workload);
	cout << endl << flush;

	// all graph (wrapper) objects should be destroyed before calling finishUp()
	DataSet::finishUp();

	return 0;
}
