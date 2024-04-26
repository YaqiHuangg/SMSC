#ifndef __WORKLOAD_H
#define __WORKLOAD_H

#include <vector>

#include "dataset.h"
#include "coregraph.h"

class Workload
{
private:
    vector<coregraph*> workload;
    bool loaded_from_a_custom_file;

public:
    Workload(){
        loaded_from_a_custom_file = false;
    }

    ~Workload(){
        if(!loaded_from_a_custom_file) return;

        for(unsigned i = 0; i < workload.size(); i++)
            delete workload[i];
    }

    void generateWorkload(){
        loaded_from_a_custom_file = false;
        DataSet* dataset = DataSet::getInstance();
        // vector<coregraph*> data = dataset->getData();

        for(unsigned i = 0; i < dataset->numGraphs() && workload.size() < 100; i++){
            coregraph* g = dataset->graphAt(i);
            if(g->vsize() >= 60 && g->vsize() < 70 && g != NULL)
                workload.push_back(g);
        }
        
        // srand(3123);
        // srand(3120);
        // srand(99);
        // while(workload.size() < 100){
        //     int gid = rand() % dataset->numGraphs();
        //     coregraph* g = dataset->graphAt(gid);
        //     // coregraph* g = data[gid];
        //     if(g != NULL) workload.push_back(g);
        // }
    }

    void readWorkload(char* filename){
        loaded_from_a_custom_file = true;

        unsigned sz_v = vlabel_map.size();
        unsigned sz_e = elabel_map.size();

        DataFile* file = new SyntheticFile(filename);
        coregraph* g = file->getNextGraph();
        while(g != NULL){
            workload.push_back(g);
            g = file->getNextGraph();
        }

        file->close();
        delete file;

        for(unsigned i = 0; i < workload.size(); i++)
            workload[i]->labelScan();

        if(sz_v == vlabel_map.size() && sz_e == elabel_map.size()) return;

		// if new labels are added from queries, rescan label frequenices for each data.
        DataSet* dataset = DataSet::getInstance();
        // vector<coregraph*> data = dataset->getData();
        // dataset->resort(data);
        for(unsigned i = 0; i < dataset->numGraphs(); i++){
            coregraph* g = dataset->graphAt(i);
            // coregraph* g = data[i];
            g->labelRescan(sz_v, sz_e);
        }
    }

    coregraph* operator[](unsigned idx){
        return workload[idx];
    }

    unsigned size(){ return workload.size(); }
};

#endif // __WORKLOAD_H
