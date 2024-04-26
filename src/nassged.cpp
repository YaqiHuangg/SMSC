#include <cassert>

#include "nassged.h"
#include "inves.h"

#include <unistd.h>
#include <limits.h>

long long* PriorityQueue::total_size;

LabelFreq* NassGED::freqx;
LabelFreq* NassGED::freqy;

vertex_t* NassGED::ordering;
vertex_t** NassGED::orderingx;
vertex_t** NassGED::orderingy;

bool NassGED::interrupted = false;
int NassGED::victim = 0;

bool* NassGED::wait = NULL; // for debugging

void NassGED::initialize()
{
	NassGED::freqx = new LabelFreq[NTHREADS];
	NassGED::freqy = new LabelFreq[NTHREADS];

	DataSet* dataset = DataSet::getInstance();
	NassGED::ordering = new vertex_t[dataset->vmax()];
	for(vertex_t i = 0; i < (vertex_t)dataset->vmax(); i++) NassGED::ordering[i] = i;

	orderingx = new vertex_t*[NTHREADS];
	orderingy = new vertex_t*[NTHREADS];
	for(int i = 0; i < NTHREADS; i++){
		NassGED::orderingx[i] = new vertex_t[dataset->vmax()];
		NassGED::orderingy[i] = new vertex_t[dataset->vmax()];
		memcpy(NassGED::orderingx[i], NassGED::ordering, sizeof(vertex_t)*dataset->vmax());
		memcpy(NassGED::orderingy[i], NassGED::ordering, sizeof(vertex_t)*dataset->vmax());
	}

	PriorityQueue::total_size = new long long[NTHREADS];

	NassGED::wait = new bool[NTHREADS];
	for(int i = 0; i < NTHREADS; i++) NassGED::wait[i] = false;
}

void NassGED::finalize()
{
	delete [] NassGED::freqx;
	delete [] NassGED::freqy;

	delete [] NassGED::ordering;

	for(int i = 0; i < NTHREADS; i++){ 
		delete [] NassGED::orderingx[i];
		delete [] NassGED::orderingy[i];
	}
	delete [] NassGED::orderingx;
	delete [] NassGED::orderingy;

	delete [] PriorityQueue::total_size;

	delete [] NassGED::wait;
}

void NassGED::saveUnmappedOf(mapping* m)
{
	uint64* key = m->bitmap;
	map_iter it = mtree->find(key);
	if(it != mtree->end()){
		m->unmapped = it->second;
		m->udist = m->unmapped->distance(); // update distance
	}
	else{
		m->unmapped = new unmappedDistance();
		m->unmapped->dist[LF] = m->udist;
		m->unmapped->method = LF;

		uint64* keydup = new uint64[bitmap_size];
		memcpy(keydup, key, bitmap_size * sizeof(keydup));
		mtree->insert(pair<uint64*, unmappedDistance*>(keydup, m->unmapped));
	}
}

unsigned NassGED::calculateDistance(mapping* m, mapping* pm)
{
	m->mdist = pm->mdist + gx->mappedErrors(pm->size, *gy);
	if(m->distance() > threshold) return m->distance();

/*
	gx->countBridgeFrequencies();
	gy->countBridgeFrequencies();
	m->bdist = gx->bridgeErrors(*gy);
*/
	gx->countBridgesIncremental(false);
	gy->countBridgesIncremental(false);

	m->bdist = gx->bridgeErrors(*gy);

	gx->countBridgesIncremental(true);
	gy->countBridgesIncremental(true);

	if(m->distance() > threshold) return m->distance();

/*
	gx->countLabelFrequencies();
	gy->countLabelFrequencies();
	m->udist = gx->unmappedErrors(*gy);
*/
	gx->countLabelsIncremental(false);
	gy->countLabelsIncremental(false);

	m->udist = gx->unmappedErrors(*gy);

	gx->countLabelsIncremental(true);
	gy->countLabelsIncremental(true);

	return m->distance();
}

void NassGED::expandState(mapping* pm, PriorityQueue& queue)
{
	// for incremental counting
	gx->countLabelFrequencies();
	gy->countLabelFrequencies();
	gx->countBridgeFrequencies();
	gy->countBridgeFrequencies();

	// int mindist = INT_MAX;

	for(int next = 0; next < int(gx->vsize()); next++){
		mapping* m = gx->generateNextState(next, pm, bitmap_size);

		// set current mapping
		gx->setState(m, &freqx[wid]);
		gy->setState(m->size, &freqy[wid]);

		calculateDistance(m, pm);

		if(real_ged > m->distance())
			real_ged = m->distance(); 
		
		// real_ged = m->distance();
		// cout<<"real_ged = "<<real_ged<<", ";
		if(m->distance() <= threshold) saveUnmappedOf(m);
		if(m->distance() > threshold){ 
			delete m;
		}
		else queue.push(m);

		// restore parent mapping
		gx->setState(pm, &freqx[wid]);
		gy->setState(pm->size, &freqy[wid]);
	}
	// real_ged = mindist;
	// cout<<"5: real_ged = "<<real_ged<<endl;
	// cout<<endl;
}

void NassGED::expandMCSState(mapping* pm, PriorityQueue& queue)
{
	// for incremental counting
	gx->countLabelFrequencies();
	gy->countLabelFrequencies();
	gx->countBridgeFrequencies();
	gy->countBridgeFrequencies();

	int mindist = INT_MAX;
	// mapping* mp = NULL;

	for(int next = 0; next < int(gx->vsize()); next++){
		mapping* m = gx->generateNextState(next, pm, bitmap_size);

		// set current mapping
		gx->setState(m, &freqx[wid]);
		gy->setState(m->size, &freqy[wid]);

		calculateDistance(m, pm); 
		
		// queue.push(m);
		if(mindist > m->distance()){
			mindist = m->distance();
		// 	// if(mp != NULL) delete mp;
		// 	// mp = new mapping(m);
		// 	saveUnmappedOf(m);
		// 	queue.push(m);
		}
		// real_ged = m->distance();
		// cout<<"5: real_ged = "<<m->distance()<<", ";
		if(m->distance() <= threshold) saveUnmappedOf(m);
		if(m->distance() > threshold){ 
			delete m;
		}
		else queue.push(m);

		// restore parent mapping
		gx->setState(pm, &freqx[wid]);
		gy->setState(pm->size, &freqy[wid]);
	}
	// real_ged = mindist;
	// cout<<"5: real_ged = "<<real_ged<<", ";
	// if(mindist > threshold){
	// 	saveUnmappedOf(mp);
	// 	queue.push(mp);
	// }
}

NassGED::NassGED(coregraph* x, coregraph* y, unsigned threshold, loop_timer* timer, int wid)
	: threshold(threshold), timer(timer)
{
	if(x->vsize() < y->vsize()){ coregraph* t = x; x = y; y = t; }

	DataSet* dataset = DataSet::getInstance();
	memcpy(NassGED::orderingx[wid], NassGED::ordering, sizeof(vertex_t)*dataset->vmax());
	memcpy(NassGED::orderingy[wid], NassGED::ordering, sizeof(vertex_t)*dataset->vmax());

	this->wid = wid;
	full_mapping_size = (x->vsize() > y->vsize() ? x->vsize() : y->vsize());

	bitmap_size = full_mapping_size/64 + (full_mapping_size%64 ? 1 : 0);
	mtree = new maptree(mappingCompare(bitmap_size));//自定义比较函数,升序

	real_ged = INT_MAX;

	gx = new graph(*x, orderingx[wid]);
	gy = new graph(*y, orderingy[wid]);
}

NassGED::~NassGED()
{
	delete gx;
	delete gy;

	for(map_iter it = mtree->begin(); it != mtree->end(); it++){
		delete [] it->first;
		delete it->second;
	}

	delete mtree;
}

int NassGED::victimSelection(mapping* m, PriorityQueue& queue)
{
	queue.shrink_to_fit();

	while(NassGED::interrupted){
		usleep(1); // busy waiting with 1 microsecond sleep;
		wait[wid] =true; // for synchronization
	}
	wait[wid] = false;

	if(NassGED::victim == wid){
		queue.clear(wid);
		int distance = (int)m->distance();
		delete m;

		if(distance == 0) return INEXACT_ZERO; // negative zero
		return -1 * distance; // negative sign indicats it's inexact
	}

	return 1; // not a victim
}

bool NassGED::applyFilter(mapping* m)
{
	if(m->size == 0) return false;
	
	switch(m->unmapped->method){
	case NONE:
		cout << "THE CONTROL NEVER REACHES TO THIS CASE" << endl;
		break;
	case LF: // apply BF
	{
		gx->countVertexEdgeFrequencies();
		gy->countVertexEdgeFrequencies();

		int lb10 = compactBranchFilter();
		m->unmapped->dist[BF] = lb10/10;
		m->unmapped->dist[BF] += (lb10 % 10 ? 1 : 0);
		m->unmapped->method = BF;
		m->udist = m->unmapped->distance();

		break;
	}

	case BF: // apply PF
		m->unmapped->partition = new Inves(gx, gy); // should partition gx
		m->unmapped->method = PF;

	case PF:
		Inves* inves = m->unmapped->partition;

		if(inves != NULL){
			gx->countVertexEdgeFrequencies();
			gy->countVertexEdgeFrequencies();

			unsigned tau = threshold - (m->mdist + m->bdist);
			m->unmapped->dist[PF] = inves->lowerBound(tau);//incrementalLB(tau);
			m->udist = m->unmapped->distance(); // update udist
			if(inves->tightLB){
				delete inves;
				m->unmapped->partition = NULL;
			}
		}

		break;
	}

	return (m->distance() > threshold);
}

extern long long cands_final;
extern long long push_count;
int NassGED::computeGED()
{
	assert(wait[wid] == false);
	PriorityQueue queue;
	// real_ged = 0;

	unsigned lowerbound = compactBranchFilter();
	// cout<<"1: lower bound = "<<float(lowerbound)/10<<endl;
	if(float(lowerbound)/10 > (float)threshold){ 
		real_ged = lowerbound/10;
		return threshold+1;
	}

	Inves inves(gx, gy);
	lowerbound = inves.lowerBound(threshold);
	// cout<<"2: lower bound = "<<lowerbound<<endl;
	if(lowerbound > threshold){ 
		real_ged = lowerbound;
		return threshold+1;
	}

	cands_final++;//所有workload的候选集

	gy->reorder(*gx);

	mapping* m = new mapping(full_mapping_size, bitmap_size, true);
	queue.push(m);
	int num = 1;

	while(!queue.empty()){
		m = queue.pop();
		
		if(complete(m)){ 
			queue.clear(wid);
			int distance = (int)m->distance();
			delete m;

			push_count += queue.push_count;	
			return distance;
		}

		if(NassGED::interrupted){ //index building: reached to memory limit
			int inexact = victimSelection(m, queue);
			if(inexact < 0) return inexact;
		}

		gx->setState(m, &freqx[wid]);
		gy->setState(m->size, &freqy[wid]);

		if(applyFilter(m)){ 
			// cout<<"num = "<<num++<<", ";
			if(real_ged > m->distance())
				real_ged = m->distance();
			// cout<<"1: real_ged = "<<real_ged<<endl;
			delete m;  
			continue; 
		}
// cout<<endl;
		if(gy->vsize() == 0){ // case: deletions of remaining vertices in gx
			for(unsigned i = 0; i < gx->vsize(); i++)
				m->mdist = m->mdist + gx->mappedErrors(m->size + i);
			m->size += gx->vsize();
			m->bdist = m->udist = 0;

			if(m->distance() > threshold){ 
				delete m;
			}
			else queue.push(m);
		}
		else{
			expandState(m, queue);
			// cout<<"m->distance = "<<m->distance()<<endl;
			delete m;

			PriorityQueue::total_size[wid] = queue.size();
		}
	}

	queue.clear(wid);  // for reset status

	push_count += queue.push_count;	
	return threshold+1;
}

int NassGED::computeMCSGED()
{
	assert(wait[wid] == false);
	PriorityQueue queue;
	// real_ged = 0;

	unsigned lowerbound = compactBranchFilter();
	// cout<<"1: lower bound = "<<float(lowerbound)/10<<endl;
	// if(float(lowerbound)/10 > (float)threshold){ 
	// 	real_ged = lowerbound/10;
	// 	return threshold+1;
	// }

	Inves inves(gx, gy);
	lowerbound = inves.lowerBound(threshold);
	// cout<<"2: lower bound = "<<lowerbound<<endl;
	// if(lowerbound > threshold){ 
	// 	real_ged = lowerbound;
	// 	return threshold+1;
	// }

	// cands_final++;//所有workload的候选集

	gy->reorder(*gx);

	mapping* m = new mapping(full_mapping_size, bitmap_size, true);
	queue.push(m);
	int num = 1;

	while(!queue.empty()){
		m = queue.pop();
		
		if(complete(m)){ 
			queue.clear(wid);
			int distance = (int)m->distance();
			delete m;

			push_count += queue.push_count;	
			return distance;
		}

		if(NassGED::interrupted){ //index building: reached to memory limit
			int inexact = victimSelection(m, queue);
			if(inexact < 0) return inexact;
		}

		gx->setState(m, &freqx[wid]);
		gy->setState(m->size, &freqy[wid]);

		if(applyFilter(m)){ 
			// cout<<"num = "<<num++<<", ";
			if(real_ged > m->distance())
				real_ged = m->distance();
			// cout<<"1: real_ged = "<<real_ged<<endl;
			delete m;  
			continue; 
		}
// cout<<endl;
		if(gy->vsize() == 0){ // case: deletions of remaining vertices in gx
			for(unsigned i = 0; i < gx->vsize(); i++)
				m->mdist = m->mdist + gx->mappedErrors(m->size + i);
			m->size += gx->vsize();
			m->bdist = m->udist = 0;

			if(m->distance() > threshold){ 
				delete m;
			}
			else queue.push(m);
		}
		else{
			expandState(m, queue);
			// cout<<"m->distance = "<<m->distance()<<endl;
			delete m;

			PriorityQueue::total_size[wid] = queue.size();
		}
	}

	queue.clear(wid);  // for reset status

	push_count += queue.push_count;	
	return threshold+1;
}
