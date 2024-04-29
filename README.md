# SMCS
This code of our paper "SMCS: Accelerating Graph Similarity Search with Maximum Common Subgraph".

# Datasets
We provide AIDS and PubChem datasets in [data]. (data/AIDS and data/CHEM)

# How to build binaries
> make

# How to run(search)
> nass threshold data_file tg_file mcs_file mcs_index_file [index_file1 index_file2 query_file]

# How to run an index
> nass-index threshold data_file index_file thresholdGroup_file resort_graph_file mcs_priginal_file mcs_file heuristic

nass-index has the following options.<br>
-M memory limit in MB (default 1000, i.e., 1GB)<br>
-p number of threads (default 8)<br>
--coordinator (to play a role of a coordinator for distributed index building, see below)<br>
-s sampling_rate (experimental, see the description at the end) <br>
