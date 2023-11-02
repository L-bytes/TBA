#!/bin/bash

# Hierarchical HotNet is parallelizable, but this script runs each Hierarchical
# HotNet sequentially.  Please see the example_commands_parallel.sh script for a
# parallelized example.

# Compile Fortran module.
#cd ../src
#f2py -c fortran_module.f95 -m fortran_module > /dev/null
#cd ..

################################################################################
#
#   Prepare data.
#
################################################################################

### Create data, intermediate data and results, and results directories.
echo "Creating directory topology"
data=$PWD/output/$3/$4/output/hotnet/HotNet_input
intermediate=$PWD/output/$3/$4/output/hotnet/HotNet_intermediate
results=$PWD/output/$3/$4/output/hotnet/HotNet_results
scripts=$PWD/hierarchical-hotnet/src #Replace with directory where Hierarchical HotNet is installed

mkdir -p $intermediate
mkdir -p $results

### Set parameters
modules=($1)
thresholds=($2)
methods=(log2)
num_permutations=100
num_cores=7

### Run procedure for each module individually, all edge thresholds and for the different score types
for module in ${modules[@]}
do
	echo $module
	start=$(date +%s)
	if [ $module == grey ]
  then
    continue
  else
  	for thr in ${thresholds[@]}
  	do
  		echo $thr
  		################################################################################
  		#
  		#   Construct similarity matrices.
  		#
  		################################################################################
  
  		echo "Construct similarity matrices..."
  
  		#echo $scripts/construct_similarity_matrix.py
  		python3.5 $scripts/construct_similarity_matrix.py \
  			-i   $data/edge_list_"$module".tsv \
  			-o   $intermediate/similarity_matrix_"$module".h5 \
  			-bof $intermediate/beta_"$module".txt
  		
  		################################################################################
  		#
  		#   Permute data.
  		#
  		################################################################################
  
  		echo "Permuting scores..."
  		for method in ${methods[@]}
  		do
  			echo $method
  			python3.5 $scripts/find_permutation_bins.py \
  				-gsf $data/g2s_"$method"_"$module".tsv \
  				-igf $data/i2g_"$module".tsv \
  				-elf $data/edge_list_"$module".tsv \
  				-ms  5 \
  				-o   $intermediate/score_bins_"$method"_"$module".tsv
  
  			parallel -j $num_cores \
  				python3.5 $scripts/permute_scores.py \
  					-i  $data/g2s_"$method"_"$module".tsv \
  					-bf $intermediate/score_bins_"$method"_"$module".tsv \
  					-s  {} \
  					-o  $intermediate/scores_"$method"_"$module"_{}.tsv \
  			  ::: `seq $num_permutations`
  
  			################################################################################
  			#
  			#   Construct hierarchies.
  			#
  			################################################################################
  
  			echo "Constructing hierarchies..."
  			cp $data/g2s_"$method"_"$module".tsv $intermediate/scores_"$method"_"$module"_0.tsv
  
  			parallel -j $num_cores \
  				python3.5  $scripts/construct_hierarchy.py \
  					-smf  $intermediate/similarity_matrix_"$module".h5 \
  					-igf  $data/i2g_"$module".tsv \
  					-gsf  $intermediate/scores_"$method"_"$module"_{}.tsv \
  					-helf $intermediate/hierarchy_edge_list_"$method"_"$module"_{}.tsv \
  					-higf $intermediate/hierarchy_index_"$method"_"$module"_gene_{}.tsv \
          ::: `seq 0 $num_permutations`
  
  			################################################################################
  			#
  			#   Process hierarchies.
  			#
  			################################################################################
  
  			echo "Processing hierarchies..."
  			python3.5  $scripts/process_hierarchies.py \
  				-oelf $intermediate/hierarchy_edge_list_"$method"_"$module"_0.tsv \
  				-oigf $intermediate/hierarchy_index_"$method"_"$module"_gene_0.tsv \
  				-pelf $(for i in `seq $num_permutations`; do echo " $intermediate/hierarchy_edge_list_"$method"_"$module"_"$i".tsv "; done) \
  				-pigf $(for i in `seq $num_permutations`; do echo " $intermediate/hierarchy_index_"$method"_"$module"_gene_"$i".tsv "; done) \
  				-lsb   10 \
  				-cf   $results/clusters_hierarchies_"$method"_"$module".tsv \
  				-pl   network_"$method"_"$module" \
  				-pf   $results/sizes_network_"$method"_"$module".pdf \
          -nc   $num_cores
  
  			################################################################################
  			#
  			#   Perform consensus.
  			#
  			################################################################################
  
  			echo "Performing consensus..."
  			python3.5  $scripts/perform_consensus.py \
  				-cf  $results/clusters_hierarchies_"$method"_"$module".tsv \
  				-igf $data/i2g_"$module".tsv \
  				-elf $data/edge_list_"$module".tsv \
  				-n   network_"$method"_"$module" \
  				-s   g2s_"$method"_"$module" \
  				-t   1 \
  				-cnf $results/consensus_nodes_"$method"_"$module".tsv \
  				-cef $results/consensus_edges_"$method"_"$module".tsv 
  
  			################################################################################
  			#
  			#   Create output graph.
  			#
  			################################################################################
        
  			echo "Creating expression graph..."
  			python3.5  $scripts/HotNet_graph_consensus.py \
  				-igf $data/i2g_"$module".tsv \
  				-elf $results/consensus_edges_"$method"_"$module".tsv \
  				-gsf $data/g2s_"$method"_"$module".tsv \
  				-pf $results/subnetworks_"$method"_"$module".png
      done
  	done
  fi
  end=$(date +%s)
  echo "Elapsed Time: $(($end-$start)) seconds"
done