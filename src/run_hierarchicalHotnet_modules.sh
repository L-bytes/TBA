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
data=$PWD/output/hotnet/HotNet_input
intermediate=$PWD/output/hotnet/HotNet_intermediate
results=$PWD/output/hotnet/HotNet_results
scripts=$PWD/hierarchical-hotnet/src #Replace with directory where Hierarchical HotNet is installed

mkdir -p $intermediate
mkdir -p $results

### Set parameters
modules=($1)
thresholds=(003)
methods=(log2)
num_permutations=100
: '
### Run procedure on full dataset
for thr in ${thresholds[@]}
do
	echo $thr
	################################################################################
	#
	#   Create relevant input files
	#
	################################################################################
	#echo "Creating input files..."
	#python3.53 ../create_hotnet_input.py $data/edges_expression_"$thr".tsv $data/i2g_"$thr"_all.tsv $data/edge_list_"$thr"_all.tsv
	################################################################################
	#
	#   Construct similarity matrices.
	#
	################################################################################
	echo "Construct similarity matrices..."
	#echo $scripts/construct_similarity_matrix.py
	python3.5 $scripts/construct_similarity_matrix.py \
		-i   $data/edge_list_"$thr"_all.tsv \
		-o   $intermediate/similarity_matrix_"$thr"_all.h5 \
		-bof $intermediate/beta_"$thr"_all.txt
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
			-gsf $data/g2s_"$method"_all.tsv \
			-igf $data/i2g_"$thr"_all.tsv \
			-elf $data/edge_list_"$thr"_all.tsv \
			-ms  1000 \
			-o   $intermediate/score_bins_"$method"_"$thr"_all.tsv
		for i in `seq $num_permutations`
		do
			python3.5 $scripts/permute_scores.py \
				-i  $data/g2s_"$method"_all.tsv \
				-bf $intermediate/score_bins_"$method"_"$thr"_all.tsv \
				-s  "$i" \
				-o  $intermediate/scores_"$method"_"$thr"_all_"$i".tsv
		done
		################################################################################
		#
		#   Construct hierarchies.
		#
		################################################################################
		echo "Constructing hierarchies..."
		cp $data/g2s_"$method"_all.tsv $intermediate/scores_"$method"_"$thr"_all_0.tsv
		for i in `seq 0 $num_permutations`
		do
			python3.5  $scripts/construct_hierarchy.py \
				-smf  $intermediate/similarity_matrix_"$thr"_all.h5 \
				-igf  $data/i2g_"$thr"_all.tsv \
				-gsf  $intermediate/scores_"$method"_"$thr"_all_"$i".tsv \
				-helf $intermediate/hierarchy_edge_list_"$method"_"$thr"_all_"$i".tsv \
				-higf $intermediate/hierarchy_index_"$method"_"$thr"_all_gene_"$i".tsv
			break
		done
		break
		################################################################################
		#
		#   Process hierarchies.
		#
		################################################################################
		echo "Processing hierarchies..."
		# This example uses -lsb/--lower_size_bound 1 because it is a small toy example
		# with 25 vertices.  Use larger value (default is 10) for larger graphs.
		python3.5  $scripts/process_hierarchies.py \
			-oelf $intermediate/hierarchy_edge_list_"$method"_"$thr"_all_0.tsv \
			-oigf $intermediate/hierarchy_index_"$method"_"$thr"_all_gene_0.tsv \
			-pelf $(for i in `seq $num_permutations`; do echo " $intermediate/hierarchy_edge_list_"$method"_"$thr"_all_"$i".tsv "; done) \
			-pigf $(for i in `seq $num_permutations`; do echo " $intermediate/hierarchy_index_"$method"_"$thr"_all_gene_"$i".tsv "; done) \
			-lsb  5 \
			-cf   $results/clusters_hierarchies_"$method"_"$thr"_all.tsv \
			-pl   network_"$method"_"$thr"_all \
			-pf   $results/sizes_network_"$method"_"$thr"_all.pdf
		################################################################################
		#
		#   Perform consensus.
		#
		################################################################################
		echo "Performing consensus..."
		python3.5  $scripts/perform_consensus.py \
			-cf  $results/clusters_hierarchies_"$method"_"$thr"_all.tsv \
			-igf $data/i2g_"$thr"_all.tsv \
			-elf $data/edge_list_"$thr"_all.tsv \
			-n   network_"$method"_"$thr"_all \
			-s   g2s_"$method"_all \
			-t   1 \
			-cnf $results/consensus_nodes_"$method"_"$thr"_all.tsv \
			-cef $results/consensus_edges_"$method"_"$thr"_all.tsv
		################################################################################
		#
		#   Create output graph.
		#
		################################################################################
		echo "Creating expression graph"
		python3.5  $scripts/HotNet_graph_consensus.py \
			-igf $data/i2g_"$thr"_"$module".tsv \
			-elf $results/consensus_edges_"$method"_"$thr"_all.tsv \
			-gsf $data/g2s_"$method"_"$module".tsv \
			-pf $results/subnetworks_"$method"_"$thr"_all.png
		
		break
	done
	break
done
'
#: '
### Run procedure for each module individually, all edge thresholds and for the different score types
for module in ${modules[@]}
do
	echo $module
	if [ $module == grey ]
  then
    continue
  else
  	for thr in ${thresholds[@]}
  	do
  		echo $thr
  		################################################################################
  		#
  		#   Create relevant input files
  		#
  		################################################################################
  		#: '
  
  		#python3.5 ../create_hotnet_input.py $data/name_edges_expression_"$thr"_"$module".tsv $data/i2g_"$thr"_"$module".tsv $data/edge_list_"$thr"_"$module".tsv
  		################################################################################
  		#
  		#   Construct similarity matrices.
  		#
  		################################################################################
  
  		echo "Construct similarity matrices..."
  
  		#echo $scripts/construct_similarity_matrix.py
  		python3.5 $scripts/construct_similarity_matrix.py \
  			-i   $data/edge_list_"$thr"_"$module".tsv \
  			-o   $intermediate/similarity_matrix_"$thr"_"$module".h5 \
  			-bof $intermediate/beta_"$thr"_"$module".txt
  		
  		################################################################################
  		#
  		#   Permute data.
  		#
  		################################################################################
  
  		echo "Permuting scores..."
  		#'
  		for method in ${methods[@]}
  		do
  			echo $method
  			#: '
  			python3.5 $scripts/find_permutation_bins.py \
  				-gsf $data/g2s_"$method"_"$module".tsv \
  				-igf $data/i2g_"$thr"_"$module".tsv \
  				-elf $data/edge_list_"$thr"_"$module".tsv \
  				-ms  10 \
  				-o   $intermediate/score_bins_"$method"_"$thr"_"$module".tsv
  
  			for i in `seq $num_permutations`
  			do
  				python3.5 $scripts/permute_scores.py \
  					-i  $data/g2s_"$method"_"$module".tsv \
  					-bf $intermediate/score_bins_"$method"_"$thr"_"$module".tsv \
  					-s  "$i" \
  					-o  $intermediate/scores_"$method"_"$thr"_"$module"_"$i".tsv
  			done
  
  			################################################################################
  			#
  			#   Construct hierarchies.
  			#
  			################################################################################
  
  			echo "Constructing hierarchies..."
  			cp $data/g2s_"$method"_"$module".tsv $intermediate/scores_"$method"_"$thr"_"$module"_0.tsv
  
  			for i in `seq 0 $num_permutations`
  			do
  				python3.5  $scripts/construct_hierarchy.py \
  					-smf  $intermediate/similarity_matrix_"$thr"_"$module".h5 \
  					-igf  $data/i2g_"$thr"_"$module".tsv \
  					-gsf  $intermediate/scores_"$method"_"$thr"_"$module"_"$i".tsv \
  					-helf $intermediate/hierarchy_edge_list_"$method"_"$thr"_"$module"_"$i".tsv \
  					-higf $intermediate/hierarchy_index_"$method"_"$thr"_"$module"_gene_"$i".tsv
  			done
  
  			################################################################################
  			#
  			#   Process hierarchies.
  			#
  			################################################################################
  
  			echo "Processing hierarchies..."
  			# This example uses -lsb/--lower_size_bound 1 because it is a small toy example
  			# with 25 vertices.  Use larger value (default is 10) for larger graphs.
  			python3.5  $scripts/process_hierarchies.py \
  				-oelf $intermediate/hierarchy_edge_list_"$method"_"$thr"_"$module"_0.tsv \
  				-oigf $intermediate/hierarchy_index_"$method"_"$thr"_"$module"_gene_0.tsv \
  				-pelf $(for i in `seq $num_permutations`; do echo " $intermediate/hierarchy_edge_list_"$method"_"$thr"_"$module"_"$i".tsv "; done) \
  				-pigf $(for i in `seq $num_permutations`; do echo " $intermediate/hierarchy_index_"$method"_"$thr"_"$module"_gene_"$i".tsv "; done) \
  				-lsb  5 \
  				-cf   $results/clusters_hierarchies_"$method"_"$thr"_"$module".tsv \
  				-pl   network_"$method"_"$thr"_"$module" \
  				-pf   $results/sizes_network_"$method"_"$thr"_"$module".pdf
  
  			################################################################################
  			#
  			#   Perform consensus.
  			#
  			################################################################################
  
  			echo "Performing consensus..."
  			python3.5  $scripts/perform_consensus.py \
  				-cf  $results/clusters_hierarchies_"$method"_"$thr"_"$module".tsv \
  				-igf $data/i2g_"$thr"_"$module".tsv \
  				-elf $data/edge_list_"$thr"_"$module".tsv \
  				-n   network_"$method"_"$thr"_"$module" \
  				-s   g2s_"$method"_"$module" \
  				-t   1 \
  				-cnf $results/consensus_nodes_"$method"_"$thr"_"$module".tsv \
  				-cef $results/consensus_edges_"$method"_"$thr"_"$module".tsv
  
  			################################################################################
  			#
  			#   Create output graph.
  			#
  			################################################################################
  			#'
  			echo "Creating expression graph..."
  			python3.5  $scripts/HotNet_graph_consensus.py \
  				-igf $data/i2g_"$thr"_"$module".tsv \
  				-elf $results/consensus_edges_"$method"_"$thr"_"$module".tsv \
  				-gsf $data/g2s_"$method"_"$module".tsv \
  				-pf $results/subnetworks_"$method"_"$thr"_"$module".png
  		done
  	done
  fi
done
#'