#!/usr/bin/python3.5

# Load modules.
import math
import sys, argparse
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

from hhio import save_index_gene, save_edge_list, save_gene_score

# Parse arguments.
def get_parser():
	description = 'Generate example vertex-weighted graphs.'
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-igf', '--index_gene_file', type=str, required=True, help='Index-gene filename')
	parser.add_argument('-elf', '--edge_list_file', type=str, required=True, help='Edge list filename')
	parser.add_argument('-gsf', '--gene_score_file', type=str, required=True, help='Gene-score filename')
	parser.add_argument('-pf', '--plot_file', type=str, required=False, help='Plot filename')
	return parser
"""
# Run script.
def run(args):
	# Initialize nodes, etc.
	node_list = list()
	edge_list = list()
	score_list_a = list()
	score_list_b = list()
	pos = dict()

	# Add a clique.
	letters = 'abcdefgh'
	node_list += [letters[i] for i in range(8)]
	for i in range(8):
		for j in range(i+1, 8):
			edge_list += [(letters[i], letters[j])]
	score_list_a += [5]*8
	score_list_b += [0.1]*8
	for i in range(8):
		pos[letters[i]] = (math.cos(2*math.pi*i/8)-2, math.sin(2*math.pi*i/8))


	# Add an approximate clique.
	letters = 'ijklmnopij'
	node_list += [letters[i] for i in range(8)]
	for i in range(8):
		edge_list += [(letters[i], letters[i+1])]
		edge_list += [(letters[i], letters[i+2])]
	score_list_a += [10]*8
	score_list_b += [0.7, 0.4, 0.5, 0.6, 0.6, 0.6, 0.5, 0.4]
	for i in range(8):
		pos[letters[i]] = (math.cos(2*math.pi*i/8)+2, math.sin(2*math.pi*i/8))

	# Add a cycle.
	letters = 'qrstuvwxq'
	node_list += [letters[i] for i in range(8)]
	for i in range(8):
		edge_list += [(letters[i], letters[i+1])]
	score_list_a += [3]*8
	score_list_b += [0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2]
	for i in range(8):
		pos[letters[i]] = (math.cos(2*math.pi*i/8), math.sin(2*math.pi*i/8)-2)

	# Add a point connecting the clique, approximate clique, and cycle, add an isolated point.
	node_list += ['y']
	edge_list += [('a', 'y'), ('m', 'y'), ('s', 'y')]
	score_list_a += [5]
	score_list_b += [1.0]
	pos['y'] = (0.0, 0.0)

	# Save Hierarchical HotNet input.
	index_to_gene = dict((i+1, gene) for i, gene in enumerate(node_list))
	gene_to_index = dict((gene, i) for i, gene in index_to_gene.items())
	gene_to_score_a = dict((gene, score) for gene, score in zip(node_list, score_list_a))
	gene_to_score_b = dict((gene, score) for gene, score in zip(node_list, score_list_b))
"""
def read_indeces(index_file):
	""" Read nodes from index file """
	indeces = {}
	indeces_inv = {}
	with open(index_file, 'r') as f:
		lines = f.readlines()
		for l in lines:
			l = l.strip().split('\t')
			indeces[l[0]] = l[1]
			indeces_inv[l[1]] = l[0]
	return indeces, indeces_inv
	
def read_edges(edge_list_file, indeces):
	""" Read edges from index file """
	nodes = []
	edges = []
	with open(edge_list_file, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			i1 = line[0]
			i2 = line[1]
			edges.append((i1, i2))
			if i1 not in nodes:
				nodes.append(i1)
			if i2 not in nodes:
				nodes.append(i2)
	return nodes, edges

def read_scores(nodes, score_file):
	""" Read scores from score file """
	scores = []
	with open(score_file, 'r') as f:
		n = []
		s = []
		for line in f:
			line = line.strip().split('\t')
			n.append(line[0])
			s.append(float(line[1]))
	for node in nodes:
		if node in n:
			scores.append(s[n.index(node)])
		else:
			scores.append(0)
	
	return scores
	
def create_pos(nodes):
	""" Create circular positions """
	l = len(nodes)
	pos = {}
	
	for i in range(l):
		pos[nodes[i]] = (math.cos(2*math.pi*i/l), math.sin(2*math.pi*i/l))
	return pos

def main():
	args = get_parser().parse_args(sys.argv[1:])
	
	# Draw the graph for illustration.
	indeces, indeces_inv = read_indeces(args.index_gene_file)
	node_list, edge_list = read_edges(args.edge_list_file, indeces)
	for x in edge_list:
		if 'RNF14' in x:
			print (x)
	score_list = read_scores(node_list, args.gene_score_file)
	pos = create_pos(node_list)
	
	G = nx.Graph()
	G.add_nodes_from(node_list)
	G.add_edges_from(edge_list)

	plt.figure(figsize=(10, 10))
	#plt.subplot(121)
	node_labels = dict((v, r'${}$'.format(v)) for v in node_list)
	#nx.draw_networkx_nodes(G, pos=pos, nodelist=node_list, node_color=score_list, alpha=0.4, cmap=plt.cm.YlOrRd)
	nx.draw_networkx_nodes(G, pos=pos, nodelist=node_list, node_color=score_list, alpha=0.9, cmap=plt.cm.bwr, vmin=-2, vmax=2)
	nx.draw_networkx_edges(G, pos=pos, alpha=0.3)
	nx.draw_networkx_labels(G, pos=pos, labels=node_labels, alpha=1.0)
	plt.axis('off')
	"""
	plt.subplot(122)
	node_labels = dict((v, r'${}$'.format(v)) for v in node_list)
	nx.draw_networkx_nodes(G, pos=pos, nodelist=node_list, node_color=score_list_b, alpha=0.2, cmap=plt.cm.YlOrRd)
	nx.draw_networkx_edges(G, pos=pos, alpha=0.2)
	nx.draw_networkx_labels(G, pos=pos, labels=node_labels, alpha=0.0)
	plt.axis('off')
	"""
	if args.plot_file:
		plt.savefig(args.plot_file, bbox_inches='tight')
	#plt.show()

if __name__ == "__main__":
	main()
