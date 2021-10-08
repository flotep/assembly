#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib 
import matplotlib.pyplot as plt
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

import itertools as it
from collections import Counter

__author__ = "Florian TEP"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Florian TEP"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Florian TEP"
__email__ = "tepflorian@gmail.com"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    content = []
    with open (fastq_file, "r") as file:
        for line in file : 
            new_content = line.strip()
            content.append(new_content)
    for i in range(1,(len(content)),4): 
        yield content[i]
    
   
            


def cut_kmer(read, kmer_size):
    for i in range (0,len(read)) :
        if len(read)-i >= kmer_size-1 : 
            yield read[i:i+kmer_size]

    
def build_kmer_dict(fastq_file, kmer_size):
    kmers_list =[]
    sequence = read_fastq(fastq_file)
    for i in it.chain(sequence):
        kmers=cut_kmer(i, kmer_size)
        for kmer in it.chain(kmers):
            if len(kmer) == kmer_size : #liste le kmer que si sa longueur est égale au kmer_size
                kmers_list.append(kmer)
    #print (kmers_list)

    counter = Counter(kmers_list)
   
    return dict(counter)

def build_graph(kmer_dict):
    G = nx.DiGraph()
    for item in kmer_dict.items() :
        kmer = item[0]
        weight = item[1]
        prefixe = kmer[:-1]
        suffixe = kmer[1:]
        G.add_edge(prefixe, suffixe, weight=weight)
    return (G)

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    pass

def std(data):
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):
    pass

def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    pass

def get_sink_nodes(graph):
    pass

def get_contigs(graph, starting_nodes, ending_nodes):
    pass

def save_contigs(contigs_list, output_file):
    pass


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    #with open(graph_file, "wt") as save:
            #pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    #sequences = list(read_fastq(args.fastq_file))
    sequences = list(read_fastq(args.fastq_file))
    for i in sequences :
        print (i)

    dict = build_kmer_dict(args.fastq_file,args.kmer_size)
    print (dict)

    graph = build_graph(dict)
    draw_graph(graph, "G.pdf") 

    #G.add_nodes_from(dict)

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()

