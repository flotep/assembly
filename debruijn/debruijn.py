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

from networkx.algorithms.shortest_paths.unweighted import predecessor
from networkx.utils.misc import graphs_equal
random.seed(9001)
from random import randint


import itertools as it
from collections import Counter
import statistics

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
            if len(kmer) == kmer_size : 
                kmers_list.append(kmer)

    #Répertorie dans un dict chaque k-mer avec leur nombre d'occurences (à partir de la liste de k-mers)
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
    for path in path_list : 

        if delete_entry_node == True and delete_sink_node == True: 
            graph.remove_nodes_from(path) #On supprime tous les noeuds

        elif delete_entry_node == True and delete_sink_node == False: 
            keep_sink = [node for node in path[:-1]]
            graph.remove_nodes_from(keep_sink) #On supprime tous les noeuds sauf le dernier

        elif delete_entry_node == False and delete_sink_node == True:
            keep_entry = [node for node in path[1:]]
            graph.remove_nodes_from(keep_entry) #On supprime tous les noeuds sauf le premier

        elif delete_entry_node == False and delete_sink_node == False:
            keep_both= [node for node in path[1:-1]]
            graph.remove_nodes_from(keep_both) #On supprime tous les noeuds sauf le premier et le dernier

    return graph

def std(data):
    return statistics.stdev(data)

  
def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):

    #On regarde si les fréquences des chemins sont différentes
    if std(weight_avg_list) > 0 : 
        max_weight = max(weight_avg_list)
        max_weight_index = weight_avg_list.index(max_weight)
        best_path_index = max_weight_index #Le meilleur chemin est celui le plus fréquent

    #Sinon on regarde si les longueurs des chemins sont différentes
    else : 
        if std(path_length) > 0 : 
            max_length = max(path_length)
            max_length_index = path_length.index(max_length)
            best_path_index = max_length_index #Le meilleur chemin est celui le plus long
        
        #Sinon on tire au sort aléatoirement le meilleur chemin 
        else : 
            random_index = randint(0,len(path_list)-1)
            best_path_index = random_index

    #On retire le meilleur chemin de la liste des chemins à supprimer
    del path_list[best_path_index]

    #On supprime les autres chemins
    graph = remove_paths(graph, path_list,delete_entry_node, delete_sink_node)
    return graph


def path_average_weight(graph, path):
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])
   

def solve_bubble(graph, ancestor_node, descendant_node):
    weight_avg_list = []
    path_length =[]

    #On liste les chemins simples possibles entre le noeud ancêtre et descendant
    path_list = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    
    #Pour chaque chemin possible
    for path in path_list :
        #On calcule le poid moyen des noeuds du chemin 
        weight_avg = path_average_weight(graph,path)
        weight_avg_list.append(weight_avg)

        #On calcule le nombre de noeuds dans le chemin
        length = len(path)
        path_length.append(length) 
                    
    #On garde seulement le meilleur chemin dans le graphe
    graph = select_best_path(graph, path_list, path_length, weight_avg_list)
    return graph
    
def simplify_bubbles(graph):
    bubble = False
    for node in graph : 
        if graph.has_node(node):
            predecessors_list = list(graph.predecessors(node))
            if len(predecessors_list) > 1:
                combinaisons = list(it.combinations(predecessors_list, 2))
                for combinaison in it.chain(combinaisons) :
                    #print (node)
                    #print (combinaison)
                    node_ancestor = nx.lowest_common_ancestor(graph, combinaison[0], combinaison[1])
                    #print (node_ancestor)
                    if node_ancestor != None:
                        bubble = True
                        #print(bubble)
                        break
               
    if bubble == True :
        graph = simplify_bubbles(solve_bubble(graph, node_ancestor, node))
    return graph

def solve_entry_tips(graph, starting_nodes):
    entry_path_list =[]

    #Pour chaque noeud du graphe avec 2 prédecesseurs minimum
    for node in graph : 
        predecessors_list = list(graph.predecessors(node))
        if len(predecessors_list) > 1:

            #On récupère tous les chemins simples possibles entre le noeud et les noeuds d'entrées
            for starting_node in starting_nodes :
                entry_path_list = entry_path_list + list(nx.all_simple_paths(graph,starting_node, node))

    #S'il existe plus d'un chemin d'entrée possible : 
    if len(entry_path_list) > 1:
        weight_avg_list = []
        path_length = []
        for path in entry_path_list :
            weight_avg = path_average_weight(graph,path)
            weight_avg_list.append(weight_avg)
            length = len(path)
            path_length.append(length) 
        
        #On sélectionne le meilleur chemin (le + fréquent, le + long ou random)
        #On supprime les autres chemins en gardant le dernier noeud (car on retire juste les pointes)
        graph = select_best_path(graph, entry_path_list, path_length, weight_avg_list,delete_entry_node=True, delete_sink_node=False)
    return graph


def solve_out_tips(graph, ending_nodes):
    ending_path_list =[]

    #Pour chaque noeud du graphe avec 2 successeurs minimum
    for node in graph : 
        successors_list = list(graph.successors(node))
        if len(successors_list) > 1:

            #On récupère tous les chemins simples possibles entre le noeud et les noeuds de sorties
            for ending_node in ending_nodes :
                ending_path_list = ending_path_list + list(nx.all_simple_paths(graph,node, ending_node))

    #S'il existe plus d'un chemin de sortie possible : 
    if len(ending_path_list) > 1:
        weight_avg_list = []
        path_length = []
        for path in ending_path_list :
            weight_avg = path_average_weight(graph,path)
            weight_avg_list.append(weight_avg)
            length = len(path)
            path_length.append(length) 
        
        #On sélectionne le meilleur chemin (le + fréquent, le + long ou random)
        #On supprime les autres chemins en gardant le premier noeud (car on retire juste les pointes)
        graph = select_best_path(graph, ending_path_list, path_length, weight_avg_list,delete_entry_node=False, delete_sink_node=True)
    return graph


def get_starting_nodes(graph):
    starting_nodes = []
    for node in graph :
        predecessors = list(graph.predecessors(node))
        if len(predecessors)==0:
            starting_nodes.append(node)
    return starting_nodes

def get_sink_nodes(graph):
    sink_nodes = []
    for node in graph :
        successors = list(graph.successors(node))
        if len(successors)==0:
            sink_nodes.append(node)
    return sink_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs=()

    for source in starting_nodes : 
        for end in ending_nodes : 
            #S'il y existe un chemin entre source et end dans le graphe
            if nx.has_path (graph, source,end) == True : 
                paths = list(nx.all_simple_paths(graph, source,end)) #Génère les chemins simples entre ces 2 noeuds

                #Pour chaque chemin possible  
                for path in it.chain(paths): 
                    contig = source
                    path_length = len(source)+(len(path)-1) #Longueur du contig = taille 1er kmer(source) + 1 nt par noeud supplémentaire
                    
                    #Pour chaque noeud du chemin 
                    for i in range (1,len(path)) : 
                        nucleotide = path[i][-1] 
                        contig = contig + nucleotide #Chaque noeud du chemin donne le nt suivant (= le dernier nt du noeud)

                    contigs = contigs + ((contig,path_length),) #Tuples contenant la composition et la longueur de chaque chemin/contig
    return (contigs)


def save_contigs(contigs_list, output_file):
    output =  open (output_file, "w")
    id = 0
    for contig in contigs_list :
        output.write(">contig_{} len={}".format(id,contig[1]) + '\n') #contig[1] dans le tuple = longueur du contig 
        output.write(fill(contig[0]) + '\n')
        id += 1
    output.close()

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
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
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
    #Récupère les arguments
    args = get_arguments()

    #Lecture du fichier de séquences et construction du dictionnaire recensant les k-mers 
    kmer_dict = build_kmer_dict(args.fastq_file,args.kmer_size)

    #Construction du graphe à partir du dictionnaire de k-mers
    graph = build_graph(kmer_dict)
    print(graph)

    #Résolution des bulles
    graph = simplify_bubbles(graph)
    print(graph)

    #Résolution des pointes d'entrée
    starting_nodes = get_starting_nodes(graph)
    graph = solve_entry_tips(graph,starting_nodes)
    print(graph)

    #Résolution des pointes de de sortie
    sink_nodes = get_sink_nodes(graph)
    graph = solve_out_tips(graph,sink_nodes)
    print(graph)

    #Ecriture du/des contig(s) dans un fichier au format fasta
    contigs = get_contigs(graph, starting_nodes, sink_nodes)
    save_contigs(contigs, args.output_file)

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

