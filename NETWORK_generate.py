#!/usr/bin/python
# -*- coding: UTF-8 -*-
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from netFunction import *
def NW(n,k,p,folder):
    plt.figure(figsize=(7.5, 2.25))
    # Generate the graph
    # G = nx.newman_watts_strogatz_graph(n, k, p)
    G = nx.random_graphs.barabasi_albert_graph(n, k,seed=n+k)
    # G = nx.watts_strogatz_graph(n, k, p,seed=n+k)
    A = np.array(nx.adjacency_matrix(G).todense())
    ebc = nx.edge_betweenness_centrality(G)
    print(nx.info(G))
    print(matrix_file(A,folder+"Adjacencymatrix.txt"))
    print(Net_algorthm(folder,'nx.degree_centrality', G))
    print(Net_algorthm(folder,'nx.closeness_centrality',G))
    print(Net_algorthm(folder,'nx.betweenness_centrality',G))
    print(Net_algorthm(folder, 'nx.eigenvector_centrality', G))
    print(edge_betweenness_centrality2file(ebc,folder+'edge_betweenness_centrality2file.txt'))
    edge=nx.edges(G)
    print(edge2file(edge,folder+"edgelist.txt"))
    Du=list(G.degree(range(n)))
    print(degreefile(Du, folder+"node_degree.txt"))
    # Create layout and draw
    # plt.subplot(1, 3, 1)

    #
    # pos = nx.kamada_kawai_layout(G)
    # nx.draw_networkx(G, pos=pos,with_labels=False,node_size=5)
    # plt.title("NW_networks")    # plt.show()
if __name__ == '__main__':
    NW(1000,8,0.8,'F:\\py\\Mtb\\16_BA1000_8\\')

    # NW(500,7,0.34,'F:\\py\\Mtb\\8_NW500_7_034\\')8_SW1000_8_10
    # NW(50,6,0.5,'F:\\py\\Mtb\\6_SW50_6_05\\')
    # NETWORK_NW(500, 3, 1)  # 随机
    # NETWORK_BA(500, 2)
    # NETWORK_SW(500, 4, 0.8)
    # NETWORK_NW(500, 4, 0.02)  # 规则
    # NETWORK_BA(5000,4)
    # NETWORK_SW(5000,8,0.5)
    # NETWORK_NW(5000,7,0.34)






