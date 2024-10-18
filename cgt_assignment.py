# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 12:38:26 2024

@author: DELL
"""

import networkx as nx
import matplotlib.pyplot as plt
import random
import copy
import heapq
from itertools import combinations

#Question 1
def havel_hakimi(seq):
    while True: 
        seq.sort(reverse=True)
        if seq[0]== 0 and seq[len(seq)-1]== 0:
            return True
        
        v = seq[0]
        if(v<0):
            return False
        seq = seq[1:]
        
        if v>len(seq): 
            return False

        for i in range(v):
            seq[i]-= 1
            if seq[i]<0:
                return False

def draw_graph(sequence):
    if not havel_hakimi(sequence):
        print("\nThis is a not graphical sequence")
        return 0,0,0
    print("\nThis is a graphical sequence")
    G=nx.Graph()
    n=len(sequence)
    G.add_nodes_from(range(n))
    remaining_degrees=[(sequence[i],i) for i in range(n)]
    
    adjacency_and_cost=[]
    while remaining_degrees:
        remaining_degrees.sort(reverse=True)
        d,node_index=remaining_degrees.pop(0)
        if d==0:
            continue
        
        temp=[]
        for i in range(d):
            next_d, next_index=remaining_degrees.pop(0)
            cost=random.randint(1,7)
            G.add_edge(node_index, next_index, weight= cost)
            adjacency_and_cost.append((node_index, next_index, cost))
            temp.append((next_d-1, next_index))
            
        for item in temp:
            if item[0]>0:
                remaining_degrees.append(item)
    layout= nx.spring_layout(G)
    nx.draw(G, layout, with_labels=True, node_color='lightblue', node_size=800, font_size=10,font_color='black')
    edge_labels= nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edge_labels(G, layout, edge_labels)
    plt.title("Generated Graph for the Sequence")
    plt.show()
    return 1,G,adjacency_and_cost

#Question 2
def Euler(G, sequence, adjacency_matrix,vertices,edges,visited):
    if len(list(nx.connected_components(G))) > 1:
        print("\nNot a Euler Graph as the graph is disconnected")
        return
    count=0
    for i in range(vertices):
        if(sequence[i]%2!=0) :
            odd_node=i
            count+=1
    if not (count==0 or count==2) :
        print("\nNot a Euler Graph as number of vertices with odd degrees is greater than 2")
        return
    else:
        if count==2:
            print("\nEuler path exists as there are exactly 2 vertices of odd degree:")
            current_vertex=odd_node
        if count==0:
            print("\nEuler circuit exists as all vertices are of even degree:")
            current_vertex=0
        edge_count=0
        H=G.copy()
        while(edge_count<edges):
            for i in range(vertices):
                if adjacency_matrix[current_vertex][i]==1 and visited[current_vertex][i]==0:
                    H.remove_edge(current_vertex,i)
                    if(len(H.nodes)==2):
                        print(current_vertex,"->",i)
                        edge_count+=1
                    if(dict(H.degree)[current_vertex]==0):
                        H.remove_node(current_vertex)
                    if(dict(H.degree)[i]==0 and (G.number_of_nodes())>1):
                        H.add_edge(current_vertex,i)
                        continue
                    if len(list(nx.connected_components(H))) > 1:
                        H.add_edge(current_vertex,i)
                        continue
                    print(current_vertex,"->",i)
                    visited[i][current_vertex]=1
                    visited[current_vertex][i]=1
                    current_vertex=i
                    edge_count+=1
                    break
    
#Question 3
def dijkstra(graph, src, vertices):
    dist = {i: float('inf') for i in range(vertices)}
    dist[src] = 0
    pq = [(0, src)]
    visited = []

    while pq:
        current_dist, current_vertex = heapq.heappop(pq)
        if current_vertex in visited:
            continue
        visited.append(current_vertex)
        for neighbor, attributes in graph[current_vertex].items():
            weight = attributes['weight']
            distance = current_dist + weight
            if distance < dist[neighbor]:
                dist[neighbor] = distance
                heapq.heappush(pq, (distance, neighbor))
    return dist
    
#Question 4
#Kruskal Algorithm
#For minimal spanning tree
def kruskal_mst(G, vertices):
    edges = []
    for u, v, data in G.edges(data=True):
        edges.append((data['weight'], u, v))
    edges.sort()
    
    # Union-Find data structure to detect cycles
    parent = list(range(vertices))
    rank = [0] * vertices
    
    def find(u):
        if parent[u] != u:
            parent[u] = find(parent[u])
        return parent[u]
    
    def union(u, v):
        root_u = find(u)
        root_v = find(v)
        if root_u != root_v:
            if rank[root_u] > rank[root_v]:
                parent[root_v] = root_u
            elif rank[root_u] < rank[root_v]:
                parent[root_u] = root_v
            else:
                parent[root_v] = root_u
                rank[root_u] += 1
    
    # Kruskal's algorithm
    mst = []
    plot_mst=nx.Graph()
    for weight, u, v in edges:
        if find(u) != find(v):
            plot_mst.add_edge(u,v,weight=weight)
            mst.append((u, v, weight))
            union(u, v)
    return mst, plot_mst

#For fundamental cutsets
def find_fundamental_cutsets(G, mst, vertices):
    cutsets = []
    mst_graph = nx.Graph()
    mst_graph.add_weighted_edges_from(mst)

    for edge in mst:
        u, v, _ = edge
        mst_graph.remove_edge(u, v)
        components = list(nx.connected_components(mst_graph))
        mst_graph.add_edge(u, v)
        
        cutset = []
        for node_u in components[0]:
            for node_v in components[1]:
                if G.has_edge(node_u, node_v):
                    cutset.append((node_u, node_v))
        cutsets.append(cutset)
    return cutsets

#For fundamental circuits
def find_fundamental_circuits(G, mst, vertices):
    mst_graph = nx.Graph()
    mst_graph.add_weighted_edges_from(mst)

    circuits = []
    for u, v in G.edges():
        if (u, v) not in mst_graph.edges() and (v, u) not in mst_graph.edges():
            mst_graph.add_edge(u, v)
            cycle = list(nx.find_cycle(mst_graph))
            circuits.append(cycle)
            mst_graph.remove_edge(u, v)

    return circuits

#Question 5
#To find edge-connectivity
def edge_connectivity(G,edges):
    for i in range(1,edges):
        G_copy=copy.deepcopy(G)
        removed_count=0
        for p in combinations(G_copy.edges(), i):
            G_copy=copy.deepcopy(G)
            removed_count=0
            for u,v in p:
                G_copy.remove_edge(u,v)
                removed_count+=1
                if (len(list(nx.connected_components(G_copy)))>1):
                    return removed_count,p

#To find vertex connectivity
def vertex_connectivity(G,vertices):
    for i in range(1,vertices):
        G_copy=copy.deepcopy(G)
        removed_count=0
        for p in combinations(G_copy.nodes(), i):
            G_copy=copy.deepcopy(G)
            for j in p:
                G_copy.remove_node(j)
                removed_count+=1
                if (len(list(nx.connected_components(G_copy)))>1) or len(G_copy.nodes())==1:
                    return removed_count,p
                
#Main program code  
degree_sequence=[]
n=int(input("Enter number of elements in graphical sequence : "))
print("Enter elements : ")
degree_sum=0
for i in range(n):
    element=int(input())
    degree_sequence.append(element)
    degree_sum+=element
vertices= len(degree_sequence)
edges= int(degree_sum/2)

#For Havel-Hakimi algorithm (Q1)
print("\nQuestion 1")
graphical, G, adjacency_and_cost= draw_graph(degree_sequence)

u=[]
adjacency_matrix=[]
visited=[]
cost_adjacency=[]

#For Euler Graph (Q2) , if and only if it is graphical sequence
if graphical :
    for i in range(vertices):
        u=[0]*vertices
        adjacency_matrix.append(u)
    visited=copy.deepcopy(adjacency_matrix)
    cost_adjacency=copy.deepcopy(adjacency_matrix)
    for item in adjacency_and_cost:
        adjacency_matrix[item[0]][item[1]]=1
        adjacency_matrix[item[1]][item[0]]=1
        cost_adjacency[item[0]][item[1]]=item[2]
        cost_adjacency[item[1]][item[0]]=item[2]
        
    print("\nAdjacency matrix:")
    for i in range(vertices):
        print(adjacency_matrix[i])
    print("\nQuestion 2")
    Euler(G,degree_sequence,adjacency_matrix,vertices,edges,visited)
    print("\nCost of edges as adjacency matrix:")
    for i in range(vertices):
        print(cost_adjacency[i])
    
    #For Dijkstra Algorithm (Q3)
    print("\nQueston 3")
    src = int(input("\nEnter the source vertex for Dijkstra's algorithm: "))
    distances = dijkstra(G, src, vertices)
    print("Shortest distances from vertex" , src,":")
    for vertex, distance in distances.items():
        print("Vertex",vertex, "-> Distance: ",distance)
        
    #For Kruskal Algorithm (Q4)
    print("\nQuestion 4")
    if len(list(nx.connected_components(G))) > 1:
        print("\nAs graph is disconnected, it doesn't have minimum spanning tree.")
    else:
        mst, plot_mst = kruskal_mst(G, vertices)
        print("\nMinimum Spanning Tree (MST):")
        total_cost=0
        for edge in mst:
            print("Edge : (", edge[0],",",edge[1],") | Cost :",edge[2])
            total_cost+=edge[2]
        print("Total cost : ",total_cost)
        layout_mst= nx.spring_layout(plot_mst)
        nx.draw(plot_mst, layout_mst, with_labels=True, node_color='lightblue', node_size=300, font_size=10,font_color='black')
        edge_labels_mst= nx.get_edge_attributes(plot_mst, 'weight')
        nx.draw_networkx_edge_labels(plot_mst, layout_mst, edge_labels_mst)
        plt.title("Minimum Spanning Tree from Kruskal's Algorithm")
        plt.show()
    
        #For Fundamental Cutsets
        cutsets = find_fundamental_cutsets(G, mst, vertices)
        print("\nFundamental Cutsets:")
        for i, cutset in enumerate(cutsets):
            print("Cutset", i + 1, ":", cutset)

        #For Fundamental Circuits
        circuits = find_fundamental_circuits(G, mst, vertices)
        print("\nFundamental Circuits:")
        for i, circuit in enumerate(circuits):
            print("Circuit", i + 1,":", circuit)
        
    #For Edge-Connectivity and Vertex-Connectivity and K-Connectivity(Q5)
    print(list(combinations(G.nodes(),2)))
    print("\nQuestion 5")
    if len(list(nx.connected_components(G))) >1 :
        print("Edge connectivity of the graph is 0")
        print("Vertex connectivity of the graph is 0")
        print("K-connectivity of the graph is 0")
    else:
        edge_conn, removed_edges= edge_connectivity(G,edges)
        print("\nEdge connectivity of the graph is", edge_conn)
        print("Edges to remove for edge connectivity are :",removed_edges)
        vertex_conn, removed_vertices= vertex_connectivity(G,vertices)
        print("\nVertex connectivity of the graph is", vertex_conn)
        print("Vertices to remove for vertex connectivity are :",removed_vertices)
        #K-connectivity : Largest number of vertices that can be removed without disconnecting the graph
        #So K-connectivity of this sequence should be equal to its vertex connectivity
        print("\nK-connectivity of graph is :",vertex_conn)
        
else:
    print("Question 2, 3, 4, 5:")
    print("As sequence is not graphic, euler graph, mst, fundamental cut-sets and circuits, graph connectivity can't be evaluated")        