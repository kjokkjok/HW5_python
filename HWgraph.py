
# coding: utf-8

# In[22]:


from Bio import SeqIO
import Bio
from collections import defaultdict
from graphviz import Digraph
import numpy as np

class Vertex:
    
    def __init__(self, seq):
        self.seq = seq
        self.coverage = 1
        self.in_edges = {}
        self.out_edges = {}
        
    def increase_coverage(self):
        self.coverage += 1

class Edge:
    
    def __init__(self,k1,k2):
        self.seq = k1 + k2[-1]
        self.n = 2
        self.coverage = 0
    
    def calc_coverage(self,c1,c2):
        self.coverage = (c1+c2)/2


class Graph:

    def __init__(self,k):
        self.vertices = {}
        self.k = k
        
    def add_read(self,read):
        read_lng = len(read)
        if read_lng < self.k:
            return
            
        kmer = read[:k]
        if kmer in self.vertices:
            self.vertices[kmer].increase_coverage()
        else:
            self.vertices[kmer] = Vertex(kmer)
        
        for next_kmer_indx in range(1,read_lng-k+1,1):
            next_kmer = read[next_kmer_indx:(next_kmer_indx+k)]
            if next_kmer in self.vertices:
                self.vertices[next_kmer].increase_coverage()
            else:
                self.vertices[next_kmer] = Vertex(next_kmer)
            
            new_edge = Edge(kmer,next_kmer)
            
            self.vertices[next_kmer].in_edges[kmer]  = [new_edge]
            
            self.vertices[kmer].out_edges[next_kmer] = [new_edge]

            kmer = next_kmer
    
    def calc_init_edge_coverage(self):
        
        for current_vertex in self.vertices.keys():
            for next_vertex in self.vertices[current_vertex].out_edges.keys():
                self.vertices[current_vertex].out_edges[next_vertex][0].calc_coverage(self.vertices[current_vertex].coverage,self.vertices[next_vertex].coverage)
        
    def Graph_viz(self, PathToWrite, full):
        
        self.graph = Digraph(comment='assembly')

        if full:
            for v, e in self.vertices.items():
                self.graph.node(v, label='{}'.format(v))
                
                for i, j in e.out_edges.items():
                    self.graph.edge(v, i,label='{}'.format(j[0].seq))

        else:
            for v, e in self.vertices.items():
                self.graph.node(v, label='coverage={}'.format(my_graph.vertices[v].coverage))
                
                for i, j in e.out_edges.items():
                    self.graph.edge(v, i,label='len={},cov={}'.format(len(j[0].seq), j[0].coverage))
        
        print(self.graph)
        with open (PathToWrite, 'w') as handle:
            handle.write(self.graph.source)
        
        
    def Graph_compressing(self):
        self.for_compressing = [(vertex, edge) for vertex, edge in self.vertices.items() if (len(edge.in_edges) == 1 and len(edge.out_edges) == 1)]

        for vertex, edge in self.for_compressing:
            if vertex in self.vertices and len(self.vertices) > 2:
                input_vertex = list(edge.in_edges.keys())[0]
                output_vertex = list(edge.out_edges.keys())[0]
                new_edge = Edge(edge.in_edges[input_vertex][0].seq, edge.out_edges[output_vertex][0].seq[-1])

                new_edge.n = self.vertices[output_vertex].in_edges[vertex][0].n + self.vertices[input_vertex].out_edges[vertex][0].n - 1

                new_edge.coverage = np.mean([self.vertices[input_vertex].out_edges[vertex][0].coverage, self.vertices[output_vertex].in_edges[vertex][0].coverage])

                self.vertices[input_vertex].out_edges[output_vertex] = [new_edge]
                self.vertices[output_vertex].in_edges[input_vertex] = [new_edge]                

                del self.vertices[output_vertex].in_edges[vertex]
                del self.vertices[input_vertex].out_edges[vertex]
                del self.vertices[vertex]
                
    
        

if __name__ == '__main__':
    
    dataset = 'Загрузки/hw_4_5_dataset.fasta'

    k = 16
    
    my_graph = Graph(k)
    direction = 'reverse'
    with open(dataset, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if direction == 'reverse':
                read = str(record.reverse_complement().seq)
                my_graph.add_read(read)
            else:
                read = str(record.seq)
                my_graph.add_read(read)

    my_graph.calc_init_edge_coverage()
    
    
#     for v in my_graph.vertices:
#         print('Vertex: {}, coverage: {}'.format(v,my_graph.vertices[v].coverage))
#         for e in my_graph.vertices[v].out_edges:
#             print('-> Out edge: {}'.format(e))
#         for e in my_graph.vertices[v].in_edges:
#             print('-> In edge: {}'.format(e))            
    
    my_graph.Graph_compressing()
    
    my_graph.Graph_viz("Graph.dot", full=True)
        


