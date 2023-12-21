import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
from collections import Counter

class AberrationMultigraph:
    def __init__(self, chromatin_edges, dsb_edges, rejoin_edges, name=''):
        self.graph = nx.Graph()
        self.chromatin_edges = tuple(sorted(tuple(sorted(edge)) for edge in chromatin_edges))
        self.dsb_edges = tuple(sorted(tuple(sorted(edge)) for edge in dsb_edges))
        self.rejoin_edges = tuple(sorted(tuple(sorted(edge)) for edge in rejoin_edges))
        self.name = name
        self.graph.add_edges_from(self.chromatin_edges, color='chromatin')
        self.graph.add_edges_from(self.dsb_edges, color='dsb')
        self.graph.add_edges_from(self.rejoin_edges, color='misrejoining')
        chromosome = 0
        dsb_ends = 0
        free_ends = set(u for (u,v) in self.dsb_edges).union(set(v for (u,v) in self.dsb_edges))
        for i in self.graph.nodes:
            if i not in free_ends:
                dsb_ends += 1
            self.graph.nodes[i]['chromosome'] = chromosome
            if dsb_ends == 2:
                chromosome += 1
                dsb_ends = 0
        self.chromosome = chromosome

    def diameter(self):
        return np.inf if not self.is_connected() else nx.distance_measures.diameter(self.graph)
    
    def cycles(self):
        return nx.cycles.cycle_basis(nx.Graph([(u,v) for u,v,color in self.graph.edges.data('color') if color!='chromatin']))
    
    def cycle_structure(self):
        return Counter([len(cycle) for cycle in self.cycles()])
    
    def is_connected(self):
        return nx.is_connected(self.graph)
    
    def display_config(self, view):
        next_dict = {u: v for u,v in view.edges}
        prev_dict = {v: u for u,v in view.edges}
        print(next_dict)
        print(prev_dict)
        components = []
        vertices = set(v for v in self.graph.nodes if (v in next_dict) ^ (v in prev_dict))
        print(vertices)
        while len(vertices)>0:
            v = min(vertices)
            vertices.discard(v)
            print(f'popped {v}')
            component = ''
            if v in next_dict:
                while v in next_dict:
                    component += str(v)
                    vertices.discard(v)
                    v = next_dict[v]
                component += str(v)
                vertices.discard(v)
                if component[0] in prev_dict:
                    v = component[0]
                    while v in prev_dict:
                        v = prev_dict[v]
                        component = str(v) + component
                        vertices.discard(v)
            else:
                while v in prev_dict:
                    component = str(v) + component
                    vertices.discard(v)
                    v = prev_dict[v]
                component += str(v)
                vertices.discard(v)
                # v = prev_dict[v]
            components.append(component)
        print('done')
        return ', '.join(components)

    def init_config(self):
        def filter_edge_init(u,v):
            valid_edges = set(self.chromatin_edges).union(set(self.dsb_edges))
            return (u,v) in valid_edges or (v,u) in valid_edges
        init_view = nx.subgraph_view(self.graph, filter_edge=filter_edge_init)
        # return self.display_config(init_view)
        return init_view

    def final_config(self):
        # print('final')
        def filter_edge_final(u,v):
            valid_edges = set(self.chromatin_edges).union(set(self.rejoin_edges))
            return (u,v) in valid_edges or (v,u) in valid_edges
        final_view = nx.subgraph_view(self.graph, filter_edge=filter_edge_final)
        # print(final_view.edges)
        # return self.display_config(final_view)
        return final_view

    def draw(self):
        pos = nx.multipartite_layout(self.graph, subset_key='chromosome', align='horizontal')
        nx.draw_networkx_nodes(self.graph, pos, node_color = 'white', node_size=50)
        nx.draw_networkx_edges(self.graph, pos, width=1.5, edgelist=self.chromatin_edges, edge_color='black')
        nx.draw_networkx_edges(self.graph, pos, width=1.5, edgelist=self.dsb_edges, edge_color='red')
        nx.draw_networkx_edges(self.graph, pos, width=1.5, edgelist=self.rejoin_edges, arrows=True, 
                               connectionstyle=f'arc3, rad = 0.2', edge_color='green')
        nx.draw_networkx_labels(self.graph, pos, font_size=8)

    def chromatin_edge_twist(self, edge):
        pass

    def total_twist(self, chromosome):
        if chromosome >= self.chromosome:
            return self
        chrom_nodes = set(i for i in self.graph.nodes if self.graph.nodes[i]['chromosome'] == chromosome)
        start, stop = min(chrom_nodes), max(chrom_nodes)
        chromatin_edges = self.chromatin_edges_total_swap(start, stop)
        dsb_edges = self.dsb_edges_total_twist(start, stop)
        rejoin_edges = self.rejoin_edges_total_twist(start, stop)
        return AberrationMultigraph(chromatin_edges, dsb_edges, rejoin_edges)
        
    def chromatin_edges_total_twist(self, start, stop):
        return tuple(sorted([tuple(sorted((stop-u+start, stop-v+start))) 
                                if start <= u <= stop else (u,v) 
                                for (u,v) in self.chromatin_edges]))
    
    def dsb_edges_total_twist(self, start, stop):
        return tuple(sorted([tuple(sorted((stop-u+start, stop-v+start)))
                            if start <= u <=stop else (u,v)
                            for (u,v) in self.dsb_edges]))
    
    def rejoin_edges_total_twist(self, start, stop):
        rejoin_edges = []
        for (u,v) in self.rejoin_edges:
            if start <= u <= stop and start <= v <= stop:
                new_u, new_v = stop-u+start, stop-v+start
            elif start <= u <= stop:
                new_u, new_v = stop-u+start, v
            elif start <= v <= stop:
                new_u, new_v = u, stop-v+start
            else:
                new_u, new_v = u,v
            rejoin_edges.append(tuple(sorted((new_u,new_v))))
        return tuple(sorted(rejoin_edges))
        
    def total_swap(self, chrom_1, chrom_2):
        if chrom_1 >= self.chromosome or chrom_2 >= self.chromosome or chrom_1 == chrom_2:
            return self
        chrom_1, chrom_2 = (chrom_2, chrom_1) if chrom_2 < chrom_1 else (chrom_1, chrom_2)
        chrom_1_nodes = set(i for i in self.graph.nodes if self.graph.nodes[i]['chromosome'] == chrom_1)
        chrom_2_nodes = set(i for i in self.graph.nodes if self.graph.nodes[i]['chromosome'] == chrom_2)
        if len(chrom_1) != len(chrom_2):
            return self
        start_1, stop_1 = min(chrom_1_nodes), max(chrom_1_nodes)
        start_2, stop_2 = min(chrom_2_nodes), max(chrom_2_nodes)
        chromatin_edges = self.chromatin_edges_total_swap(start_1, start_2, stop_1, stop_2)
        dsb_edges = self.dsb_edges_total_swap(start_1, start_2, stop_1, stop_2)
        rejoin_edges = self.rejoin_edges_total_swap(start_1, start_2, stop_1, stop_2)
        return AberrationMultigraph(chromatin_edges, dsb_edges, rejoin_edges)

    def chromatin_edges_total_swap(self, start_1, start_2, stop_1, stop_2):
        chromatin_edges = []
        for (u,v) in self.chromatin_edges:
            if u < start_1 or v > stop_2:
                new_u, new_v = u,v
            elif start_1 <= u <= stop_1:
                new_u, new_v = u+stop_2-stop_1,v+stop_2-stop_1
            elif stop_1 < u < start_2:
                new_u, new_v = u-(stop_1-start_1)+(stop_2-start_2), v-(stop_1-start_1)+(stop_2-start_2)
            else:
                new_u, new_v = u-start_2+start_1, v-start_2+start_1
            chromatin_edges.append(tuple(sorted((new_u, new_v))))
        return tuple(sorted(chromatin_edges))
    
    def dsb_edges_total_swap(self, start_1, start_2, stop_1, stop_2):
        dsb_edges = []
        for (u,v) in self.dsb_edges:
            if u < start_1 or v > stop_2:
                new_u, new_v = u,v
            elif start_1 <= u <= stop_1:
                new_u, new_v = u+stop_2-stop_1,v+stop_2-stop_1
            elif stop_1 < u < start_2:
                new_u, new_v = u-(stop_1-start_1)+(stop_2-start_2), v-(stop_1-start_1)+(stop_2-start_2)
            else:
                new_u, new_v = u-start_2+start_1, v-start_2+start_1
            dsb_edges.append(tuple(sorted((new_u, new_v))))
        return tuple(sorted(dsb_edges))
    
    def rejoin_edges_total_swap(self, start_1, start_2, stop_1, stop_2):
        rejoin_edges = []
        for (u,v) in self.rejoin_edges:
            if u < start_1 or u > stop_2:
                new_u = u
            elif start_1 <= u <= stop_1:
                new_u = u+stop_2-stop_1
            elif stop_1 < u < start_2:
                new_u = u-(stop_1-start_1)+(stop_2-start_2)
            else:
                new_u = u-start_2+start_1
            if v < start_1 or v > stop_2:
                new_v = v
            elif start_1 <= v <= stop_1:
                new_v = v+stop_2-stop_1
            elif stop_1 < v < start_2:
                new_v = v-(stop_1-start_1)+(stop_2-start_2)
            else:
                new_v = v-start_2+start_1
            rejoin_edges.append(tuple(sorted((new_u, new_v))))
        return tuple(sorted(rejoin_edges))

    def save_to_file(self):
        #TODO: add method to save an AMG to a file.
        pass

    @staticmethod
    def load_from_file(file):
        #TODO: add method to load AMG from a file.
        pass
    
    def __hash__(self) -> int:
        return hash(self.dsb_edges+self.rejoin_edges)
    
    def __eq__(self, other):
        return True if self.chromatin_edges == other.chromatin_edges and \
                        self.dsb_edges == other.dsb_edges and \
                        self.rejoin_edges == other.rejoin_edges \
                    else False

if __name__ == '__main__':
    # amg_1 = AberrationMultigraph((('A','B'), ('C','D'), ('E','F'), ('G','H'), ('I','J'), ('K','L')), (('B','C'), ('D','E'), ('H','I'), ('J','K')), (('B', 'H'), ('D', 'J'), ('C', 'I'), ('E', 'K')))
    # amg_2 = AberrationMultigraph((('A','B'), ('C','D'), ('E','F'), ('G','H'), ('I','J'), ('K','L')), (('B','D'), ('C','E'), ('H','I'), ('J','K')), (('B', 'H'), ('D', 'J'), ('C', 'I'), ('E', 'K')))
    # amg_1 = AberrationMultigraph([(1,2), (3,4), (5,6), (7,8), (9,10), (11,12)], [(2,3),(4,5),(8,9),(10,11)], [(2,10), (3,11), (4,9), (5,8)])
    # amg_1 = AberrationMultigraph([(1,2),(3,4),(5,6),(7,8),(9,10)], [(2,3),(4,5),(6,7),(8,9)], [(2,9), (3,4),(5,7),(6,8)])
    amg_1 = AberrationMultigraph([(1,2),(3,4),(5,6),(7,8),(9,10),(11,12),(13,14),(15,16)], [(2,3),(6,7),(10,11),(14,15)], [(2,6),(3,11),(7,14),(10,15)])
    amg_2 = amg_1.total_swap(0,1)
    amg_3 = amg_1.total_swap(0,2)
    amg_4 = amg_2.total_swap(2,3)
    plt.subplot(2,2,1)
    amg_1.draw()
    plt.subplot(2,2,2)
    amg_2.draw()
    plt.subplot(2,2,3)
    amg_3.draw()
    plt.subplot(2,2,4)
    amg_4.draw()
    plt.show()
    # amg = AberrationMultigraph([1,2,3,4], [(2,3)], [(1,2), (3,4)])
    
    # print(amg.initial)
    # print(amg.final)
    # print(amg.diameter())