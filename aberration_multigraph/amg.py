import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
from collections import Counter
import pickle

class AberrationMultigraph:
    """
    A class to store aberration multigraphs as defined by Sachs et al.
    An AMG corresponding to a chromosome aberration is defined by
     chromatin edges, double-strand break edges and rejoin edges.
    The graph is stored as a networkx.Graph object.
    This class includes methods to compute basic properties of AMGs such as 
     diameter, girth and cycle structure.
    There is also a method to draw an AMG, but it should only be used for simple
     AMGs.
    Finally, there are also methods to perform combinatorial transformations on
     AMGs such as chromosome edge swaps, total swaps and total twists.
    Additionally, AMGs can be saved to files and recovered from previously saved
     files.

    Attributes
    ----------
    chromatin_edges : tuple
        Collection of chromatin edges in the AMG.
    dsb_edges : tuple
        Collection of edges corresponding to double-strand breaks.
    rejoin_edges : tuple
        Collection of edges corresponding to rejoins.
    graph : networkx.Graph
        The aberration multigraph defined by the above edges.
    num_chromosome : int
        The number of chromosomes in the AMG.
    """
    def __init__(self, chromatin_edges, dsb_edges, rejoin_edges, name=''):
        """
        Parameters
        ----------
        chromatin_edges : iterable
            Collection of chromatin edges in the AMG.
        dsb_edges : iterable
            Collection of edges corresponding to double-strand breaks.
        rejoin_edges : iterable
            Collection of edges corresponding to rejoins.
        name : str, optional
            A name for the aberration multigraph, by default ''
        """
        self.graph = nx.Graph()
        self.chromatin_edges = tuple(sorted(
                            tuple(sorted(edge)) for edge in chromatin_edges))
        self.dsb_edges = tuple(sorted(
                            tuple(sorted(edge)) for edge in dsb_edges))
        self.rejoin_edges = tuple(sorted(
                            tuple(sorted(edge)) for edge in rejoin_edges))
        self.name = name
        self.graph.add_edges_from(self.chromatin_edges, color='chromatin')
        self.graph.add_edges_from(self.dsb_edges, color='dsb')
        self.graph.add_edges_from(self.rejoin_edges, color='misrejoining')
        chromosome = 0
        dsb_ends = 0
        dsb_vertices = set(u for (u,_) in self.dsb_edges).union(
                                            set(v for (_,v) in self.dsb_edges)
                                            )
        for i in self.graph.nodes:
            if i not in dsb_vertices:
                dsb_ends += 1
            self.graph.nodes[i]['chromosome'] = chromosome
            if dsb_ends == 2:
                chromosome += 1
                dsb_ends = 0
        self.num_chromosome = chromosome

    def diameter(self):
        """Method to compute diameter of the AMG.

        Returns
        -------
        int
            Diameter of the AMG.
            Returns numpy.inf if the graph is not connected.
        """
        return (np.inf if not self.is_connected()
                        else nx.distance_measures.diameter(self.graph))
    
    def girth(self):
        """Method to compute girth of the AMG.

        Returns
        -------
        int
            Girth of the AMG.
            Returns numpy.inf if there are no cycles in case of an incomplete AMG.
        """
        return min([len(c) for c in nx.simple_cycles(self.graph)])
        # As per networkx documentation the following should work, but doesn't.
        # return nx.cycles.girth(self.graph)

    def cycles(self):
        """Method to compute the cycles in the AMG.

        Returns
        -------
        iterable
            Collection of cycles present in the AMG.
        """
        return nx.cycles.cycle_basis(nx.Graph([
            (u,v) for u,v,color in self.graph.edges.data('color')
                    if color!='chromatin' ]))
    
    def cycle_structure(self):
        """Method to compute cycle structure of the AMG.

        Returns
        -------
        Counter
            Counts the number of cycles by length in the AMG.
        """
        return Counter(sorted([len(cycle) for cycle in self.cycles()]))
    
    def is_connected(self):
        """Method to check whether the AMG is connected or not.

        Returns
        -------
        bool
            True if the AMG is connected, False otherwise.
        """
        return nx.is_connected(self.graph)
    
    #TODO: Display cycle structure in a pretty way.

    # def display_config(self, view):
    #     next_dict = {u: v for u,v in view.edges}
    #     prev_dict = {v: u for u,v in view.edges}
    #     print(next_dict)
    #     print(prev_dict)
    #     components = []
    #     vertices = set(v for v in self.graph.nodes if (v in next_dict)^(v in prev_dict))
    #     print(vertices)
    #     while len(vertices)>0:
    #         v = min(vertices)
    #         vertices.discard(v)
    #         print(f'popped {v}')
    #         component = ''
    #         if v in next_dict:
    #             while v in next_dict:
    #                 component += str(v)
    #                 vertices.discard(v)
    #                 v = next_dict[v]
    #             component += str(v)
    #             vertices.discard(v)
    #             if component[0] in prev_dict:
    #                 v = component[0]
    #                 while v in prev_dict:
    #                     v = prev_dict[v]
    #                     component = str(v) + component
    #                     vertices.discard(v)
    #         else:
    #             while v in prev_dict:
    #                 component = str(v) + component
    #                 vertices.discard(v)
    #                 v = prev_dict[v]
    #             component += str(v)
    #             vertices.discard(v)
    #             # v = prev_dict[v]
    #         components.append(component)
    #     print('done')
    #     return ', '.join(components)

    # def init_config(self):
    #     def filter_edge_init(u,v):
    #         valid_edges = set(self.chromatin_edges).union(set(self.dsb_edges))
    #         return (u,v) in valid_edges or (v,u) in valid_edges
    #     init_view = nx.subgraph_view(self.graph, filter_edge=filter_edge_init)
    #     # return self.display_config(init_view)
    #     return init_view

    # def final_config(self):
    #     # print('final')
    #     def filter_edge_final(u,v):
    #         valid_edges = set(self.chromatin_edges).union(set(self.rejoin_edges))
    #         return (u,v) in valid_edges or (v,u) in valid_edges
    #     final_view = nx.subgraph_view(self.graph, filter_edge=filter_edge_final)
    #     # print(final_view.edges)
    #     # return self.display_config(final_view)
    #     return final_view

    def draw(self):
        """Method to draw an AMG.

        This method should be called after instantiating a
         matplotlib.pyplot.figure object.
        """
        pos = nx.multipartite_layout(self.graph,
                                        subset_key='chromosome',
                                        align='horizontal')
        nx.draw_networkx_nodes(self.graph,
                                pos,
                                node_color = 'white',
                                node_size=50)
        nx.draw_networkx_edges(self.graph,
                                pos,
                                width=1.5,
                                edgelist=self.chromatin_edges,
                                edge_color='black')
        nx.draw_networkx_edges(self.graph, 
                                pos,
                                width=1.5,
                                edgelist=self.dsb_edges,
                                edge_color='red')
        nx.draw_networkx_edges(self.graph,
                                pos,
                                width=1.5,
                                edgelist=self.rejoin_edges,
                                arrows=True, 
                                connectionstyle=f'arc3, rad = 0.2',
                                edge_color='green')
        nx.draw_networkx_labels(self.graph, pos, font_size=8)

    def chromatin_edge_twist(self, edge):
        #TODO
        pass

    def total_twist(self, chromosome):
        """Method to generate an AMG by reversing the orientation of a
         chromosome in this AMG.

        Parameters
        ----------
        chromosome : int
            A number describing the chromosome that is to be twisted.

        Returns
        -------
        AberrationMultigraph
            An AMG where the input chromosome is reversed in orientation.
        """
        if chromosome >= self.num_chromosome:
            return self
        chrom_nodes = set(i for i in self.graph.nodes \
                          if self.graph.nodes[i]['chromosome'] == chromosome)
        start, stop = min(chrom_nodes), max(chrom_nodes)
        chromatin_edges = self._chromatin_edges_total_twist(start, stop)
        dsb_edges = self._dsb_edges_total_twist(start, stop)
        rejoin_edges = self._rejoin_edges_total_twist(start, stop)
        return AberrationMultigraph(chromatin_edges, dsb_edges, rejoin_edges)
        
    def _chromatin_edges_total_twist(self, start, stop):
        """Helper method for total_twist to manage chromatin edges.

        Parameters
        ----------
        start : int
            The first node of the chromosome to be twisted.
        stop : int
            The last node of the chromosome to be twisted.

        Returns
        -------
        tuple
            A tuple containing the chromatin edges of the new AMG with
             chromatin edges from the input chromosome reversed.
        """
        return tuple(sorted([tuple(sorted((stop-u+start, stop-v+start))) 
                                if start <= u <= stop else (u,v) 
                                for (u,v) in self.chromatin_edges]))
    
    def _dsb_edges_total_twist(self, start, stop):
        """Helper method for total_twist to manage DSB edges.

        Parameters
        ----------
        start : int
            The first node of the chromosome to be twisted.
        stop : int
            The last node of the chromosome to be twisted.

        Returns
        -------
        tuple
            A tuple containing the DSB edges of the new AMG with
             DSB edges from the input chromosome reversed.
        """
        return tuple(sorted([tuple(sorted((stop-u+start, stop-v+start)))
                                if start <= u <=stop else (u,v)
                                for (u,v) in self.dsb_edges]))
    
    def _rejoin_edges_total_twist(self, start, stop):
        """Helper method for total_twist to manage rejoin edges.

        Parameters
        ----------
        start : int
            The first node of the chromosome to be twisted.
        stop : int
            The last node of the chromosome to be twisted.

        Returns
        -------
        tuple
            A tuple containing updated rejoin edges if the input chromosome is
             reversed.
        """
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
        """Method to swap locations of two chromosomes in an AMG.

        Parameters
        ----------
        chrom_1 : int
            An integer representing the first chromosome to be swapped.
        chrom_2 : int
            An integer representing the second chromosome to be swapped.

        Returns
        -------
        AberrationMultigraph
            The method returns a new AMG where the locations of chrom1 and
             chrom2 are swapped.
        """
        if (chrom_1 >= self.num_chromosome 
                or chrom_2 >= self.num_chromosome 
                or chrom_1 == chrom_2):
            return self
        chrom_1, chrom_2 = ((chrom_2, chrom_1) if chrom_2 < chrom_1
                                                else (chrom_1, chrom_2))
        chrom_1_nodes = set(i for i in self.graph.nodes
                                if self.graph.nodes[i]['chromosome'] == chrom_1)
        chrom_2_nodes = set(i for i in self.graph.nodes
                                if self.graph.nodes[i]['chromosome'] == chrom_2)
        if len(chrom_1_nodes) != len(chrom_2_nodes):
            return self
        start_1, stop_1 = min(chrom_1_nodes), max(chrom_1_nodes)
        start_2, stop_2 = min(chrom_2_nodes), max(chrom_2_nodes)
        chromatin_edges = self._chromatin_edges_total_swap(start_1,
                                                            start_2,
                                                            stop_1,
                                                            stop_2)
        dsb_edges = self._dsb_edges_total_swap(start_1,
                                                start_2,
                                                stop_1,
                                                stop_2)
        rejoin_edges = self._rejoin_edges_total_swap(start_1,
                                                    start_2,
                                                    stop_1,
                                                    stop_2)
        return AberrationMultigraph(chromatin_edges, dsb_edges, rejoin_edges)

    def _chromatin_edges_total_swap(self, start_1, start_2, stop_1, stop_2):
        """Helper method for total_swap to manage chromatin edges.

        Parameters
        ----------
        start_1 : int
            The first node of the first chromosome to be swapped.
        start_2 : int
            The first node of the second chromosome to be swapped.
        stop_1 : int
            The last node of the first chromosome to be swapped.
        stop_2 : int
            The last node of the second chromosome to be swapped.

        Returns
        -------
        tuple
            A tuple containing the chromatin edges of the new AMG where
             chromatin edges from chrom1 and chrom2 are swapped.
        """
        chromatin_edges = []
        for (u,v) in self.chromatin_edges:
            if u < start_1 or v > stop_2:
                new_u, new_v = u,v
            elif start_1 <= u <= stop_1:
                new_u, new_v = u+stop_2-stop_1,v+stop_2-stop_1
            elif stop_1 < u < start_2:
                new_u = u-(stop_1-start_1)+(stop_2-start_2)
                new_v = v-(stop_1-start_1)+(stop_2-start_2)
            else:
                new_u, new_v = u-start_2+start_1, v-start_2+start_1
            chromatin_edges.append(tuple(sorted((new_u, new_v))))
        return tuple(sorted(chromatin_edges))
    
    def _dsb_edges_total_swap(self, start_1, start_2, stop_1, stop_2):
        """Helper method for total_swap to manage DSB edges.

        Parameters
        ----------
        start_1 : int
            The first node of the first chromosome to be swapped.
        start_2 : int
            The first node of the second chromosome to be swapped.
        stop_1 : int
            The last node of the first chromosome to be swapped.
        stop_2 : int
            The last node of the second chromosome to be swapped.

        Returns
        -------
        tuple
            A tuple containing the DSB edges of the new AMG where
             DSB edges from chrom1 and chrom2 are swapped.
        """
        dsb_edges = []
        for (u,v) in self.dsb_edges:
            if u < start_1 or v > stop_2:
                new_u, new_v = u,v
            elif start_1 <= u <= stop_1:
                new_u, new_v = u+stop_2-stop_1,v+stop_2-stop_1
            elif stop_1 < u < start_2:
                new_u = u-(stop_1-start_1)+(stop_2-start_2)
                new_v = v-(stop_1-start_1)+(stop_2-start_2)
            else:
                new_u, new_v = u-start_2+start_1, v-start_2+start_1
            dsb_edges.append(tuple(sorted((new_u, new_v))))
        return tuple(sorted(dsb_edges))
    
    def _rejoin_edges_total_swap(self, start_1, start_2, stop_1, stop_2):
        """Helper method for total_swap to manage rejoin edges.

        Parameters
        ----------
        start_1 : int
            The first node of the first chromosome to be swapped.
        start_2 : int
            The first node of the second chromosome to be swapped.
        stop_1 : int
            The last node of the first chromosome to be swapped.
        stop_2 : int
            The last node of the second chromosome to be swapped.

        Returns
        -------
        tuple
            A tuple containing the updated rejoin edges resulting from the swap
             of chrom1 and chrom2.
        """
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

    def save_to_file(self, filename='', path=''):
        if filename == '':
            filename = str(self.name)
        with open(path+filename+'.pickle', 'wb') as file:
            pickle.dump(self, file)

    @staticmethod
    def load_from_file(filename, path=''):
        return pickle.load(open(path+filename, 'rb'))
    
    def __hash__(self) -> int:
        return hash(self.dsb_edges+self.rejoin_edges)
    
    def __eq__(self, other):
        return (True if self.chromatin_edges == other.chromatin_edges
                        and self.dsb_edges == other.dsb_edges 
                        and self.rejoin_edges == other.rejoin_edges
                    else False)