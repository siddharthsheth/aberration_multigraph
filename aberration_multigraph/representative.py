from aberration_multigraph.generator import AMGGenerator
from matplotlib import pyplot as plt
import matplotlib as mpl
import networkx as nx

class AMGRepresentative():
    def __init__(self, num_chromosomes, num_dsbs):
        self.chromosomes = num_chromosomes
        self.dsbs = num_dsbs
        self.amg_generator = AMGGenerator(num_chromosomes, num_dsbs, None)
        self.num_operations = num_chromosomes*(num_chromosomes+1)/2
        # self.amgs = self.amg_generator.generate_amgs()

    def compute_amg_adjacency_list(self):
        self.adjacency_list = dict()
        self.rejoin_name_map = dict()
        for amg in self.amg_generator.generate_amgs():
            # self.rejoin_name_map[(amg.dsb_edges,amg.rejoin_edges)] = amg.name
            self.rejoin_name_map[amg.rejoin_edges] = amg.name
            for k in range(self.chromosomes):
                chrom_nodes = set(i for i in amg.graph.nodes if amg.graph.nodes[i]['chromosome'] == k)
                start, stop = min(chrom_nodes), max(chrom_nodes)
                if amg in self.adjacency_list:
                    # self.adjacency_list[amg].add((amg.dsb_edges_total_twist(start,stop), amg.rejoin_edges_total_twist(start, stop), k+1))
                    self.adjacency_list[amg].add((amg._rejoin_edges_total_twist(start, stop), k+1))
                else:
                    # self.adjacency_list[amg] = {(amg.dsb_edges_total_twist(start,stop), amg.rejoin_edges_total_twist(start, stop), k+1)}
                    self.adjacency_list[amg] = {(amg._rejoin_edges_total_twist(start, stop), k+1)}

            operation_number = self.chromosomes+1
            for i in range(self.chromosomes):
                for j in range(i, self.chromosomes):
                    # if chrom_1 >= amg.chromosome or chrom_2 >= amg.chromosome or chrom_1 == chrom_2:
                        # return amg.dsb_edges+amg.rejoin_edges
                    chrom_1, chrom_2 = (j, i) if j < i else (i, j)
                    chrom_1_nodes = set(k for k in amg.graph.nodes if amg.graph.nodes[k]['chromosome'] == chrom_1)
                    chrom_2_nodes = set(k for k in amg.graph.nodes if amg.graph.nodes[k]['chromosome'] == chrom_2)
                    if len(chrom_1_nodes) == len(chrom_2_nodes):
                        start_1, stop_1 = min(chrom_1_nodes), max(chrom_1_nodes)
                        start_2, stop_2 = min(chrom_2_nodes), max(chrom_2_nodes)
                        # self.adjacency_list[amg].add((amg.dsb_edges_total_swap(start_1,start_2,stop_1,stop_2), amg.rejoin_edges_total_swap(start_1,start_2,stop_1,stop_2), operation_number))
                        self.adjacency_list[amg].add((amg._rejoin_edges_total_swap(start_1,start_2,stop_1,stop_2), operation_number))
                        operation_number += 1
    
    def compute_representative_graph(self):
        self.rep_graph = nx.Graph()
        for amg in self.adjacency_list:
            # for nbr_dsb_edge, nbr_rejoin_edge, operation in self.adjacency_list[amg]:
            for nbr_rejoin_edge, operation in self.adjacency_list[amg]:
                # nbr = self.rejoin_name_map[(nbr_dsb_edge, nbr_rejoin_edge)]
                nbr = self.rejoin_name_map[nbr_rejoin_edge]
                if nbr != amg.name:
                    self.rep_graph.add_edge(amg.name, nbr, color=operation)

    def draw_representative_graph(self):
        pos = nx.spring_layout(self.rep_graph)
        cmap = mpl.colormaps['tab20']
        # color_palette = [i for i in range(self.num_operations)]
        colors = [cmap(self.rep_graph[u][v]['color']) for u,v in self.rep_graph.edges]
        nx.draw_networkx_nodes(self.rep_graph, pos, node_color = 'white', node_size=50)
        nx.draw_networkx_labels(self.rep_graph, pos, font_size=8)
        nx.draw_networkx_edges(self.rep_graph, pos, width=1, edgelist=self.rep_graph.edges, edge_color=colors)

    # def compute_representatives(self):
    #     self.groups = dict()
    #     computed_amgs = set()
    #     for amg in self.amg_generator.generate_amgs():
    #         print(f'Computed AMGs = {computed_amgs}')
    #         print(f'AMG has rejoin edges = {amg.rejoin_edges}')
    #         if amg.rejoin_edges in computed_amgs:
    #             print('repeat AMG found')
    #             continue
    #         generated_amgs = {amg.twist(i) for i in range(self.chromosomes)}
    #         print(f'There are {len(generated_amgs)} after twists')
    #         # print(f'Twisted AMGs for {amg} are {generated_amgs}')
    #         for i in range(self.chromosomes):
    #             for j in range(i, self.chromosomes):
    #                 generated_amgs.add(amg.swap(i,j))
    #             # generated_amgs = generated_amgs.union({amg.swap(i,j) for j in range(i, self.chromosomes)})
    #         print(f'There are {len(generated_amgs)} in total')
    #         generated_amgs.add(amg)
    #         # print(f'Generated AMGs for {amg} are {generated_amgs}')
    #         self.groups[amg] = generated_amgs
    #         # other_set = {(amg.swap(i,j) for j in range(i, self.chromosomes)) for i in range(self.chromosomes)}
    #         computed_amgs = computed_amgs.union(set(tuple(sorted(g.rejoin_edges)) for g in generated_amgs))
    #         # computed_amgs.add(amg.rejoin_edges)


if __name__ == '__main__':
    amg_rep = AMGRepresentative(2, (2,2,1))
    amg_rep.compute_amg_adjacency_list()
    # print(amg_rep.adjacency_list)
    # print(amg_rep.rejoin_name_map)
    amg_rep.compute_representative_graph()
    amg_rep.draw_representative_graph()
    # nx.draw(amg_rep.rep_graph)
    # print(len(amg_rep.groups))
    # for graph in amg_rep.groups.keys():
    #     plt.figure()
    #     i = 1
    #     for amg in amg_rep.groups[graph]:
    #         plt.subplot(4,2,i)
    #         amg.draw()
    #         i += 1
    plt.show()
