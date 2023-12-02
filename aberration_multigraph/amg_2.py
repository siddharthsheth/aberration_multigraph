import networkx as nx
# from matplotlib import pyplot as plt

class AberrationMultigraph:
    def __init__(self, vertices, dsbs, rejoinings):
        self.initial = nx.Graph()
        self.final = nx.Graph()
        self.exchange = nx.Graph()
        for i in vertices[:-1]:
            if (i, i+1) not in dsbs:
                self.initial.add_edge(i, i+1, color='chromatin')
                self.final.add_edge(i, i+1, color='chromatin')

        for i, j in dsbs:
            self.initial.add_edge(i,j, color='dsb')
            self.exchange.add_edge(i,j, color='dsb')

        for i, j in rejoinings:
            self.exchange.add_edge(i,j, color='misrejoining')
            self.final.add_edge(i,j, color='misrejoining')

    def diameter(self):
        master_graph = nx.Graph(self.initial)
        master_graph.add_edges_from(self.final.edges)
        return nx.distance_measures.diameter(master_graph)
    
    def cycles(self):
        return nx.cycles.cycle_basis(self.exchange)
    
    def draw(self, colormap={'chromatin':'black', 'dsb':'red', 'misrejoining':'green'}):
        master_graph = nx.Graph(self.initial)
        master_graph.add_edges_from(self.final.edges)
        chromatin_edges = [(u,v) for u,v,color in master_graph.edges.data('color') if color=='chromatin']
        dsb_edges = [(u,v) for u,v,color in master_graph.edges.data('color') if color=='dsb']
        rejoin_edges = [(u,v) for u,v,color in master_graph.edges.data('color') if color=='misrejoining']
        pos = nx.spring_layout(master_graph, seed=3113794652)
        nx.draw_networkx_nodes(master_graph, pos)
        nx.draw_networkx_edges(master_graph, pos, edgelist=chromatin_edges, edge_color='black')
        nx.draw_networkx_edges(master_graph, pos, edgelist=dsb_edges, edge_color='red')
        nx.draw_networkx_edges(master_graph, pos, edgelist=rejoin_edges, edge_color='blue')
        nx.draw_networkx_labels(master_graph, pos)
        # draw(master_graph, with_labels=True, edge_color=colormap)

if __name__ == '__main__':
    amg = AberrationMultigraph([1,2,3,4], [(2,3)], [(1,2), (3,4)])

    
    # print(amg.initial)
    # print(amg.final)
    # print(amg.diameter())