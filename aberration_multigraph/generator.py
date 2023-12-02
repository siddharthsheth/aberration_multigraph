import heapq as hq
from aberration_multigraph.amg import AberrationMultigraph
from collections import defaultdict
from matplotlib import pyplot as plt

class AMGGenerator:
    def __init__(self, num_chromosomes, num_dsbs, labels=None):
        self.num_chromosomes = num_chromosomes
        self.num_dsbs = num_dsbs
        self.vertices = self._generate_labels() if labels is None else labels
        self.dsbs = self._generate_dsbs()
        self.chromatins = self._generate_chromatin_edges()
        # self.free_ends = set(i for i in self.vertices if (i,i+1) in self.dsbs or (i-1, i) in self.dsbs)
        self.free_ends = set(u for (u,v) in self.dsbs).union(set(v for (u,v) in self.dsbs))
        self.dsb_pair = dict()
        self.amg_counter = 1
        # self.labels = labels
        for i, j in self.dsbs:
            self.dsb_pair[i] = j
            self.dsb_pair[j] = i
        
    def generate_amgs(self):
        # print(self.vertices)
        # print(self.dsbs)
        # print(self.free_ends)
        unpaired_vertices = [(len(self.free_ends)-2, i) for i in self.free_ends]
        hq.heapify(unpaired_vertices)
        rejoinings = []
        return self._generate_permutations(rejoinings, unpaired_vertices)

    def _generate_permutations(self, rejoinings, unpaired_vertices):
        _, vertex = hq.heappop(unpaired_vertices)
        if len(unpaired_vertices) == 1:
            _, last_vertex = hq.heappop(unpaired_vertices)
            # print(f'generating AMG: {rejoinings+[(vertex, last_vertex)]}')
            amg = AberrationMultigraph(self.chromatins, self.dsbs, rejoinings+[(vertex, last_vertex)], str(self.amg_counter))
            self.amg_counter += 1
            if amg.is_connected():
                yield amg
        else:
            for _, vertex_to_pair in unpaired_vertices:
                if vertex_to_pair != self.dsb_pair[vertex]:
                    new_unpaired_vertices = self._generate_new_unpaired_vertices(unpaired_vertices, vertex, vertex_to_pair)
                    # print(f'nested call: {new_unpaired_vertices}, {rejoinings+[(vertex, vertex_to_pair)]}')
                    yield from self._generate_permutations(rejoinings+[(vertex, vertex_to_pair)], new_unpaired_vertices)

    def _generate_new_unpaired_vertices(self, unpaired_vertices, vertex, vertex_to_pair):
        new_unpaired_vertices = []
        for priority, unpaired_vertex in unpaired_vertices:
            if self.dsb_pair[unpaired_vertex] == vertex or self.dsb_pair[unpaired_vertex] == vertex_to_pair:
                new_unpaired_vertices.append((priority-1, unpaired_vertex))
            elif unpaired_vertex != vertex_to_pair:
                new_unpaired_vertices.append((priority-2, unpaired_vertex))
        hq.heapify(new_unpaired_vertices)
        return new_unpaired_vertices
    
    def _generate_dsbs(self):
        dsbs = set()
        label = iter(self.vertices)
        for i in range(self.num_chromosomes):
            next(label)                                      # first telomere in a chromosome
            for _ in range(self.num_dsbs[i]):
                dsbs.add((next(label), next(label)))
                # label += 2                                  # the two free ends created by each DSB
            next(label)                                      # second telomere in the chromosome
        return dsbs

    def _generate_labels(self):
        return list(range(2*sum(self.num_dsbs)+2*self.num_chromosomes))
    
    def _generate_chromatin_edges(self):
        edges = []
        label = iter(self.vertices)
        for i in range(self.num_chromosomes):
            edges.append((next(label), next(label)))
            for j in range(self.num_dsbs[i]):
                edges.append((next(label), next(label)))
        return edges
    
    def summarize(self):
        cycles = defaultdict(int)
        diameters = defaultdict(int)
        for amg in self.generate_amgs():
            c = [len(cycle) for cycle in amg.cycles()]
            diameters[amg.diameter()] += 1
            cycles['+'.join(str(i) for i in sorted(c))] += 1
        print('DISTRIBUTION BY CYCLE STRUCTURE')
        for c in cycles:
            print(f'{c}: {cycles[c]}')
        print('\nDISTRIBUTION BY DIAMETER')
        for d in sorted(diameters):
            print(f'{d}: {diameters[d]}')

    def full_report(self, file):
        with open(file, 'w') as output:
            output.write('DIAMETER\tCYCLE STRUCTURE\t\tCYCLES\t\t\tINIT CONFIG\t\t\t\t\t\tFINAL CONFIG\n')
            for amg in self.generate_amgs():
                diameter = amg.diameter()
                cycles = ','.join('('+','.join(str(c) for c in cycle)+')' for cycle in amg.cycles())
                cycle_struct = '+'.join(str(i) for i in sorted([len(cycle) for cycle in amg.cycles()]))
                init_config = ','.join(f'({u},{v})' for u,v in amg.init_config().edges)
                # init_config = amg.init_config()
                final_config = ','.join(f'({u},{v})' for u,v in amg.final_config().edges)
                # final_config = amg.final_config()
                output.write(f'{diameter}\t\t{cycle_struct}\t\t\t{cycles}\t\t{init_config}\t\t{final_config}\n')

    def full_csv_report(self, file):
        with open(file, 'w') as output:
            for amg in self.generate_amgs():
                diameter = amg.diameter()
                cycle_struct = '+'.join(str(i) for i in sorted([len(cycle) for cycle in amg.cycles()]))
                output.write(f'{diameter},{cycle_struct}\n')



if __name__ == '__main__':
    # amg_gen = AMGGenerator(1, (2,), ('ABCDEFGHIJKL'))
    amg_gen = AMGGenerator(2, (2,2), ('ABCDEFGHIJKLMN'))
    i = 1
    for graph in amg_gen.generate_amgs():
        # d = graph.diameter()
        # c = len(graph.cycles())
        # print(graph.rejoin_edges)
        if graph.rejoin_edges == (('B', 'H'), ('D', 'I'), ('J', 'C'), ('E', 'K')) or graph.rejoin_edges == (('B', 'H'), ('D', 'J'), ('C', 'I'), ('E', 'K')):
            plt.subplot(2,2,i)
            # ax = plt.gca()
            # ax.set_xlim([xmin, xmax])
            # ax.set_ylim([0, 2])
            graph.draw()
            i += 1
    plt.show()
    # amg_gen.full_report('dsb_2_2_1.txt')


    # configs = [(1,(2,)),(1,(3,)),(2,(1,1)),(2,(2,1)),(2,(2,2)),(2,(3,2)),(3,(1,1,1)),(3,(2,1,1)),(3,(3,1,1)),(3,(2,2,1)),(3,(2,2,2))]
    # configs = [(4,(1,1,1,1)), ((4,(2,2,2,2)))]
    # configs = [(3, (4,1,1)), (2, (5,1)), (2, (4,2)), (2, (3,3)), (1, (6,)), (4, (2,2,1,1))]
    # configs = [(3,(3,3,2))]
    # for chromosome, dsbs in configs:
    #     print(str(chromosome)+str(dsbs))
    #     amg_gen = AMGGenerator(chromosome, dsbs)
    #     str_dsb = '_'.join(str(i) for i in dsbs)
    #     amg_gen.summarize()
    #     amg_gen.full_csv_report(f'dsb_{chromosome}_{str_dsb}.csv')