import heapq as hq
from aberration_multigraph.amg import AberrationMultigraph
from collections import defaultdict

class AMGGenerator:
    """
    This class represents an object that generates all possible AMGs with
     a specified number of chromosomes and an iterable specifying the number of
     double-strand breaks per chromosome.
    It uses a brute-force approach.
    """
    def __init__(self, num_chromosomes, num_dsbs, labels=None):
        """
        Parameters
        ----------
        num_chromosomes : int
            The number of chromosomes in the generated AMGs.
        num_dsbs : iterable
            The number of DSBs per chromosome.
        labels : str, optional
            Names of AMG vertices, by default None
        """
        if len(num_dsbs) != num_chromosomes:
            print('Mismatch')               #TODO raise ValueError
        self.num_chromosomes = num_chromosomes
        self.num_dsbs = num_dsbs
        self.vertices = self._get_labels() if labels is None else labels
        #TODO check length of labels matches number of vertices and all labels
        # are unique
        self.dsbs = self._get_dsbs()
        self.chromatins = self._get_chromatins()
        self.dsb_pair = dict()
        self.amg_counter = 1
        for i, j in self.dsbs:
            self.dsb_pair[i] = j
            self.dsb_pair[j] = i
        
    def generate_amgs(self):
        """This method constructs all possible AMGs for this DSB distribution.

        It returns a generator that builds the AMGs.
        
        Returns
        -------
        generator
            A generator that iterates over the AMGs with this DSB distribution.
        """
        dsb_vertices = set(u for u,_ in self.dsbs).union(set(v for _,v in self.dsbs))
        free_vertices = [(len(dsb_vertices)-2, i) for i in dsb_vertices]
        self.amg_counter = 0
        hq.heapify(free_vertices)
        return self._gen_rejoins([], free_vertices)

    def count_amgs(self):
        """Counts the number of AMGs with this DSB distribution.

        Returns
        -------
        int
            The number of AMGs wit this DSB distribution.
        """
        self.amg_counter = 0
        # TODO: count the AMGs without constructing them.
        return self.amg_counter

    def _gen_rejoins(self, rejoins, free_verts):
        """A helper method that recursively iterates over all possible rejoin 
         edge combinations and generates AMGs.

        Parameters
        ----------
        rejoins : list
            A list of pairs of rejoin edges created so far.
        free_verts : list
            Heap of free vertices ordered by number of free vertices at time of
             insertion into the heap.

        Yields
        ------
        AberrationMultigraph
            An AMG with this DSB distribution.
        """
        _, v = hq.heappop(free_verts)
        # If there is only one unmatched vertex left, pair it to this vertex 
        # and generate an AMG.
        # Else, pair this vertex with all remaining vertices and generate AMGs
        # recursively.
        if len(free_verts) == 1:
            _, u = hq.heappop(free_verts)
            amg = AberrationMultigraph(self.chromatins,
                                        self.dsbs,
                                        rejoins+[(v, u)],
                                        str(self.amg_counter))
            if amg.is_connected():
                self.amg_counter += 1
                yield amg
        else:
            for _, w in free_verts:
                if w != self.dsb_pair[v]:
                    new_free_verts = self._remaining_verts(free_verts, v, w)
                    new_rejoins = rejoins + [(v, w)]
                    yield from self._gen_rejoins(new_rejoins, new_free_verts)

    def _remaining_verts(self, free_verts, v, w):
        new_free_verts = []
        for priority, u in free_verts:
            if self.dsb_pair[u] == v or self.dsb_pair[u] == w:
                new_free_verts.append((priority-1, u))
            elif u != w:
                new_free_verts.append((priority-2, u))
        hq.heapify(new_free_verts)
        return new_free_verts
    
    def _get_dsbs(self):
        """Helper method to generate all DSB edges.

        Returns
        -------
        list
            List of DSB edges for this DSB distribution.
        """
        dsbs = []
        label = iter(self.vertices)
        for i in range(self.num_chromosomes):
            next(label)                      # first telomere in a chromosome
            for _ in range(self.num_dsbs[i]):
                dsbs.append((next(label), next(label)))
            next(label)                      # second telomere in the chromosome
        return dsbs

    def _get_labels(self):
        """Method to generate vertex names if these are not provided.

        Returns
        -------
        list
            List of vertex labels.
        """
        return list(range(2*sum(self.num_dsbs)+2*self.num_chromosomes))
    
    def _get_chromatins(self):
        """Method to generate chromatin edges.

        Returns
        -------
        list
            List of chromatin edges for this DSB distribution.
        """
        edges = []
        label = iter(self.vertices)
        for i in range(self.num_chromosomes):
            edges.append((next(label), next(label)))
            for j in range(self.num_dsbs[i]):
                edges.append((next(label), next(label)))
        return edges
    
    def summarize(self):
        """Prints summarized information of all AMGs with this DSB distribution.

        This method iterates over all AMGs with this DSB distribution and prints
         the frequency distribution of AMGs by cycle structure, diameter and
         girth.
        """
        cycles = defaultdict(int)
        diameters = defaultdict(int)
        # girths = defaultdict(int)         # Method not supported by networkx
        for amg in self.generate_amgs():
            cs = amg.cycle_structure()
            s_cs = [f'{i}*{cs[i]}' if cs[i]!=1 else str(i) for i in sorted(cs)]
            cycles['+'.join(s_cs)] += 1
            diameters[amg.diameter()] += 1
            # girths[amg.girth()] += 1
        print(f'TOTAL NUMBER OF AMGS: {self.amg_counter}')
        print('\nDISTRIBUTION BY CYCLE STRUCTURE')
        for c in cycles:
            print(f'{c}: {cycles[c]}')
        print('\nDISTRIBUTION BY DIAMETER')
        for d in sorted(diameters):
            print(f'{d}: {diameters[d]}')
        # print('\nDISTRIBUTION BY GIRTH')
        # for g in sorted(girths):
        #     print(f'{g}: {girths[g]}')

    def full_report(self, file):
        """Writes a detailed report containing the statistics of each AMG with
         this DSB distribution.

        Parameters
        ----------
        file : file
            The file where the report is written.
        """
        # TODO: make this a CSV file?
        with open(file, 'w') as output:
            output.write("""DIAMETER\t
                            GIRTH\t
                            CYCLE STRUCTURE\t\t
                            CYCLES\t\t\t
                            INIT CONFIG\t\t\t\t\t\t
                            FINAL CONFIG\n""")
            for amg in self.generate_amgs():
                diameter = amg.diameter()
                girth = amg.girth()
                cs = amg.cycle_structure()
                cycles = ','.join('('+ ','.join(str(c) for c in cycle)
                                    + ')' for cycle in amg.cycles())
                cycle_struct = '+'.join(str(cs[i]) for i in sorted(cs))
                init_config = ','.join(f'({u},{v})'
                                       for u,v in amg.init_config().edges)
                # init_config = amg.init_config()
                final_config = ','.join(f'({u},{v})'
                                        for u,v in amg.final_config().edges)
                # final_config = amg.final_config()
                output.write(f"""{diameter}\t\t
                                {girth}\t\t
                                {cycle_struct}\t\t\t
                                {cycles}\t\t
                                {init_config}\t\t
                                {final_config}\n""")

    def full_csv_report(self, file):
        #TODO: deprecated
        with open(file, 'w') as output:
            for amg in self.generate_amgs():
                diameter = amg.diameter()
                cs = amg.cycle_structure()
                cycle_struct = '+'.join(str(cs[i]) for i in sorted(cs))
                output.write(f'{diameter}, {cycle_struct}\n')