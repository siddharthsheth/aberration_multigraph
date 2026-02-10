"""
AMG generation utilities.

This module provides a brute-force generator for aberration multigraphs (AMGs)
with a fixed number of chromosomes and a fixed distribution of double-strand
breaks (DSBs) per chromosome.

Model
-----
Each chromosome contributes two telomeres and a specified number of DSBs.
Each DSB introduces two free ends, which must be paired via *rejoin edges*.
An AMG corresponds to a perfect matching of all free ends such that no end is
rejoined with its original DSB partner. Only connected multigraphs are retained.

Due to the combinatorial explosion in the number of perfect matchings, this
implementation is intended for small instances and exploratory use.
"""

import heapq as hq
from aberration_multigraph.amg import AberrationMultigraph
from collections import defaultdict

class AMGGenerator:
    """
    Generator for all aberration multigraphs (AMGs) with a fixed DSB distribution.

    The generator enumerates all valid rejoinings of DSB ends using a recursive
    backtracking algorithm with simple constraint propagation. At each step, the
    most constrained free vertex is paired first to reduce branching.
    """
    def __init__(self, num_chromosomes, num_dsbs, labels=None):
        """
        Initialize an AMG generator.

        Parameters
        ----------
        num_chromosomes : int
            The number of chromosomes in the generated AMGs.
        num_dsbs : iterable
            The number of DSBs per chromosome.
        labels : str, optional
            Vertex labels. If ``None``, vertices are labeled consecutively
            by integers.

        Notes
        -----
        The total number of vertices is

        ``2 * sum(num_dsbs) + 2 * num_chromosomes``.
        """

        if len(num_dsbs) != num_chromosomes:
            print('Mismatch')               #TODO raise ValueError
        self.num_chromosomes = num_chromosomes
        self.num_dsbs = num_dsbs
        self.vertices = self._get_labels() if labels is None else labels
        #TODO check length of labels matches number of vertices and all labels
        # are unique

        # Structural edges
        self.dsbs = self._get_dsbs()
        self.chromatins = self._get_chromatins()

        # Map each DSB end to its original partner
        self.dsb_pair = dict()
        for i, j in self.dsbs:
            self.dsb_pair[i] = j
            self.dsb_pair[j] = i
        
        # Counter for generated AMGs
        self.amg_counter = 1
        
    def generate_amgs(self):
        """
        Generate all connected AMGs consistent with this DSB distribution.

        Returns
        -------
        generator of AberrationMultigraph
            A generator yielding AMGs one at a time.

        Side Effects
        ------------
        Resets and updates ``self.amg_counter``.
        """
        dsb_vertices = set(u for u,_ in self.dsbs).union(set(v for _,v in self.dsbs))
        free_vertices = [(len(dsb_vertices)-2, i) for i in dsb_vertices]
        self.amg_counter = 0
        hq.heapify(free_vertices)
        return self._gen_rejoins([], free_vertices)

    def count_amgs(self):
        """
        Count the number of AMGs with this DSB distribution.

        Returns
        -------
        int
            Number of AMGs.

        Notes
        -----
        This method is currently a placeholder. Counting is performed implicitly
        during generation.
        """
        self.amg_counter = 0
        # TODO: count the AMGs without constructing them.
        return self.amg_counter

    def _gen_rejoins(self, rejoins, free_verts):
        """
        Recursively enumerate all valid rejoin matchings.

        Parameters
        ----------
        rejoins : list of tuple
            Rejoin edges fixed so far.
        free_verts : list
            Heap of unmatched vertices, prioritized by remaining matching
            flexibility.

        Yields
        ------
        AberrationMultigraph
            A connected AMG obtained by completing the current partial matching.
        """
        if len(free_verts) == 0:
            return
        _, v = hq.heappop(free_verts)
        # Base case: if there is only one unmatched vertex left,
        #  pair it to this vertex and generate an AMG.
        if len(free_verts) == 1:
            _, u = hq.heappop(free_verts)
            amg = AberrationMultigraph(self.chromatins,
                                        self.dsbs,
                                        rejoins+[(v, u)],
                                        str(self.amg_counter))
            if amg.is_connected():
                self.amg_counter += 1
                yield amg
        # Else, pair this vertex with all remaining vertices and generate AMGs
        #  recursively.
        else:
            for _, w in free_verts:
                if w != self.dsb_pair[v]:
                    new_free_verts = self._remaining_verts(free_verts, v, w)
                    new_rejoins = rejoins + [(v, w)]
                    yield from self._gen_rejoins(new_rejoins, new_free_verts)

    def _remaining_verts(self, free_verts, v, w):
        """
        Update the heap of free vertices after pairing ``v`` and ``w``.

        Parameters
        ----------
        free_verts : list
            Current heap of free vertices.
        v, w : int
            Vertices paired at the current recursion step.

        Returns
        -------
        list
            Updated heap reflecting reduced matching options.
        """
        new_free_verts = []
        for priority, u in free_verts:
            if self.dsb_pair[u] == v or self.dsb_pair[u] == w:
                new_free_verts.append((priority-1, u))
            elif u != w:
                new_free_verts.append((priority-2, u))
        hq.heapify(new_free_verts)
        return new_free_verts
    
    def _get_dsbs(self):
        """
        Generate all DSB edges implied by the chromosome specification.

        Returns
        -------
        list of tuple
            DSB edges.
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
        """
        Generate default vertex labels.

        Returns
        -------
        list
            Consecutive integer labels.
        """
        return list(range(2*sum(self.num_dsbs)+2*self.num_chromosomes))
    
    def _get_chromatins(self):
        """
        Generate chromatin edges encoding chromosome structure.

        Returns
        -------
        list of tuple
            Chromatin edges.
        """
        edges = []
        label = iter(self.vertices)
        for i in range(self.num_chromosomes):
            edges.append((next(label), next(label)))
            for j in range(self.num_dsbs[i]):
                edges.append((next(label), next(label)))
        return edges
    
    def summarize(self):
        """
        Print summary statistics over all generated AMGs.

        Reported statistics include:
        - Total number of AMGs
        - Distribution by cycle structure
        - Distribution by diameter
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
        """
        Write a detailed per-AMG report to a file.

        Parameters
        ----------
        file : str
            Path to the output file.
        """
        # TODO: make this a CSV file?
        with open(file, 'w') as output:
            output.write("""DIAMETER\t
                            CYCLE STRUCTURE\t\t
                            CYCLES\t\t\t
                            INIT CONFIG\t\t\t\t\t\t
                            FINAL CONFIG\n""")
            for amg in self.generate_amgs():
                diameter = amg.diameter()
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
                                {cycle_struct}\t\t\t
                                {cycles}\t\t
                                {init_config}\t\t
                                {final_config}\n""")