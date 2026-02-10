"""
Incomplete aberration multigraphs (IAMGs).

This module defines the concept of an **incomplete aberration multigraph**,
which represents a partially specified aberration multigraph (AMG) in which
only a subset of rejoin edges has been fixed.

Incomplete AMGs arise naturally during recursive generation and enumeration
procedures, where rejoinings are constructed incrementally before a full
perfect matching is obtained.

An incomplete AMG consists of:

- A fixed chromatin structure
- A fixed set of DSB edges
- A partial set of rejoin edges

and supports basic graph-theoretic queries and extension operations.

.. note::
   Incomplete AMGs are intended as intermediate combinatorial objects and are
   not required to satisfy connectivity or completeness constraints.
"""

from aberration_multigraph.amg import AberrationMultigraph
import heapq as hq

class IncompleteAMG(AberrationMultigraph):
    """
    Representation of an incomplete aberration multigraph.

    An incomplete aberration multigraph (IAMG) encodes a partial rejoining of DSB
    ends on a fixed chromosomal backbone. Unlike a complete
    :class:`~aberration_multigraph.amg.AberrationMultigraph`, the rejoin edges do not
    form a perfect matching.

    This class is primarily used as an intermediate object during AMG generation
    and backtracking algorithms.
    """
    def __init__(self, chromatins, dsbs, rejoins, name=''):
        """
        Initialize an incomplete aberration multigraph.

        Parameters
        ----------
        chromatins : iterable of tuple
            Chromatin edges encoding chromosome structure.
        dsbs : iterable of tuple
            DSB edges defining breakpoints.
        rejoins : iterable of tuple
            Partial set of rejoin edges.
        name : str
            Identifier for the incomplete AMG.
        """
        super().__init__(chromatins, dsbs, rejoins, name)
        rejoin_verts = set([u for (u,v) in self.rejoins]+
                           [v for (u,v) in self.rejoins])
        self.free = [v for v in self.graph.nodes
                        if self.graph.degree[v] > 1 and v not in rejoin_verts]
        self.count = None
        if len(self.free) % 2 != 0:
            raise ValueError('There are an odd number of free vertices!')
    
    def complete_amgs(self):
        """
        Generate all complete aberration multigraphs extending this incomplete AMG.

        Returns
        -------
        generator of AberrationMultigraph
            Generator over all possible completions obtained by pairing the remaining
            free vertices into valid rejoin edges.

        Side Effects
        ------------
        Sets the attribute ``self.count``, which tracks the number of
        completed AMGs generated.
        """
        self.count = 0
        priorities = {v: len(self.free)-1 for v in self.free}
        for (a,b) in self.dsbs:
            if a in self.free and b in self.free:
                priorities[a] -= 1
                priorities[b] -= 1
        free_verts = [(priorities[v], v) for v in priorities]
        hq.heapify(free_verts)
        return self._gen_rejoins([], free_verts)
    
    def count_amgs(self):
        """
        Count the number of complete aberration multigraphs extending this incomplete AMG.

        Returns
        -------
        int
            Number of complete AMGs consistent with the current partial rejoining.

        Notes
        -----
        The count is cached after the first computation and stored in ``self.count``.
        """
        if self.count is not None:
            return self.count
        priorities = {v: len(self.free)-1 for v in self.free}
        for (a,b) in self.dsbs:
            if a in self.free and b in self.free:
                priorities[a] -= 1
                priorities[b] -= 1
        free_verts = [(priorities[v], v) for v in priorities]
        hq.heapify(free_verts)
        self.count = self._count_rejoins([], free_verts)
        return self.count
        
    def _count_rejoins(self, rejoins, free_verts):
        """
        Recursively count valid completions of the rejoin set.

        Parameters
        ----------
        rejoins : list of tuple
            Rejoin edges fixed so far in the recursion.
        free_verts : list
            Heap of currently unmatched vertices, prioritized by remaining flexibility.

        Returns
        -------
        int
            Number of valid completions extending the current partial rejoining.

        Notes
        -----
        This is an internal helper method and not part of the public API.
        """
        loc_count = 0
        _, v = hq.heappop(free_verts)
        if len(free_verts) == 1:
            return loc_count+1
        else:
            for _, w in free_verts:
                if (v,w) not in self.dsbs and (w,v) not in self.dsbs:
                    new_free_verts = self._remaining_vertices(free_verts, v, w)
                    new_rejoins = rejoins + [(v, w)]
                    loc_count += self._count_rejoins(new_rejoins, new_free_verts)
        return loc_count
        
    def _gen_rejoins(self, rejoins, free_verts):
        """
        Recursively generate valid rejoin completions.

        Parameters
        ----------
        rejoins : list of tuple
            Rejoin edges fixed so far in the recursion.
        free_verts : list
            Heap of currently unmatched vertices.

        Yields
        ------
        AberrationMultigraph
            A complete aberration multigraph obtained by extending the current
            incomplete AMG.

        Notes
        -----
        This is an internal helper method used by :meth:`complete_amgs`.
        """
        _, v = hq.heappop(free_verts)
        # If there is only one unmatched vertex left, pair it to this vertex 
        # and generate an AMG.
        # Else, pair this vertex with all remaining vertices and generate AMGs
        # recursively.
        if len(free_verts) == 1:
            _, u = hq.heappop(free_verts)
            self.count += 1
            yield AberrationMultigraph(self.chromatins,
                                       self.dsbs,
                                       list(self.rejoins)+rejoins+[(v,u)],
                                       self.name+'_'+str(self.count))
        else:
            for _, w in free_verts:
                if (v,w) not in self.dsbs and (w,v) not in self.dsbs:
                    new_free_verts = self._remaining_vertices(free_verts, v, w)
                    new_rejoins = rejoins + [(v, w)]
                    yield from self._gen_rejoins(new_rejoins, new_free_verts)

    def _remaining_vertices(self, free_verts, v, w):
        """
        Update the heap of free vertices after pairing two vertices.

        Parameters
        ----------
        free_verts : list
            Heap of currently unmatched vertices.
        v, w : int
            Vertices paired at the current recursion step.

        Returns
        -------
        list
            Updated heap of free vertices with adjusted priorities.

        Notes
        -----
        Priority updates reflect reduced pairing options induced by DSB constraints.
        """
        new_free_verts = []
        for priority, u in free_verts:
            if ((u,v) in self.dsbs or (v,u) in self.dsbs or 
                (u,w) in self.dsbs or (w,u) in self.dsbs):
                new_free_verts.append((priority-1, u))
            # if (self.dsb_pair[u] == v
            #         or self.dsb_pair[u] == w):
            #     new_free_verts.append((priority-1, u))
            elif u != w:
                new_free_verts.append((priority-2, u))
        hq.heapify(new_free_verts)
        return new_free_verts

    # def save_all_complete(self, filename):
    #     pass


if __name__ == '__main__':
    chromatin = [(1,2), (3,4), (5,6), (7,8), (9,10), (11,12), (13,14), (15,16), (17,18), (19,20), (21,22), (23,24), (25,26)]
    dsb = [(2,3), (4,5), (8,9), (10,11), (12,13), (14,15), (18,19), (20,21), (22,23), (24,25)]
    rejoin = [(3,5), (9,10), (13,14), (19,20), (23,24)]
    inc = IncompleteAMG(chromatin, dsb, rejoin)
    print(inc.count_amgs())
    # for amg in inc.complete_amgs():
    #     print(amg.rejoins)
