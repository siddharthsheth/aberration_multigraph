from aberration_multigraph.amg import AberrationMultigraph
import heapq as hq

class IncompleteAMG(AberrationMultigraph):
    def __init__(self, chromatin_edges, dsb_edges, rejoin_edges, name=''):
        super().__init__(chromatin_edges, dsb_edges, rejoin_edges, name)
        rejoin_verts = set([u for (u,v) in self.rejoin_edges]+
                           [v for (u,v) in self.rejoin_edges])
        self.free = [v for v in self.graph.nodes
                        if self.graph.degree[v] > 1 and v not in rejoin_verts]
        self.count = None
        if len(self.free) % 2 != 0:
            print('There are an odd number of free vertices!')
            raise ValueError
    
    def complete_amgs(self):
        self.count = 0
        priorities = {v: len(self.free)-1 for v in self.free}
        for (a,b) in self.dsb_edges:
            if a in self.free and b in self.free:
                priorities[a] -= 1
                priorities[b] -= 1
        free_verts = [(priorities[v], v) for v in priorities]
        hq.heapify(free_verts)
        return self._gen_rejoins([], free_verts)
    
    def count_amgs(self):
        if self.count is not None:
            return self.count
        priorities = {v: len(self.free)-1 for v in self.free}
        for (a,b) in self.dsb_edges:
            if a in self.free and b in self.free:
                priorities[a] -= 1
                priorities[b] -= 1
        free_verts = [(priorities[v], v) for v in priorities]
        hq.heapify(free_verts)
        self.count = self._count_rejoins([], free_verts)
        return self.count
        
    def _count_rejoins(self, rejoins, free_verts):
        loc_count = 0
        _, v = hq.heappop(free_verts)
        if len(free_verts) == 1:
            return loc_count+1
        else:
            for _, w in free_verts:
                if (v,w) not in self.dsb_edges and (w,v) not in self.dsb_edges:
                    new_free_verts = self._remaining_vertices(free_verts, v, w)
                    new_rejoins = rejoins + [(v, w)]
                    loc_count += self._count_rejoins(new_rejoins, new_free_verts)
        return loc_count
        
    def _gen_rejoins(self, rejoins, free_verts):
        _, v = hq.heappop(free_verts)
        # If there is only one unmatched vertex left, pair it to this vertex 
        # and generate an AMG.
        # Else, pair this vertex with all remaining vertices and generate AMGs
        # recursively.
        if len(free_verts) == 1:
            _, u = hq.heappop(free_verts)
            self.count += 1
            yield AberrationMultigraph(self.chromatin_edges,
                                       self.dsb_edges,
                                       list(self.rejoin_edges)+rejoins+[(v,u)],
                                       self.name+'_'+str(self.count))
        else:
            for _, w in free_verts:
                if (v,w) not in self.dsb_edges and (w,v) not in self.dsb_edges:
                    new_free_verts = self._remaining_vertices(free_verts, v, w)
                    new_rejoins = rejoins + [(v, w)]
                    yield from self._gen_rejoins(new_rejoins, new_free_verts)

    def _remaining_vertices(self, free_verts, v, w):
        new_free_verts = []
        for priority, u in free_verts:
            if ((u,v) in self.dsb_edges or (v,u) in self.dsb_edges or 
                (u,w) in self.dsb_edges or (w,u) in self.dsb_edges):
                new_free_verts.append((priority-1, u))
            # if (self.dsb_pair[u] == v
            #         or self.dsb_pair[u] == w):
            #     new_free_verts.append((priority-1, u))
            elif u != w:
                new_free_verts.append((priority-2, u))
        hq.heapify(new_free_verts)
        return new_free_verts

    def save_all_complete(self, filename):
        pass


if __name__ == '__main__':
    chromatin = [(1,2), (3,4), (5,6), (7,8), (9,10), (11,12), (13,14), (15,16), (17,18), (19,20), (21,22), (23,24), (25,26)]
    dsb = [(2,3), (4,5), (8,9), (10,11), (12,13), (14,15), (18,19), (20,21), (22,23), (24,25)]
    rejoin = [(3,5), (9,10), (13,14), (19,20), (23,24)]
    inc = IncompleteAMG(chromatin, dsb, rejoin)
    print(inc.count_amgs())
    # for amg in inc.complete_amgs():
    #     print(amg.rejoin_edges)
