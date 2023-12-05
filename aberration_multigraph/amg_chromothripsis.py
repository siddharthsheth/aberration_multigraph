from aberration_multigraph.amg import AberrationMultigraph
from collections import defaultdict, namedtuple
from os import chdir
from numpy import inf
from matplotlib import pyplot as plt

bp_digits_truncate = 12

class SVVertex:
    def __init__(self, chrom, bp, dsb=None, rejoin=None, chromatin=None, digits_truncate=4):
        self.chrom = chrom
        self.bp = bp
        self.dsb = dsb
        self.rejoin = rejoin
        self.chromatin = chromatin
        self.bp_digits_truncate = digits_truncate

    def amg_vertex(self):
        return (self.chrom, self.bp)
    
    def __str__(self) -> str:
        return f'({self.chrom}, {self.bp%10**self.bp_digits_truncate})'
    
    def __hash__(self):
        return hash((self.chrom, self.bp))

# SVVertex = namedtuple('SVVertex', ['chrom', 'bp', 'dsb_edge', 'rejoin_edge', 'chromatin_edge'], defaults=[None, None, None])

def amg_from_pcawg(file, complete_strat=0, chromatin_strat=0):
    dsb_edges, rejoin_edges, chromatin_edges = set(), set(), set()
    chrom_to_dsb_map, vertex_to_dsb_map, vertex_to_rejoin_map = defaultdict(set), {}, {}
    
    with open(file, 'r') as amg_file:
        headers = amg_file.readline()
        for entry in amg_file:
            entry = entry.split()
            entry[0] = 23 if entry[0] == 'X' else entry[0]
            entry[0] = 24 if entry[0] == 'Y' else entry[0]
            entry[3] = 23 if entry[3] == 'X' else entry[3]
            entry[3] = 24 if entry[3] == 'Y' else entry[3]
            # u_1, v_1 = (int(entry[0]), int(entry[1])%10**bp_digits_truncate), (int(entry[0]), int(entry[2])%10**bp_digits_truncate)
            # u_2, v_2 = (int(entry[3]), int(entry[4])%10**bp_digits_truncate), (int(entry[3]), int(entry[5])%10**bp_digits_truncate)
            u_1, v_1, u_2, v_2 = SVVertex(int(entry[0]), int(entry[1])), SVVertex(int(entry[0]), int(entry[2])), SVVertex(int(entry[3]), int(entry[4])), SVVertex(int(entry[3]), int(entry[5]))
            # u_1 = SVVertex(int(entry[0]), int(entry[1]))
            # v_1 = SVVertex(int(entry[0]), int(entry[2]))
            # u_2 = SVVertex(int(entry[3]), int(entry[4]))
            # v_2 = SVVertex(int(entry[3]), int(entry[5]))

            # store locations of DSBs per chromosome
            c_1, c_2 = u_1.chrom, u_2.chrom
            chrom_to_dsb_map[c_1].update({u_1, v_1})
            chrom_to_dsb_map[c_2].update({u_2, v_2})

            # chrom_to_dsb_map[u_1[0]].add(u_1[1])
            # chrom_to_dsb_map[v_1[0]].add(v_1[1])
            # chrom_to_dsb_map[u_2[0]].add(u_2[1])
            # chrom_to_dsb_map[v_2[0]].add(v_2[1])

            
            # add dsb edges in this entry
            dsb_edges.add((u_1, v_1))
            dsb_edges.add((u_2, v_2))
            u_1.dsb = v_1
            v_1.dsb = u_1
            u_2.dsb = v_2
            v_2.dsb = u_2
            # vertex_to_dsb_map[u_1] = (u_1, v_1)
            # vertex_to_dsb_map[v_1] = (u_1, v_1)
            # vertex_to_dsb_map[u_2] = (u_2, v_2)
            # vertex_to_dsb_map[v_2] = (u_2, v_2)

            # add rejoin edges in this entry
            if entry[8] == '-' and entry[9] == '-':
                rejoin_edges.add((u_1, u_2))
                u_1.rejoin = u_2
                u_2.rejoin = u_1
            elif entry[8] == '-' and entry[9] == '+':
                rejoin_edges.add((u_1, v_2))
                u_1.rejoin = v_2
                v_2.rejoin = u_1
            elif entry[8] == '+' and entry[9] == '-':
                rejoin_edges.add((v_1, u_2))
                v_1.rejoin = u_2
                u_2.rejoin = v_1
            else:
                rejoin_edges.add((v_1, v_2))
                v_1.rejoin = v_2
                v_2.rejoin = v_1
            # if entry[8] == '-' and entry[9] == '-':
            #     rejoin_edges.add((u_1, u_2))
            #     vertex_to_rejoin_map[u_1] = (u_1, u_2)
            #     vertex_to_rejoin_map[u_2] = (u_1, u_2)
            # elif entry[8] == '-' and entry[9] == '+':
            #     rejoin_edges.add((u_1, v_2))
            #     vertex_to_rejoin_map[u_1] = (u_1, v_2)
            #     vertex_to_rejoin_map[v_2] = (u_1, v_2)
            # elif entry[8] == '+' and entry[9] == '-':
            #     rejoin_edges.add((v_1, u_2))
            #     vertex_to_rejoin_map[v_1] = (v_1, u_2)
            #     vertex_to_rejoin_map[u_2] = (v_1, u_2)
            # else:
            #     rejoin_edges.add((v_1, v_2))
            #     vertex_to_rejoin_map[v_1] = (v_1, v_2)
            #     vertex_to_rejoin_map[v_2] = (v_1, v_2)
        
        # print('Intermediate Stage')
        # print('DSB Edges')
        # print(dsb_edges)

        # print('Rejoin Edges')
        # print(rejoin_edges)

        # CREATE CHROMATIN EDGES
        for chrom in chrom_to_dsb_map:
            chrom_vertices = sorted(list(chrom_to_dsb_map[chrom]), key= lambda x: x.bp)
            v0, vn = SVVertex(chrom, 0), SVVertex(chrom, inf)
            chromatin_edges.add((v0, chrom_vertices[0]))
            chrom_vertices[0].chromatin = v0
            i = 1
            # print('stop at ' + str(len(chrom_vertices)))
            if chrom == 11:
                print([str(c) for c in chrom_vertices])
            while i < len(chrom_vertices)-1:
                u, v, w = chrom_vertices[i-1], chrom_vertices[i], chrom_vertices[i+1]
                if (u, v) not in dsb_edges:
                    print(f'Catastrophic error! {(str(u),str(v))} is not a DSB edge at i={i}! Also, w={str(w)}.')
                    break
                # print(chrom, i, u, v, w)
                if (v, w) not in dsb_edges:
                    chromatin_edges.add((v, w))
                    v.chromatin = w
                    w.chromatin = v
                    i += 2
                else:
                    # what if there are two successive DSB edges?
                    # if the intermediate vertex is free, then merge the two vertices
                    if v.rejoin is None:
                        print(f'Merging vertices in chromosome {chrom}.')
                        dsb_edges.remove((u, v))
                        dsb_edges.remove((v, w))
                        dsb_edges.add((u, w))
                        chrom_to_dsb_map[chrom].discard(v)
                        u.dsb = w
                        w.dsb = u
                        v.dsb = None
                    else:
                        # option 1: merge them nevertheless?
                        if chromatin_strat == 0:
                            # make dsb edge between u and w
                            dsb_edges.remove((u, v))
                            dsb_edges.remove((v, w))
                            dsb_edges.add((u, w))

                            # discard v as a DSB location in chrom
                            chrom_to_dsb_map[chrom].discard(v)
                            u.dsb = w
                            w.dsb = u
                            v.dsb = None
                            
                            # update vertex in rejoin edge incident to v
                            a = v.rejoin
                            if (a,v) in rejoin_edges or (v,a) in rejoin_edges:
                                rejoin_edges.discard((a,v))
                                rejoin_edges.discard((v,a))
                                if u.rejoin is not None:
                                    rejoin_edges.add((w,a) if (v,a) in rejoin_edges else (a,w))
                                else:
                                    rejoin_edges.add((u,a) if (v,a) in rejoin_edges else (a,u))

                        # option 2: add an imaginary chromatin edge at the common vertex?
                        elif chromatin_strat == 1:
                            # v_2 = (chrom, chrom_vertices[i]+0.5)                  # imaginary vertex
                            v_2 = SVVertex(chrom, chrom_vertices[i].bp+0.5)
                            chrom_to_dsb_map[chrom].add(v_2[1])
                            chromatin_edges.add((v, v_2))               # imaginary chromatin edge
                            # a, b = vertex_to_rejoin_map[v]
                            # if (a,b) in rejoin_edges:
                            #     rejoin_edges.remove((a,b))
                            #     if u in vertex_to_rejoin_map:
                            #         rejoin_edges.add((v_2,b) if v==a else (a,v_2))
                            a = v.rejoin
                            if (a,v) in rejoin_edges or (v,a) in rejoin_edges:
                                rejoin_edges.difference_update({(a,v), (v,a)})
                                if u.rejoin is not None:
                                    rejoin_edges.add((v_2,a) if (v,a) in rejoin_edges else (a,v_2))
                            else:
                                print('Something suspicious')

                        # option 3: ignore it?
                        else:
                            pass

                    i += 1
            
            chromatin_edges.add((chrom_vertices[-1], vn))


    # COMPLETE THE AMG BY ADDING MORE REJOIN EDGES
    # if complete_strat == 0:
    #     # complete the AMG by pairing complement vertices in each rejoining
    #     unmatched_dsbs = {i for i in dsb_edges}
    #     for (u,v) in dsb_edges:
    #         if (u,v) not in unmatched_dsbs:
    #             continue
    #         left, right = u, v
    #         while left in vertex_to_rejoin_map:
    #             a, b = vertex_to_rejoin_map[left]
    #             if left == a:
    #                 c, d = vertex_to_dsb_map[b]
    #                 left = d if b == c else c
    #                 unmatched_dsbs.discard((c,d))
    #             else:
    #                 c, d = vertex_to_dsb_map[a]
    #                 left = d if a == c else c
    #                 unmatched_dsbs.discard((c,d))
    #         while right in vertex_to_rejoin_map:
    #             a, b = vertex_to_rejoin_map[right]
    #             if right == a:
    #                 c, d = vertex_to_dsb_map[b]
    #                 right = d if b == c else c
    #                 unmatched_dsbs.discard((c,d))
    #             else:
    #                 c, d = vertex_to_dsb_map[a]
    #                 right = d if a == c else c
    #                 unmatched_dsbs.discard((c,d))
    #         rejoin_edges.add((left, right))


    # print('Final Stage')
    # print('DSB Edges')
    # print(dsb_edges)

    # print('Rejoin Edges')
    # print(rejoin_edges)

    # print('Chromosome Edges')
    # print(chromatin_edges)

    return AberrationMultigraph({(u.amg_vertex(), v.amg_vertex()) for (u, v) in chromatin_edges}, {(u.amg_vertex(), v.amg_vertex()) for (u,v) in dsb_edges}, {(u.amg_vertex(), v.amg_vertex()) for (u,v) in rejoin_edges})
    
if __name__ == '__main__':
    chromothripsis = '72f0a49a-aec8-47e5-846a-956c4da1507c.pcawg_consensus_1.6.161116.somatic.sv'
    simple = 'e1217ebe-1826-41a9-b6c4-702100a66f5e.pcawg_consensus_1.6.161116.somatic.sv'
    medium = '0ae2193f-0d68-485a-b8c2-7568cbcce33e.pcawg_consensus_1.6.161116.somatic.sv'
    file = medium
    file += '.txt'
    path = '/Users/siddharthsheth/Dropbox/work/research-projects/TQFT Cancer Progression/tcga/open/'
    for strat in range(3):
        print(f'Strategy: {strat}')
        amg = amg_from_pcawg(path+file, chromatin_strat=strat)
        plt.subplot(1,3,strat+1)
        amg.draw()
        print(amg.diameter())
        print(amg.cycle_structure())
        # print(amg.dsb_edges)
        # print(amg.rejoin_edges)
        # print(amg.chromatin_edges)
        print('\n')
    plt.show()