from aberration_multigraph.amg import AberrationMultigraph
from collections import defaultdict
from os import chdir
from numpy import inf
from matplotlib import pyplot as plt

def amg_from_pcawg(file, complete_strat=0, chromatin_strat=0):
    dsb_edges = set()
    rejoin_edges = set()
    chromatin_edges = set()
    chrom_to_dsb_map, vertex_to_dsb_map, vertex_to_rejoin_map = defaultdict(set), {}, {}

    bp_digits_truncate = 12
    with open(file, 'r') as amg_file:
        headers = amg_file.readline()
        for entry in amg_file:
            entry = entry.split()
            entry[0] = 23 if entry[0] == 'X' else entry[0]
            entry[0] = 24 if entry[0] == 'Y' else entry[0]
            entry[3] = 23 if entry[3] == 'X' else entry[3]
            entry[3] = 24 if entry[3] == 'Y' else entry[3]
            u_1, v_1 = (int(entry[0]), int(entry[1])%10**bp_digits_truncate), (int(entry[0]), int(entry[2])%10**bp_digits_truncate)
            u_2, v_2 = (int(entry[3]), int(entry[4])%10**bp_digits_truncate), (int(entry[3]), int(entry[5])%10**bp_digits_truncate)

            # store locations of DSBs per chromosome
            chrom_to_dsb_map[u_1[0]].add(u_1[1])
            chrom_to_dsb_map[v_1[0]].add(v_1[1])
            chrom_to_dsb_map[u_2[0]].add(u_2[1])
            chrom_to_dsb_map[v_2[0]].add(v_2[1])

            
            # add dsb edges in this entry
            dsb_edges.add((u_1, v_1))
            dsb_edges.add((u_2, v_2))
            vertex_to_dsb_map[u_1] = (u_1, v_1)
            vertex_to_dsb_map[v_1] = (u_1, v_1)
            vertex_to_dsb_map[u_2] = (u_2, v_2)
            vertex_to_dsb_map[v_2] = (u_2, v_2)

            # add rejoin edges in this entry
            if entry[8] == '-' and entry[9] == '-':
                rejoin_edges.add((u_1, u_2))
                vertex_to_rejoin_map[u_1] = (u_1, u_2)
                vertex_to_rejoin_map[u_2] = (u_1, u_2)
            elif entry[8] == '-' and entry[9] == '+':
                rejoin_edges.add((u_1, v_2))
                vertex_to_rejoin_map[u_1] = (u_1, v_2)
                vertex_to_rejoin_map[v_2] = (u_1, v_2)
            elif entry[8] == '+' and entry[9] == '-':
                rejoin_edges.add((v_1, u_2))
                vertex_to_rejoin_map[v_1] = (v_1, u_2)
                vertex_to_rejoin_map[u_2] = (v_1, u_2)
            else:
                rejoin_edges.add((v_1, v_2))
                vertex_to_rejoin_map[v_1] = (v_1, v_2)
                vertex_to_rejoin_map[v_2] = (v_1, v_2)
        
        # print('Intermediate Stage')
        # print('DSB Edges')
        # print(dsb_edges)

        # print('Rejoin Edges')
        # print(rejoin_edges)

        # CREATE CHROMATIN EDGES
        for chrom in chrom_to_dsb_map:
            dsbs = sorted(list(chrom_to_dsb_map[chrom]))
            chromatin_edges.add(((chrom, 0), (chrom, dsbs[0])))
            i = 1
            # print('stop at ' + str(len(dsbs)))
            while i < len(dsbs)-1:
                u, v, w = (chrom, dsbs[i-1]), (chrom, dsbs[i]), (chrom, dsbs[i+1])
                # print(chrom, i, u, v, w)
                if (v, w) not in dsb_edges:
                    chromatin_edges.add((v, w))
                    i += 2
                else:
                    # what if there are two successive DSB edges?
                    # if the intermediate vertex is free, then merge the two DSBs
                    if v not in vertex_to_rejoin_map:
                        dsb_edges.remove((u, v))
                        dsb_edges.remove((v, w))
                        dsb_edges.add((u, w))
                    else:
                        # option 1: merge them nevertheless?
                        if chromatin_strat == 0:
                            # make dsb edge between u and w
                            dsb_edges.remove((u, v))
                            dsb_edges.remove((v, w))
                            dsb_edges.add((u, w))
                            
                            # update vertex in rejoin edge incident to v
                            a, b = vertex_to_rejoin_map[v]
                            if (a,b) in rejoin_edges:
                                rejoin_edges.remove((a,b))
                                if u in vertex_to_rejoin_map:
                                    rejoin_edges.add((w,b) if v==a else (a,w))
                                else:
                                    rejoin_edges.add((u,b) if v==a else (a,u))

                        # option 2: add an imaginary chromatin edge at the common vertex?
                        elif chromatin_strat == 1:
                            v_2 = (chrom, dsbs[i]+0.5)                  # imaginary vertex
                            chromatin_edges.add((v, v_2))               # imaginary chromatin edge
                            a, b = vertex_to_rejoin_map[v]
                            if (a,b) in rejoin_edges:
                                rejoin_edges.remove((a,b))
                                if u in vertex_to_rejoin_map:
                                    rejoin_edges.add((v_2,b) if v==a else (a,v_2))

                        # option 3: ignore it?
                        else:
                            pass

                    i += 1
            chromatin_edges.add(((chrom, dsbs[-1]), (chrom, inf)))


    # COMPLETE THE AMG BY ADDING MORE REJOIN EDGES
    if complete_strat == 0:
        # complete the AMG by pairing complement vertices in each rejoining
        unmatched_dsbs = {i for i in dsb_edges}
        for (u,v) in dsb_edges:
            if (u,v) not in unmatched_dsbs:
                continue
            left, right = u, v
            while left in vertex_to_rejoin_map:
                a, b = vertex_to_rejoin_map[left]
                if left == a:
                    c, d = vertex_to_dsb_map[b]
                    left = d if b == c else c
                    unmatched_dsbs.discard((c,d))
                else:
                    c, d = vertex_to_dsb_map[a]
                    left = d if a == c else c
                    unmatched_dsbs.discard((c,d))
            while right in vertex_to_rejoin_map:
                a, b = vertex_to_rejoin_map[right]
                if right == a:
                    c, d = vertex_to_dsb_map[b]
                    right = d if b == c else c
                    unmatched_dsbs.discard((c,d))
                else:
                    c, d = vertex_to_dsb_map[a]
                    right = d if a == c else c
                    unmatched_dsbs.discard((c,d))
            rejoin_edges.add((left, right))


    # print('Final Stage')
    # print('DSB Edges')
    # print(dsb_edges)

    # print('Rejoin Edges')
    # print(rejoin_edges)

    # print('Chromosome Edges')
    # print(chromatin_edges)

    return AberrationMultigraph(chromatin_edges, dsb_edges, rejoin_edges)
    
if __name__ == '__main__':
    chromothripsis = '72f0a49a-aec8-47e5-846a-956c4da1507c.pcawg_consensus_1.6.161116.somatic.sv'
    simple = 'e1217ebe-1826-41a9-b6c4-702100a66f5e.pcawg_consensus_1.6.161116.somatic.sv'
    medium = '0ae2193f-0d68-485a-b8c2-7568cbcce33e.pcawg_consensus_1.6.161116.somatic.sv'
    file = chromothripsis
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