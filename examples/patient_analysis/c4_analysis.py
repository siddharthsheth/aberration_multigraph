from aberration_multigraph.generator import AMGGenerator
from aberration_multigraph.amg import AberrationMultigraph

# want those AMGs in (1,4) whose CS remains C_4 throughout

def cycle_structure_str(cs):
    cs_str = ''
    for length in cs:
        cs_str += f'{cs[length]}C{length//2}+'
    return cs_str[:-1]



# twist_edges = [(2,3), (4,5), (6,7), (8,9)]
twist_edges = [[(2,3), (4,5), (6,7), (2,3), (4,5), (6,7)], 
                [(4,5), (2,3), (6,7), (4,5), (2,3), (6,7)],
                [(2,3), (4,5), (6,7), (6,7), (4,5), (2,3)],
                [(4,5), (2,3), (6,7), (6,7), (2,3), (4,5)]]
# twist_edges = [(2,3), (4,5)]

output = {}
for i in range(1, len(twist_edges)+1):
    for edges in twist_edges:
        edge_pairs = {}
        amg_gen = AMGGenerator(1, (4,))
        constant_amgs = []
        for amg in amg_gen.generate_amgs():
            is_const = True
            orig_cs = cycle_structure_str(amg.cycle_structure())
            flip_amg = AberrationMultigraph(amg.chromatin_edges, amg.dsb_edges, amg.rejoin_edges)
            for edge in edges:
                flip_amg = flip_amg.edge_reversal(edge)
                if(cycle_structure_str(flip_amg.cycle_structure()) != orig_cs):
                    is_const = False
                    break
            if is_const:
                constant_amgs.append(amg)
        output[tuple(edges)] = constant_amgs

print(i)
for x in output:
    print(x)
    for amg in output[x]:
        print(amg.chromatin_edges, amg.dsb_edges, amg.rejoin_edges)
# print(output)

# ['graph', 'chromatin_edges', 'dsb_edges', 'rejoin_edges', 'name', 'num_chromosome', '__module__', '__doc__', '__init__', 'diameter', 'girth', 'cycles', 'cycle_structure', 'is_connected', 'draw', 'chromatin_edge_twist', 'total_twist', '_chromatin_edges_total_twist', '_dsb_edges_total_twist', '_rejoin_edges_total_twist', 'total_swap', '_chromatin_edges_total_swap', '_dsb_edges_total_swap', '_rejoin_edges_total_swap', 'save_to_file', 'load_from_file', '__hash__', '__eq__', '__dict__', '__weakref__', '__repr__', '__str__', '__getattribute__', '__setattr__', '__delattr__', '__lt__', '__le__', '__ne__', '__gt__', '__ge__', '__new__', '__reduce_ex__', '__reduce__', '__subclasshook__', '__init_subclass__', '__format__', '__sizeof__', '__dir__', '__class__']