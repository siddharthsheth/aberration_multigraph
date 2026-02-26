import os
from collections import defaultdict
from nihms_patient import NIHMSPatient
from sv_utils import BreakLocation

def cycle_structure_str(cs):
    cs_str = ''
    for length in cs:
        cs_str += f'{cs[length]}C{length//2}+'
    return cs_str[:-1]


dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path+'/../../data/')
data_file = 'P05-1657.pickle'
    
patient = NIHMSPatient.load_from_file(data_file, 'nihms_files_backup/')

subset_twist_edge_pairs = [((4,), None),
                           ((7,), ((BreakLocation(chrom=7, bp=3581948), BreakLocation(chrom=7, bp=46822440)), )),
                           ((8,12), ((BreakLocation(chrom=8, bp=40932539), BreakLocation(chrom=8, bp=40938607)),
                                     (BreakLocation(chrom=12, bp=0), BreakLocation(chrom=12, bp=28921594)))),
                           ((21,), None)]

twist_edges = [(BreakLocation(chrom=7, bp=3581948), BreakLocation(chrom=7, bp=46822440))]
                # (BreakLocation(chrom=8, bp=40932539), BreakLocation(chrom=8, bp=40938607)),
                # (BreakLocation(chrom=12, bp=0), BreakLocation(chrom=12, bp=28921594))]
output = {}
for subset, twist_edges in subset_twist_edge_pairs:
    edge_pairs = {}
    if twist_edges:
        amgs = list(patient.amg(subset).complete_amgs())
        amg_names = {amg: i for i, amg in enumerate(amgs)}
        amg_cs = {}
        for edge in twist_edges:
            related_pairs = set()
            related_cs = defaultdict(int)
            for amg in amgs:
                flip_amg = amg.edge_reversal(edge)
                if flip_amg in amg_names and amg_names[flip_amg]!= amg_names[amg]:
                    # related_pairs.add(tuple(sorted([amg_names[amg], amg_names[flip_amg]])))
                    orig_cs = cycle_structure_str(amg.cycle_structure())
                    flip_cs = cycle_structure_str(flip_amg.cycle_structure())                
                    related_cs[(orig_cs, flip_cs)] += 1
            edge_pairs[edge] = (related_pairs, related_cs)
    output[subset] = edge_pairs

for key, value in output.items():
    print(f'{key}: {value}')