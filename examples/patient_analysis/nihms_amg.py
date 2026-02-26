import os
from collections import Counter
from aberration_multigraph.amg import AberrationMultigraph
from nihms_patient import NIHMSPatient

def cycle_structure_str(cs):
    cs_str = ''
    for length in cs:
        cs_str += f'{cs[length]}C{length//2}+'
    return cs_str[:-1]

PATIENT_ID = 'P05-1657'
# PATIENT_ID = 'P08-217'
SUBSETS = [{4}, {7}, {8,12}, {21}]
# SUBSETS = [{1}, {3}, {8, 12}, {10, 14}, {13}, {17}, {22}, {21}]
# SUBSETS = [{1}, {2}, {3, 19, 22}, {5}, {14}, ]

dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path+'/../../data/nihms_patient_files/')
patient = NIHMSPatient.load_from_file(f'{PATIENT_ID}.pickle')
if not os.path.isdir('../nihms_amg/'):
    os.mkdir('../nihms_amg')

if not os.path.isdir(f'../nihms_amg/{PATIENT_ID}/'):
    os.mkdir(f'../nihms_amg/{PATIENT_ID}/')

for subset in SUBSETS:
    diams = Counter()
    css = Counter()

    inc_amg = patient.amg(subset)

    for amg in inc_amg.complete_amgs():
        diam, cs = amg.diameter(), amg.cycle_structure()
        diams[diam] += 1
        css[cycle_structure_str(cs)] += 1
        # amg.save_to_file(path=f'../nihms_amg/{PATIENT_ID}/')

    print(subset)
    print(f'No of AMGs: {inc_amg.count_amgs()}')
    print("Diameter distribution:")
    for diam, count in diams.items():
        print(f'{diam}: {count}')
    print("Cycle structure distribution")
    for cs, count in css.items():
        print(f'{cs}: {count}')
    print('\n')