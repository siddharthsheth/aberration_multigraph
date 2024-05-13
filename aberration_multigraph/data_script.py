import os
import pickle
from collections import Counter
from aberration_multigraph.amg import AberrationMultigraph
from aberration_multigraph.nihms_patient import NIHMSPatient

dir_path = os.path.dirname(os.path.realpath(__file__))

os.chdir(dir_path+'/../data/nihms_patient_files/')
patient = NIHMSPatient.load_from_file('P05-1657.pickle')
subsets = [{4}, {7}, {8,12}, {21}]
for subset in subsets:
    print(subset)
    inc_amg = patient.amg(subset)

    for amg in inc_amg.complete_amgs():
        amg.save_to_file(path='../nihms_amg/P05-1657/')

    print(inc_amg.count_amgs())

os.chdir(dir_path+'/../data/nihms_amg/P05-1657/')
diams = Counter()
css = Counter()
for file in os.listdir():
    amg = AberrationMultigraph.load_from_file(file)
    diam, cs = amg.diameter(), amg.cycle_structure()
    print(file, diam, cs)
    diams[diam] += 1
    css['+'.join(str(cs[c])+'*'+str(c) for c in cs)] += 1

print(diams)
print(css)