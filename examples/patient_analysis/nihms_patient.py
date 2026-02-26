from sv_utils import BreakLocation
from aberration_multigraph.incomplete_amg import IncompleteAMG
import os
import pickle

class NIHMSPatient:
    def __init__(self, id, breakpoints=(), rejoins=()):
        self.id = id
        self.breakpoints = list(breakpoints)
        self.rejoins = list(rejoins)
        self._vertices = {}
        self._amg = None

    def add_breakpoint(self, breakpoint):
        self.breakpoints.append(breakpoint)
    
    def add_rejoin(self, rejoin):
        self.rejoins.append(rejoin)

    def amg(self, subset):
        self._generate_vertices(subset)
        suffix = '_'.join(str(x) for x in subset)

        # generate incomplete amg
        self._amg = IncompleteAMG(self._chromatin_edges(subset),
                                  self._dsb_edges(subset),
                                  self._rejoin_edges(subset),
                                  self.id+'_'+suffix)
        return self._amg

    def _generate_vertices(self, subset):
        for chrom, bp in self.breakpoints:
            if chrom not in subset:
                continue
            self._vertices[(chrom, bp)] = BreakLocation(chrom, bp)
            self._vertices[(chrom, bp+1)] = BreakLocation(chrom, bp+1)
    
    def _dsb_edges(self, subset):
        dsbs = []
        for chrom, bp in self.breakpoints:
            if chrom not in subset:
                continue
            dsbs.append((self._vertices[(chrom, bp)],
                         self._vertices[(chrom, bp+1)]))
        return dsbs
    
    def _chromatin_edges(self, subset):
        chroms = []
        prev_chrom = prev_bp = None
        for chrom, bp in self.breakpoints:
            if chrom not in subset:
                continue
            if prev_chrom is None:
                prev_chrom = 0
            if prev_chrom != chrom:
                if prev_bp is not None:
                    chroms.append((self._vertices[(prev_chrom, prev_bp+1)],
                                    BreakLocation(prev_chrom, float('inf'))))
                chroms.append((BreakLocation(chrom, 0),
                               self._vertices[(chrom, bp)]))
            else:
                chroms.append((self._vertices[(prev_chrom, prev_bp+1)],
                               self._vertices[(chrom, bp)]))
            prev_chrom, prev_bp = chrom, bp
        if prev_chrom != None:
            chroms.append((self._vertices[(prev_chrom, prev_bp+1)],
                            BreakLocation(prev_chrom, float('inf'))))
        return chroms

    def _rejoin_edges(self, subset):
        rejoin_edges = []
        for r1, r2 in self.rejoins:
            c1, bp1, loc1 = r1
            c2, bp2, loc2 = r2
            if c1 not in subset or c2 not in subset:
                continue
            if loc1 == '+':
                bp1 += 1
            if loc2 == '+':
                bp2 += 1
            rejoin_edges.append((self._vertices[(c1, bp1)],
                                    self._vertices[(c2, bp2)]))
        return rejoin_edges

    def save_to_file(self, filename='', path=''):
        if filename == '':
            filename = str(self.id)
        with open(path+filename+'.pickle', 'wb') as file:
            pickle.dump(self, file)

    @staticmethod
    def load_from_file(filename, path=''):
        file = open(path+filename, 'rb')
        return pickle.load(file)

if __name__ == '__main__':
    dir_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(dir_path+'/../../data')
    data_file = 'nihms.csv'
    patients = {}
    with open(data_file) as file:
        headers = file.readline()
        for entry in file:
            entry = entry.split(',')
            pat_id = entry[0]
            patient = patients.setdefault(pat_id, NIHMSPatient(pat_id))
            
            patient.add_breakpoint((int(entry[2]), int(entry[4])))
            patient.add_breakpoint((int(entry[5]), int(entry[7])))

            r1 = (int(entry[2]), int(entry[4]), entry[3])
            r2 = (int(entry[5]), int(entry[7]), entry[6])
            patient.add_rejoin((r1, r2))
    print(os.getcwd())
    if not os.path.isdir('nihms_patient_files/'):
        os.mkdir('nihms_patient_files/')
    for pat_id, patient in patients.items():
        print(f"Saving patient {pat_id}")
        patient.breakpoints.sort()
        patient.save_to_file(path='nihms_patient_files/')
        
    patients = [NIHMSPatient.load_from_file(file, 'nihms_patient_files/')
                            for file in os.listdir('nihms_patient_files/')]
    print(len(patients))
    print(patients[0].id)