from collections import defaultdict, namedtuple
import os

class SVVertex:
    def __init__(self, chrom, bp, dsb=None, rejoin=None, chromatin=None, digits_truncate=4):
        self.chrom = chrom
        self.bp = bp
        self.dsb = dsb
        self.rejoin = rejoin
        self.chromatin = chromatin
        self.bp_digits_truncate = digits_truncate
        self.occurrences = 1

    def amg_vertex(self):
        return (self.chrom, self.bp)
    
    def __str__(self) -> str:
        return f'({self.chrom}, {self.bp%10**self.bp_digits_truncate})'
    
    def __hash__(self):
        return hash((self.chrom, self.bp))
    
    def __eq__(self, other):
        return True if self.chrom == other.chrom and self.bp == other.bp else False

class StructuralVariation:
    def __init__(self, c1, s1, e1, c2, s2, e2, sv_id, pe_support, str1, str2, svclass, svmethod):
        c1 = 23 if c1 == 'X' else c1
        c1 = 24 if c1 == 'Y' else c1
        c2 = 23 if c2 == 'X' else c2
        c2 = 24 if c2 == 'Y' else c2
        self.chrom1= int(c1)
        self.start1 = int(s1)
        self.end1 = int(e1)
        self.chrom2 = int(c2)
        self.start2 = int(s2)
        self.end2 = int(e2)
        self.sv_id = sv_id
        self.pe_support = pe_support
        self.strand1 = str1
        self.strand2 = str2
        self.svclass = svclass
        self.svmethod = svmethod

    def __str__(self):
        return str((self.chrom1, self.start1, self.end1, self.chrom2, self.start2, self.end2))

BreakLocation = namedtuple('BreakLocation', ['chrom', 'bp'])
chromothripsis = '72f0a49a-aec8-47e5-846a-956c4da1507c.pcawg_consensus_1.6.161116.somatic.sv'
simple = 'e1217ebe-1826-41a9-b6c4-702100a66f5e.pcawg_consensus_1.6.161116.somatic.sv'
medium = '0ae2193f-0d68-485a-b8c2-7568cbcce33e.pcawg_consensus_1.6.161116.somatic.sv'
# file = simple
# file += '.bedpe'
# path = '/Users/siddharthsheth/Dropbox/work/research-projects/TQFT Cancer Progression/tcga/open/'
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path+'/..')
data_path = 'data/tcga/open/'
log_file = 'svlog.log'
logs = []
for i, file in enumerate(os.listdir(data_path)):
    extension = file.split('.')[-1]
    if extension != 'bedpe':
        print(f'Skipped file number {i+1}: {file}')
        continue
    print(f'Working on file number {i+1}: {file}.')
    old_logs = len(logs)
    logs.append(f'File: {file}\n')
    with open(data_path+file, 'r') as svfile:
        headers = svfile.readline()
        svs = []
        bl_to_sv = defaultdict(dict)
        for entry in svfile:
            entry = entry.split()
            sv = StructuralVariation(*entry)
            svs.append(sv)
            entry[0] = 23 if entry[0] == 'X' else entry[0]
            entry[0] = 24 if entry[0] == 'Y' else entry[0]
            entry[3] = 23 if entry[3] == 'X' else entry[3]
            entry[3] = 24 if entry[3] == 'Y' else entry[3]
            u_1, v_1, u_2, v_2 = BreakLocation(int(entry[0]), int(entry[1])), BreakLocation(int(entry[0]), int(entry[2])), BreakLocation(int(entry[3]), int(entry[4])), BreakLocation(int(entry[3]), int(entry[5]))
            if u_1 in bl_to_sv:
                bl_to_sv[u_1]['svs'].append((sv, 1))
            else:
                # create new SVVertex for this location
                bl_to_sv[u_1]['vertex'] = SVVertex(*u_1)
                # update the StructuralVariation and the position in that SV this location belongs to
                bl_to_sv[u_1]['svs'] = [(sv, 1)]
            if v_1 in bl_to_sv:
                bl_to_sv[v_1]['svs'].append((sv, 2))
            else:
                # create new SVVertex for this location
                bl_to_sv[v_1]['vertex'] = SVVertex(*v_1)
                # update the StructuralVariation and the position in that SV this location belongs to
                bl_to_sv[v_1]['svs'] = [(sv, 2)]
            if u_2 in bl_to_sv:
                bl_to_sv[u_2]['svs'].append((sv, 3))
            else:
                # create new SVVertex for this location
                bl_to_sv[u_2]['vertex'] = SVVertex(*u_2)
                # update the StructuralVariation and the position in that SV this location belongs to
                bl_to_sv[u_2]['svs'] = [(sv, 3)]
            if v_2 in bl_to_sv:
                bl_to_sv[v_2]['svs'].append((sv, 4))
            else:
                # create new SVVertex for this location
                bl_to_sv[v_2]['vertex'] = SVVertex(*v_2)
                # update the StructuralVariation and the position in that SV this location belongs to
                bl_to_sv[v_2]['svs'] = [(sv, 4)]
            
        for bp_loc in bl_to_sv:
            if len(bl_to_sv[bp_loc]['svs']) > 2:
                logs.append(f"{bp_loc} appears in {len(bl_to_sv[bp_loc]['svs'])} SVs.\n")
        
        
        # for bl in bl_to_sv:
        #     print(f'{bl}: {bl_to_sv[bl]}')

        if len(logs)-old_logs == 1:
            logs.append(f'No anomalies detected.\n----------------------\n')
        else:
            logs.append('----------------------\n\n')

with open(log_file, 'a') as log:
    log.writelines(logs)