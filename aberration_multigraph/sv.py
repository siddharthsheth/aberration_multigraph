from collections import defaultdict, namedtuple
import os

class Patient:
    def __init__(self, id):
        self.id = id
        self.sequence = self.get_sequence()
        self.svs = self.get_svs()
        self.logs = []

    def get_sequence(self):
        file = self.id+cn_extension
        sequence = defaultdict(list)
        with open(cn_path+file, 'r') as cnfile:
            _ = cnfile.readline()                       # read headers
            for entry in cnfile:
                entry = entry.split()
                chrom = self.chrom_to_int(entry[0])
                segment = [x if x == 'NA' else int(x)  for x in entry[1:]]
                sequence[chrom].append(CNSegment(*segment))
        return sequence
            
    def get_svs(self):
        file = self.id+sv_extension
        svs = []
        with open(sv_path+file, 'r') as svfile:
            _ = svfile.readline()                       # read headers
            self.bl_to_sv = defaultdict(dict)
            for entry in svfile:
                entry = entry.split()
                entry[0] = self.chrom_to_int(entry[0])
                entry[3] = self.chrom_to_int(entry[3])
                sv = StructuralVariation(*entry)
                svs.append(sv)
                bls = (BreakLocation(entry[0], int(entry[1])),
                        BreakLocation(entry[0], int(entry[2])),
                        BreakLocation(entry[3], int(entry[4])),
                        BreakLocation(entry[3], int(entry[5]))
                    )
                for i, bl in enumerate(bls):
                    self.create_sv_vertex(bl, i+1, sv)
        return svs
    
    def chrom_to_int(self, chrom):
        if chrom == 'X':
            return 23
        elif chrom == 'Y':
            return 24
        else:
            return int(chrom)
        

    def create_sv_vertex(self, bl, pos, sv):
        # check if SVVertex has already been created
        if bl in self.bl_to_sv:
            self.bl_to_sv[bl]['svs'].append((sv, pos))
        else:
            # create new SVVertex for this BreakLocation
            self.bl_to_sv[bl]['vertex'] = SVVertex(*bl)
            # update the StructuralVariation and the position in that SV this location belongs to
            self.bl_to_sv[bl]['svs'] = [(sv, pos)]


    def check_valid(self):
        old_logs = len(self.logs)
        self.logs.append(f'File: {self.id}\n')
        
        for bp_loc in self.bl_to_sv:
            if len(self.bl_to_sv[bp_loc]['svs']) > 2:
                self.logs.append(f"{bp_loc} appears in {len(self.bl_to_sv[bp_loc]['svs'])} SVs.\n")
        
        
        # for bl in self.bl_to_sv:
        #     print(f'{bl}: {self.bl_to_sv[bl]}')

        if len(self.logs)-old_logs == 1:
            self.logs.append(f'No anomalies detected.\n----------------------\n')
        else:
            self.logs.append('----------------------\n\n')

    def check_sequence_complete(self):
        for chrom in self.sequence:
            for i in range(1, len(self.sequence[chrom])):
                if self.sequence[chrom][i].start - self.sequence[chrom][i-1].end > 1:
                    self.logs.append(f'Missing bps in chromosome {chrom}: {self.sequence[chrom][i-1].end} to {self.sequence[chrom][i].start}.')
                elif self.sequence[chrom][i].start - self.sequence[chrom][i-1].end < 1:
                    self.logs.append(f'Duplicate bps sequenced in chromosome {chrom}: {self.sequence[chrom][i].start} to {self.sequence[chrom][i-1].end}.')


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
CNSegment = namedtuple('CNSegment', ['start', 'end', 'total_cn', 'major_cn', 'minor_cn', 'star'])

chromothripsis = '72f0a49a-aec8-47e5-846a-956c4da1507c.pcawg_consensus_1.6.161116.somatic.sv'
simple = 'e1217ebe-1826-41a9-b6c4-702100a66f5e.pcawg_consensus_1.6.161116.somatic.sv'
medium = '0ae2193f-0d68-485a-b8c2-7568cbcce33e.pcawg_consensus_1.6.161116.somatic.sv'
files = [simple+'.bedpe']
# file += '.bedpe'
# path = '/Users/siddharthsheth/Dropbox/work/research-projects/TQFT Cancer Progression/tcga/open/'
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path+'/..')
sv_path = 'data/tcga/open/'
cn_path = 'data/consensus.20170119.somatic.cna.tcga.public/'
sv_extension = '.pcawg_consensus_1.6.161116.somatic.sv.bedpe'
cn_extension = '.consensus.20170119.somatic.cna.txt'
log_file = 'svlog.log'
logs = []
patient_ids = [file.split('.')[0] for file in files]
# patient_ids = [file.split('.')[0] for file in os.listdir(sv_path)]

for i, id in enumerate(patient_ids):
    print(f'Working on patient number {i+1}: {id}.')
    patient = Patient(id)
    patient.check_valid() 

with open(log_file, 'w') as log:
    log.writelines(logs)