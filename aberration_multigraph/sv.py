from collections import defaultdict, namedtuple
import os

class Patient:
    def __init__(self, id):
        self.id = id
        self.logs = []
        self.sequence = self.get_sequence()
        self.svs = self.get_svs()
        self.bp_slack = 5

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
        for chrom in range(1,25):
            if chrom not in sequence:
                self.logs.append(f'Chromosome {chrom} not sequenced.\n')
            # elif len(sequence[chrom]) > 20:
            #     print(f'Patient {self.id}: chromosome {chrom} has {len(sequence[chrom])} segments.')
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
        self.logs.append(f'Patient {self.id}\n')

        self.check_sequence_complete()
        self.check_rejoins_cn()

        # for bl in self.bl_to_sv:
        #     print(f'{bl}: {self.bl_to_sv[bl]}')

        self.logs.append('----------------------\n\n')

    def check_sequence_complete(self):
        old_logs = len(self.logs)
        for chrom in self.sequence:
            for i in range(1, len(self.sequence[chrom])):
                if self.sequence[chrom][i].start - self.sequence[chrom][i-1].end > 1:
                    self.logs.append(f'Patient {self.id}: Missing bps in chromosome {chrom}: {self.sequence[chrom][i-1].end} to {self.sequence[chrom][i].start}.\n')
                elif self.sequence[chrom][i].start - self.sequence[chrom][i-1].end < 1:
                    self.logs.append(f'Patient {self.id}: Duplicate bps sequenced in chromosome {chrom}: {self.sequence[chrom][i].start} to {self.sequence[chrom][i-1].end}.\n')
        if len(self.logs) == old_logs:
            self.logs.append(f'Patient {self.id}: Sequencing data is complete.\n')

    def bp_loc_cn(self, bl):
        chrom = bl.chrom
        bp = bl.bp
        for segment in self.sequence[chrom]:
            if segment.start <= bp <= segment.end:
                return (segment.total_cn, 0) if segment.total_cn != 'NA' else (0, 'NA')
        
        # log = f'Patient {self.id}: Copy number not found for break location {bl}.'
        if bp < self.sequence[chrom][0].start:
            return (0, '<')
            # log += ' Location less than sequenced range.'
        elif bp > self.sequence[chrom][-1].end:
            return (0, '>')
            # log += ' Location greater than sequenced range.'
        # self.logs.append(log+'\n')
        # return (0, 'Not sequenced')
    
    def bp_loc_rejoins(self, bl):
        rejoins = 0
        for sv, pos in self.bl_to_sv[bl]['svs']:
            if pos == 1 and sv.strand1 == '-':
                rejoins +=1
            elif pos == 2 and sv.strand1 == '+':
                rejoins += 1
            elif pos == 3 and sv.strand2 == '-':
                rejoins += 1
            elif pos == 4 and sv.strand2 == '+':
                rejoins += 1
        return rejoins

    def check_rejoins_cn(self):
        old_logs = len(self.logs)
        bp_locs = sorted([bl for bl in self.bl_to_sv])
        for bp_loc in bp_locs:
            rejoins = self.bp_loc_rejoins(bp_loc)
            copies, error = self.bp_loc_cn(bp_loc)
            if rejoins > copies:
                bp_less_slack = BreakLocation(bp_loc.chrom, bp_loc.bp-self.bp_slack)
                bp_less_slack_cn = self.bp_loc_cn(bp_less_slack)[0]
                bp_add_slack = BreakLocation(bp_loc.chrom, bp_loc.bp+self.bp_slack)
                bp_add_slack_cn = self.bp_loc_cn(bp_add_slack)[0]
                if error == 0:
                    if bp_less_slack_cn <= copies:
                        pass
                        # self.logs.append(f"Patient {self.id}: {bp_loc} is rejoined in {rejoins} SVs while CN is {copies}. But a location -{self.bp_slack} is fine.\n")
                    elif  bp_add_slack_cn <= copies:
                        pass
                        # self.logs.append(f"Patient {self.id}: {bp_loc} is rejoined in {rejoins} SVs while CN is {copies}. But a location +{self.bp_slack} is fine.\n")
                    else:
                        self.logs.append(f"Patient {self.id}: {bp_loc} is rejoined in {rejoins} SVs while CN is {copies}.\n")
                elif error == '<':
                    if bp_add_slack_cn <= copies:
                        pass
                        # self.logs.append(f'Patient {self.id}: {bp_loc} CN not found for break location. Location less than sequenced range. But a location +{self.bp_slack} is fine.\n')
                    else:
                        self.logs.append(f'Patient {self.id}: {bp_loc} CN not found for break location. Location less than sequenced range.\n')
                elif error == '>':
                    if bp_less_slack_cn  <= copies:
                        pass
                        # self.logs.append(f'Patient {self.id}: {bp_loc} CN not found for break location. Location greater than sequenced range. But a location -{self.bp_slack} is fine.\n')
                    else:
                        self.logs.append(f'Patient {self.id}: {bp_loc} CN not found for break location. Location greater than sequenced range.\n')
                elif error == 'NA':
                    self.logs.append(f'Patient {self.id}: {bp_loc} CN not found for break location. Marked NA.\n')
                else:
                    self.logs.append(f'Patient {self.id}: {bp_loc} CN not found for break location.\n')
        if len(self.logs) == old_logs:
            self.logs.append(f'Patient {self.id}: Success! For every breakpoint #(rejoins) <= #(copies).\n')

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
    
    def __eq__(self, other):
        return True if self.chrom == other.chrom and self.bp == other.bp else False

# class StructuralVariation:
#     def __init__(self, c1, s1, e1, c2, s2, e2, sv_id, pe_support, str1, str2, svclass, svmethod):
#         c1 = 23 if c1 == 'X' else c1
#         c1 = 24 if c1 == 'Y' else c1
#         c2 = 23 if c2 == 'X' else c2
#         c2 = 24 if c2 == 'Y' else c2
#         self.chrom1= int(c1)
#         self.start1 = int(s1)
#         self.end1 = int(e1)
#         self.chrom2 = int(c2)
#         self.start2 = int(s2)
#         self.end2 = int(e2)
#         self.sv_id = sv_id
#         self.pe_support = pe_support
#         self.strand1 = str1
#         self.strand2 = str2
#         self.svclass = svclass
#         self.svmethod = svmethod

#     def __str__(self):
#         return str((self.chrom1, self.start1, self.end1, self.chrom2, self.start2, self.end2))

BreakLocation = namedtuple('BreakLocation', ['chrom', 'bp'])
CNSegment = namedtuple('CNSegment', ['start', 'end', 'total_cn', 'major_cn', 'minor_cn', 'star'])
StructuralVariation = namedtuple('StructuralVariation', ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'sv_id', 'pe_support', 'strand1', 'strand2', 'svclass', 'svmethod'])

chromothripsis = '72f0a49a-aec8-47e5-846a-956c4da1507c.pcawg_consensus_1.6.161116.somatic.sv.bedpe'
simple = 'e1217ebe-1826-41a9-b6c4-702100a66f5e.pcawg_consensus_1.6.161116.somatic.sv.bedpe'
medium = '0ae2193f-0d68-485a-b8c2-7568cbcce33e.pcawg_consensus_1.6.161116.somatic.sv.bedpe'
files = [simple, medium]
# file += '.bedpe'
# path = '/Users/siddharthsheth/Dropbox/work/research-projects/TQFT Cancer Progression/tcga/open/'
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path+'/..')
sv_path = 'data/tcga/open/'
cn_path = 'data/consensus.20170119.somatic.cna.tcga.public/'
sv_extension = '.pcawg_consensus_1.6.161116.somatic.sv.bedpe'
cn_extension = '.consensus.20170119.somatic.cna.txt'
log_file = 'svlog.log'
if os.path.exists(log_file):
    os.remove(log_file)
logs = []
# patient_ids = [file.split('.')[0] for file in files]
patient_ids = [file.split('.')[0] for file in os.listdir(sv_path)]
for i, id in enumerate(patient_ids):
    print(f'Working on patient number {i+1}: {id}.')
    patient = Patient(id)
    patient.check_valid() 
    with open(log_file, 'a') as log:
        log.writelines(patient.logs)