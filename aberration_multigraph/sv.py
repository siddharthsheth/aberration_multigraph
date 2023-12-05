from collections import defaultdict

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

chromothripsis = '72f0a49a-aec8-47e5-846a-956c4da1507c.pcawg_consensus_1.6.161116.somatic.sv'
simple = 'e1217ebe-1826-41a9-b6c4-702100a66f5e.pcawg_consensus_1.6.161116.somatic.sv'
medium = '0ae2193f-0d68-485a-b8c2-7568cbcce33e.pcawg_consensus_1.6.161116.somatic.sv'
file = medium
file += '.txt'
path = '/Users/siddharthsheth/Dropbox/work/research-projects/TQFT Cancer Progression/tcga/open/'
log_file = '/Users/siddharthsheth/Dropbox/work/code-projects/aberration_multigraph/svlog.log'
logs = []
with open(file, 'r') as svfile:
    headers = svfile.readline()
    svs = []
    vertices = dict()
    for entry in svfile:
        entry = entry.split()
        sv = StructuralVariation(*entry)
        svs.append(sv)
        u_1, v_1, u_2, v_2 = SVVertex(int(entry[0]), int(entry[1])), SVVertex(int(entry[0]), int(entry[2])), SVVertex(int(entry[3]), int(entry[4])), SVVertex(int(entry[3]), int(entry[5]))
        if u_1 in vertices:
            # if vertices[u_1].strand1 == '-':
            if vertices[u_1].strand1 == entry[8]:
                logs.append(f'Duplicate rejoining \t\t\t vertex {str(u_1)} \t\t\t SV {str(sv)}')
            else:
                logs.append(f'Potential reciprocal \t\t\t vertex {str(u_1)} \t\t\t SV {str(sv)}')
        else:
            vertices[u_1] = sv
        if v_1 in vertices:
            if vertices[v_1].strand1 == entry[9]:
                logs.append(f'Duplicate rejoining \t\t\t vertex {str(u_1)} \t\t\t SV {str(sv)}')
            else:
                logs.append(f'Potential reciprocal \t\t\t vertex {str(u_1)} \t\t\t SV {str(sv)}')
        else:
            vertices[v_1] = sv
        if u_2 in vertices:
            logs.append(f'Duplicate vertex: {str(u_2)}')
        else:
            vertices[u_2] = sv
        if v_1 in vertices:
            logs.append(f'Duplicate vertex: {str(v_1)}')
        else:
            vertices[v_2] = sv
        
        