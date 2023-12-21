from collections import defaultdict, namedtuple
import os

BreakLocation = namedtuple('BreakLocation', ['chrom', 'bp'])
CNSegment = namedtuple('CNSegment', ['start', 
                                     'end', 
                                     'total_cn', 
                                     'major_cn', 
                                     'minor_cn', 
                                     'star'
                                    ])

class SVVertex:
    def __init__(self, chrom, bp, dsb=None, rejoin=None, chromatin=None):
        self.chrom = chrom
        self.bp = bp
        self.dsb = dsb
        self.rejoin = rejoin
        self.chromatin = chromatin

    def amg_vertex(self):
        return (self.chrom, self.bp)
    
    def __str__(self) -> str:
        return f'({self.chrom}, {self.bp})'
    
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

StructuralVariation = namedtuple('StructuralVariation', ['chrom1', 
                                                         'start1',
                                                         'end1', 
                                                         'chrom2', 
                                                         'start2', 
                                                         'end2', 
                                                         'sv_id', 
                                                         'pe_support', 
                                                         'strand1', 
                                                         'strand2', 
                                                         'svclass', 
                                                         'svmethod'
                                                        ])

