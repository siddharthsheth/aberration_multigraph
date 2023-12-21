from collections import defaultdict, namedtuple
import os

BreakLocation = namedtuple('BreakLocation', ['chrom', 'bp'])
BreakLocation.__doc__ = """
A breakpoint location in a chromosome.

Parameters
----------
chrom : int
    The chromosome in which the break occurred.

bp : int
    The base pair where the break occurred.
"""

CNSegment = namedtuple('CNSegment', ['start', 
                                     'end', 
                                     'total_cn', 
                                     'major_cn', 
                                     'minor_cn', 
                                     'star'
                                    ])
CNSegment.__doc__ = """
A contiguous segment in a chromosome whose copy number was sequenced.
Should be stored as a dictionary keyed by chromosome number.
The parameters follow from the PCAWG project dataset.

Parameters
----------
start : int
    The base pair where this copy number segment begins.
end : int
    The base pair where this copy number segment ends.
total_cn : int
    The total copy number of this segment.
major_cn : int
    The major copy number of this segment.
minor_cn : int
    The minor copy number of this segment.
star : int
    The confidence of the measurement.
"""

class SVVertex:
    """
    A class to represent vertices in structural variations.

    There should be only one SVVertex per BreakLocation.

    Attributes
    ----------
    chrom : int
        The chromosome this SVVertex belongs to.
    bp : int
        The base pair location of this SVVertex.
    dsb : SVVertex
        If this SVVertex is involved in a double-strand break edge.
    rejoin : SVVertex
        If this SVVertex is involved in a rejoin with another SVVertex.
    chromatin : SVVertex
        If this SVVertex is the chromatin neighbor of another SVVertex.
    """
    def __init__(self, chrom, bp, dsb=None, rejoin=None, chromatin=None):
        """_summary_

        Parameters
        ----------
        chrom : int
            The chromosome this SVVertex belongs to.
        bp : int
            The base pair location of this SVVertex.
        dsb : SVVertex, optional
            If this SVVertex is involved in a double-strand break edge,
             by default None
        rejoin : SVVertex, optional
            If this SVVertex is involved in a rejoin with another SVVertex,
             by default None
        chromatin : SVVertex, optional
            If this SVVertex is the chromatin neighbor of another SVVertex,
             by default None
        """
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

