# aberration_multigraph

This project aims to study and analyze aberration multigraphs[^1], represented as `AberrationMultigraph`, as a model for chromosome aberrations.
The data source is the PCAWG database[^3].
Each file in the database is a collection of structural variations [^2].
Each entry is represented as a `StructuralVariation` object.
Both, `AberrationMultigraph` and `StructuralVariation`, are different kinds of graphs.

## Aberration Multigraph
1. Vertices are of 2 types: telomere and non-telomere.
2. The degree of a telomere vertex is exactly 1, while that of a non-telomere vertex is exactly 3.
3. Edges are of 3 types: chromatin, DSB, and rejoin.
4. A telomere vertex is adjacent to a chromatin edge, while a non-telomere vertex is adjacent to exactly 1 edge of each type.

## Structural Variations
1. A structural variation comprises 2 DSB edges and 1 rejoin edge. <!-- Each edge has exactly 2 vertices. -->
1 vertex in every DSB has degree 2 and the other has degree 1.
2. Every SV need not correspond to a copy number variation.
Reciprocal SVs do not cause CNV. <!-- However, those SVs that do correspond to CNVs are important. -->
3. A DSB edge can appear in at most 2 SVs.
The same vertex of the DSB cannot participate in both SVs.
4. There cannot be 2 DSB edges over successive base pairs, e.g., if $(u,v)$ and $(w,y)$ are consecutive DSB edges in the sorted list of DSB edges for a particular chromosome, then $w>v+1$.


[^1]: [[https://pubmed.ncbi.nlm.nih.gov/12385633/| Using graph theory to describe and model chromosome aberrations]]
[^2]: [[https://www.nature.com/articles/s41586-019-1913-9 | Patterns of somatic structural variation in human cancer genomes]]
[^3]: [[https://dcc.icgc.org/releases/PCAWG/consensus_sv | PCAWG Consensus Callsets for Structural Variants]]