# Aberration Multigraphs
```python{cmd=true, id=setup, hide}
from aberration_multigraph.amg import AberrationMultigraph
from matplotlib import pyplot as plt
```

Cells experience breakages at a DNA-level when exposed to certain kinds of radiation.
Such a breakage of both DNA strands in a chromosome is called a _double-strand break_ (DSB).
A collection of chromosomes undergoing DSBs can rejoin, sometimes incorrectly, usually producing a different arrangement than the original one, referred to as the final configuration.
A _chromosome aberration_ is the process consisting of such misrepairs/misrejoinings.

A chromosome aberration is defined by an _initial configuration_ and a _final configuration_.
Both of these configurations can be represented by edge-colored graphs; edges in the initial cofiguration are colored _chromatin_ and _DSB_, while edges in the final configuration are colored _chromatin_ and _rejoin_.
Given graphs corresponding to an initial and final configuration, an aberration can be modeled using aberration multigraphs (AMGs).

## Using `AberrationMultigraph` to Model AMGs

The following example shows how to use `AberrationMultigraph` to represent two AMGs on the same initial configuration, but different final configurations.

The chromatin edges are colored black, DSB edges red, and rejoin edges are green.
```python{cmd=true, matplotlib=true, continue=setup}
chrom = (('A','B'), ('C','D'), ('E','F'), ('G','H'), ('I','J'), ('K','L'))
dsb = (('B','C'), ('D','E'), ('H','I'), ('J','K'))

rejoin_1 = (('B', 'H'), ('D', 'J'), ('C', 'I'), ('E', 'K'))
amg_1 = AberrationMultigraph(chrom, dsb, rejoin_1)

rejoin_2 = (('B','J'), ('C','K'), ('D','I'), ('E','H'))
amg_2 = AberrationMultigraph(chrom, dsb, rejoin_2)

plt.figure(figsize=(7.5,2))
plt.subplot(1,2,1)
amg_1.draw()
plt.subplot(1,2,2)
amg_2.draw()

plt.show()
```

## Analyzing AMGs with `AberrationMultigraph`
The following example uses `AberrationMultigraph` to model an AMG on 1 chromosome with 3 double-strand breaks.
The rejoin edges result in an AMG with diameter 4, cycle structure $C_3$ and girth 3.

```python{cmd=true, matplotlib=true, continue=setup}
amg = AberrationMultigraph((('A','B'), ('C','D'), ('E','F'), ('G','H')),
                            (('B','C'), ('D','E'), ('F','G')),
                            (('B','E'), ('C','G'), ('D','F')))

print(f"""
        Number of chromosomes: {amg.num_chromosome}
        Diameter: {amg.diameter()}
        Cycle structure: {amg.cycle_structure()}
        """)

plt.figure(figsize=(6,1))
plt.margins(y=1)
amg.draw()
plt.show()

```

