# Generating Proper AMGs
```python{cmd, hide, id=setup}
from aberration_multigraph.amg import AberrationMultigraph
from aberration_multigraph.generator import AMGGenerator
from matplotlib import pyplot as plt
```

An AMG is _proper_ if it is connected and none of the rejoin edges are parallel to any DSB edges.
Here we show this package can be used to analyze proper AMGs for a given DSB distribution.

## Summary Reports
This example shows how create an `AMGGenerator` object for a DSB distribution over 2 chromosomes with 2 DSBs each.
The summary report shows that of the 56 possible proper AMGs, 48 have cycle structure $C_4$ while 8 have cycle structure $C_2+C_2$.
Similar frequency distributions over diameter and girth are also available.

```python{cmd, continue=setup}
amg_gen = AMGGenerator(2, (2,2))
amg_gen.summarize()
```

## Accessing All Proper AMGs

This code chunk shows how to draw all possible proper AMGs over 3 chromosomes with 1 DSB each.
There are just 8 AMGs in this case.

```python{cmd, matplotlib, continue=setup}
amg_gen = AMGGenerator(3, (1,1,1), 'ABCDEFGHIJKL')
i = 1
for amg in amg_gen.generate_amgs():
    plt.subplot(2,4,i)
    amg.draw()
    i += 1
plt.show()
```