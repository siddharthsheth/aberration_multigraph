```python{cmd=true, id=setup, hide=true}
from aberration_multigraph.amg import AberrationMultigraph
from matplotlib import pyplot as plt
```

# Transformations on AMGs

Here, we illustrate the dynamics of AMGs using the following as an example.

```python{cmd=true, continue=setup, id=main}
chromatin = [(1,2),(3,4),(5,6),(7,8),(9,10),(11,12),(13,14),(15,16)]
dsb = [(2,3),(6,7),(10,11),(14,15)]
rejoin = [(2,6),(3,11),(7,14),(10,15)]
amg = AberrationMultigraph(chromatin, dsb, rejoin)
```

## Total Swaps

A total swap on an AMG takes two chromosomes and exchanges their location in the AMG.
The vertices are renamed but the edges are preserved from the original AMG.


```python{cmd=true, continue=main, matplotlib=true}
plt.figure(figsize=(6,6))
plt.subplot(2,2,1)
plt.title('amg')
amg.draw()

plt.subplot(2,2,2)
plt.title('amg.total_swap(0,1)')
amg.total_swap(0,1).draw()

plt.subplot(2,2,3)
plt.title('amg.total_swap(0,2)')
amg.total_swap(0,2).draw()

plt.subplot(2,2,4)
plt.title('amg.total_swap(0,3)')
amg.total_swap(0,3).draw()

plt.show()
```
## Total Twists

A total twist on an AMG takes a chromosome and reverses it.
The chromosome is reindexed to match the old ordering.
The edges are updated so that old neighbors are preserved in the new ordering.

```python{cmd=true, continue=main, matplotlib=true}
plt.figure(figsize=(6,6))
plt.subplot(2,2,1)
plt.title('amg')
amg.draw()
plt.subplot(2,2,2)
plt.title('amg.total_twist(0)')
amg.total_twist(0).draw()
plt.subplot(2,2,3)
plt.title('amg.total_twist(1)')
amg.total_twist(1).draw()
plt.subplot(2,2,4)
plt.title('amg.total_twist(2)')
amg.total_twist(2).draw()
plt.show()
```