```python{cmd=true, id=setup, hide=true}
from aberration_multigraph.amg import AberrationMultigraph
from matplotlib import pyplot as plt
```

# Transformations on AMGs

Here, we illustrate the dynamics of AMGs.

## Total Swaps

## Total Twists

A total twist 
```python{cmd=true, continue=setup, matplotlib=true}
chromatin = [(1,2),(3,4),(5,6),(7,8),(9,10),(11,12),(13,14),(15,16)]
dsb = [(2,3),(6,7),(10,11),(14,15)]
rejoin = [(2,6),(3,11),(7,14),(10,15)]

amg_1 = AberrationMultigraph(chromatin, dsb, rejoin)
assert amg_1.num_chromosome == 4

amg_2 = amg_1.total_twist(1)
amg_3 = amg_1.total_twist(2)
amg_4 = amg_2.total_twist(3)

plt.figure(figsize=(4,8))
plt.subplot(4,1,1)
plt.title('amg')
amg_1.draw()
plt.subplot(4,1,2)
plt.title('amg.total_twist(1)')
amg_2.draw()
plt.subplot(4,1,3)
plt.title('amg.total_twist(2)')
amg_3.draw()
plt.subplot(4,1,4)
plt.title('amg.total_twist(3)')
amg_4.draw()
plt.show()
```