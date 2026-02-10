# Completing Incomplete AMGs

```python{cmd=true, id=setup, hide}
from aberration_multigraph.incomplete_amg import IncompleteAMG
```

An incomplete aberration multigraph (`IncompleteAMG`) represents a partial rejoining of double-strand break (DSB) ends. Such objects arise naturally during recursive enumeration, where only some rejoin edges have been fixed and the remaining free ends must still be paired.

This example demonstrates how to use `IncompleteAMG` to count and enumerate all complete aberration multigraphs consistent with a given partial rejoining.

## Defining an Incomplete AMG

We specify:
	•	the chromatin backbone,
	•	the DSB edges,
	•	and a partial set of rejoin edges.

```python{cmd=true, id=main, continue=setup}
chromatin = [
    (1,2), (3,4), (5,6), (7,8), (9,10), (11,12), (13,14),
    (15,16), (17,18), (19,20), (21,22), (23,24), (25,26)
]

dsb = [
    (2,3), (4,5), (8,9), (10,11), (12,13),
    (14,15), (18,19), (20,21), (22,23), (24,25)
]

rejoin = [(3,5), (9,10), (13,14), (19,20), (23,24)]

inc = IncompleteAMG(chromatin, dsb, rejoin)
```

At this stage, the rejoin edges do not form a perfect matching, so the object represents many possible complete AMGs.

## Counting Completions

The total number of complete aberration multigraphs extending this partial configuration can be computed directly:

```python{cmd=true, id=count, continue=main}
print(inc.count_amgs())
```

This uses an internal recursive backtracking procedure that pairs the remaining free vertices in all valid ways.

## Enumerating All Completions

To generate each complete AMG explicitly, use complete_amgs():

```python{cmd=true, id=print, continue=main}
for amg in inc.complete_amgs():
    print(amg.rejoin_edges)
```

Each yielded object is a fully specified AberrationMultigraph consistent with the initial partial rejoining.

## Remarks
- Incomplete AMGs are intermediate combinatorial objects and do not enforce connectivity.
- Enumeration is combinatorial in nature and intended for small to moderate instances.
- This interface is primarily used internally by generators, but can also be useful for targeted counting problems.