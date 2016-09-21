chemistry.py
============

Forked from [https://github.com/PythonNut/chemistry.py](PythonNut/chemistry.py)
A module that adds symbolic chemistry to Python/SAGEmath. Major WIP.

```
Python 3.4.1 (default, May 19 2014, 17:23:49) 
Type "copyright", "credits" or "license" for more information.

IPython 2.1.0 -- An enhanced Interactive Python.
?         -> Introduction and overview of IPython's features.
%quickref -> Quick reference.
help      -> Python's own help system.
object?   -> Details about 'object', use 'object??' for extra details.

In [1]: from chemistry import *

In [2]: Ag
Out[2]: Ag

In [3]: Ag.mass, Ag.melting_point, Ag.electron_affinity
Out[3]: (107.8682, 1235.1, 1.30447)

In [4]: Ag.group, Ag.period, Ag.block, Ag.family
Out[4]: ('1', 5, 'd', 'Transition')

In [5]: Ag.electron_config_full
Out[5]: '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s1'

In [6]: f3s = Fe*2 + (S + O*4)*3 # Iron (III) Sulfate

In [7]: f3s.mass
Out[7]: 399.858

In [8]: f3s.masses()
Out[8]: {Fe: 55.845, O: 15.999, S: 32.06}

In [9]: f3s.mass_composition
Out[9]: {Fe: 27.932416007682725, O: 48.01404498597002, S: 24.053539006347254}

In [10]: f3s.counts
Out[10]: {Fe: 2, O: 12, S: 3}
```
