#!/usr/bin/env python3
"""K8: map the congruence-pin chain around the obstruction."""
import sys, os, json, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
from cascade import Cascade, NV, P
from circ2 import vars_of

C = Cascade()
for w in [9118, 8731, 4432, 7068, 2099, 19964, 1329, 10903, 17325, 9413, 7075, 21279]:
    print('=== x%d appears in %d atoms' % (w, len(C.var2atoms[w])))
    for i in C.var2atoms[w]:
        print('    ', C.atomnames[i][:110])
