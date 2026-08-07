#!/usr/bin/env python3
"""Write the engine input files: line 1 = base target, then the 256 ladder points."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 'ydata.json')))
lad = [(int(a), int(b)) for a, b in d['ladder']]

def w(path, P):
    with open(path, 'w') as f:
        f.write('%d %d\n' % (P[0], P[1]))
        for x, y in lad: f.write('%d %d\n' % (x, y))
    print('wrote %s  base=(%d, %d)' % (path, P[0], P[1]))

w(os.path.join(HERE, 'data_comp.txt'), (int(d['Tp'][0]), int(d['Tp'][1])))   # T'
w(os.path.join(HERE, 'data_fwd.txt'),  (int(d['T'][0]),  int(d['T'][1])))    # T  (control)
