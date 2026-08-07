#!/usr/bin/env python3
"""Q-7 VALIDATION of the whole reduction: take the deliverable assignment, read its leaf
selector bits, form the GROUP SUM of the selected leaves, and compare with the value the
circuit actually produces at the root-target wire.  If they agree the tree->group collapse is real."""
import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
sys.path.insert(0, os.path.join(HERE,'..','agentF_work'))
from qgrp import *
from fwd import Engine, NV
E=Engine()
asg=json.load(open(os.path.join(HERE,'..','best','new_instance_partial_39026.json')))
v=[0]*NV
for k,x in asg.items():
    i=int(k.split('_')[1]) if '_' in k else int(k)
    if i<NV: v[i]=int(x)
E.run(v)
pins=json.load(open(os.path.join(HERE,'..','agentF_work','pins.json')))
on=[]; Lv=leaves()
for g,vv in pins.items():
    if v[int(g)]%p==1: on.append(int(g))
print('leaf selectors ON:', len(on), sorted(on)[:20])
tot=None
for g in on:
    if g in Lv: tot=add(tot,Lv[g])
    else: print('  ON leaf %d has no decoded point'%g)
print('group sum of ON leaves =', tot)
print('circuit x24468 mod p   =', v[24468]%p, ' (target x C1 =', 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002%p,')')
print('circuit x18956 mod p   =', v[18956]%p)
print('circuit x13682 mod p   =', v[13682]%p, ' x37892 mod p =', v[37892]%p)
if tot: print('shifted->raw of sum   = x:',(tot[0]-cs)%p,' y:',tot[1])
