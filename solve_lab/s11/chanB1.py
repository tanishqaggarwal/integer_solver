"""What message does the 39,018 state use, and does it sit at a weight-2 optimum?"""
import sys, os, json, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from ip7 import load_raw
P=L.P; sys.set_int_max_str_digits(400000)
real=[r[1] for r in json.load(open(os.path.join(HERE,'data','gmp16.json')))]
tr=json.load(open(os.path.join(HERE,'data','bits_trees.json')))
tree={}
for k,v in tr.items():
    for b in v: tree[b]=k
for nm,f in [('39018 state', os.path.join(HERE,'data','finish3_named.json')),
             ('checkpoint', os.path.join(HERE,'..','best','new_instance_partial_39026.json'))]:
    v=[x%P for x in load_raw(f)]
    on=[b for b in real if v[b]%P==1]
    print(f"{nm}: message weight {len(on)} -> {[(b,tree.get(b)) for b in on]}")
    forwardp(v)
    print(f"   U={v[7715]} V={v[34554]}  channels {v[15298]}/{v[5647]}/{v[34606]}")
