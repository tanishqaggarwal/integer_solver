import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from bits5 import msg, fails
P=L.P; sys.set_int_max_str_digits(400000)
for S in [{24601,2081},{24601,4287},{24601,13195},{24601,12054},{24601,16586},{24601,24365}]:
    v=msg(S); F=fails(v)
    r=[evalp(L.polys[a],v) for a in F]
    print(f"{sorted(S)}: {len(F)} failing {F}")
    print(f"     residues differ from checkpoint: ", end='')
    print([str(x)[:12]+'..' for x in r])
