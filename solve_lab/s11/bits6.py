import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from bits5 import msg, fails
P=L.P; sys.set_int_max_str_digits(400000)
which=int(sys.argv[1]); out=sys.argv[2]
v=msg({which})
print(f"message {{x{which}}}: U={v[7715]} V={v[34554]} channels {v[15298]}/{v[5647]}/{v[34606]}")
F=fails(v); print(f"  failing checks mod p: {len(F)} {F}")
json.dump([int(x) for x in v], open(out,'w')); print("saved",out)
