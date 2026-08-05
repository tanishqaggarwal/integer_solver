"""Bit-pattern lab: apply a 256-selector pattern via pinrec, forward, count fails."""
import heal_harness as H, json
p=H.p
pin=json.load(open('pinrec.json'))
sels=sorted(set(r[1] for r in pin))
SIDX={s:i for i,s in enumerate(sels)}
from collections import defaultdict
bysel=defaultdict(list)
for r in pin: bysel[r[1]].append(r)
vA=H.loadd('best_agentA_39022.json')
BASE={v:vA.get(v,0) for v in H.freeinp}
AGENTA_BITS={s:vA.get(s,0) for s in sels}

def apply_pattern(bits, base=None, twopass=True):
    """bits: dict sel->0/1 (missing => 0). Returns list of failing eq indices."""
    b=base if base is not None else BASE
    for v in H.freeinp: H.val[v]=b.get(v,0)
    for s in sels: H.val[s]=int(bits.get(s,0))
    # first pass: on-target=CONST, off-target=0
    for (pa,se,tg,C,cf,hd) in pin:
        H.val[tg]= C if bits.get(se,0)==1 else 0
    H.forward()
    if twopass:
        # correct targets for nonzero handles
        for (pa,se,tg,C,cf,hd) in pin:
            if bits.get(se,0)==1:
                H.val[tg]=C+cf*H.val[hd]
            else:
                H.val[tg]=0
        H.forward()
    return H.fails()

def handles_nonzero():
    return [(pa,hd,H.val[hd]) for (pa,se,tg,C,cf,hd) in pin if H.val[hd]!=0]

if __name__=='__main__':
    F=apply_pattern(AGENTA_BITS)
    print('agentA pattern via apply:',len(F),'fails:',F)
    hz=handles_nonzero()
    print('nonzero handles:',len(hz))
    # sanity: all-off pattern
    F0=apply_pattern({s:0 for s in sels})
    print('all-off pattern:',len(F0),'fails')
