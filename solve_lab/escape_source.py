#!/usr/bin/env python3
"""NEW activation: at x_12779=2 (x_14402=-1) the div wire a1813 computes
x_24026 = 321447*x_38215/x_14402 = -321447*x_38215 CORRECTLY (no div-by-zero, unlike
x_12779=1). The only blocker is x_38215=x_37917*x_30077==0. So FREEZE the escape
source x_37917/x_30077 (and x_7815/x_31807 for the 3183 side) to make the twist hold,
with x_12779=2. Then a1813/a1815 stay SATISFIED (activation is consistent, not the
frozen-x_24026 hack). Count violations -- hope: fewer broken atoms than slack_active."""
import json, time
from confluent_eval5 import build5, make_forward
from slack_active import viol_atoms
from propagate import NVARS

def make_frz(kind, info, seq, bestval, FREEZE):
    seq2=[v for v in seq if v not in FREEZE]
    def run(val, frozen):
        for v in FREEZE: val[v]=frozen.get(v,val[v])
        for v in seq2:
            k=kind[v]
            if k=='gate':
                coef,terms=info[v]; rs=0
                for c,m in terms:
                    t=c
                    for x in m: t*=val[x]
                    rs+=t
                if coef and (-rs)%coef==0: val[v]=(-rs)//coef
            elif k=='load':
                bit,cbx,lt=info[v]
                if val[bit]==0: val[v]=0
                else:
                    rest=0
                    for c,m in lt:
                        t=c
                        for x in m: t*=(1 if x==bit else val[x])
                        rest+=t
                    num=-rest; den=cbx*val[bit]
                    if den and num%den==0: val[v]=num//den
            elif k=='div':
                c,u,rest=info[v]; rs=0
                for cc,m in rest:
                    t=cc
                    for x in m: t*=val[x]
                    rs+=t
                den=c*val[u]
                if den and (-rs)%den==0: val[v]=(-rs)//den
                elif den==0: val[v]=0
        return val
    return run

def main():
    t0=time.time()
    A,kind,info,seq0,bestval,ncyc=build5()
    order=json.load(open('eval_order.json'))['order']
    defset=set(v for v in kind if kind[v]!='const')
    seq=[v for v in order if v in defset and v not in (9770,3183)]
    seq+=[v for v in (9770,3183) if v in defset]
    seq+=[v for v in defset if v not in set(order) and v not in (9770,3183)]
    solve=make_forward(kind,info,seq,bestval)
    control=json.load(open('control_bits.json'))
    print(f"kinds: x_37917={kind.get(37917)}, x_30077={kind.get(30077)}, x_7815={kind.get(7815)}, x_31807={kind.get(31807)}, x_24026={kind.get(24026)}, x_14402={kind.get(14402)}", flush=True)

    # find a bit-setting giving x_12779=2
    st=5
    def rnd():
        nonlocal st; st=(st*6364136223846793005+1442695040888963407)&((1<<64)-1); return st>>33
    setbits=None
    for _ in range(600):
        k=8+rnd()%30
        S=sorted(set(control[rnd()%len(control)] for _ in range(k)))
        if solve(list(bestval),S)[12779]==2: setbits=S; break
    if setbits is None:
        print("no x_12779=2 setting found", flush=True); return
    v1=solve(list(bestval),setbits)
    print(f"x_12779=2 bits (k={len(setbits)}); x_14402={v1[14402]}, x_35186={v1[35186]!=0}, x_18274={v1[18274]!=0}", flush=True)

    # target: x_38215 = (x_35186 - x_18274)/642894 so x_9770=x_18274; x_29437-analog for 3183
    d9=v1[35186]-v1[18274]
    d3=v1[1642]-v1[17728]
    prod38=d9//642894 if d9%642894==0 else None
    prod29=d3//2 if d3%2==0 else None
    print(f"needed x_38215=(x_35186-x_18274)/642894 integer? {d9%642894==0}; x_29437=(x_1642-x_17728)/2 integer? {d3%2==0}", flush=True)
    # a1815: x_14402*x_27116 = x_29437 ; x_27116=x_29437/x_14402 = -x_29437 (x_14402=-1); x_3183=x_1642+2*x_27116
    # need x_3183=x_17728 => x_1642+2*(-x_29437) = x_17728 => x_29437=(x_1642-x_17728)/2
    run=make_frz(kind,info,seq,bestval,{37917,30077,7815,31807})
    if prod38 is not None and prod29 is not None:
        frozen={37917:prod38, 30077:1, 7815:prod29, 31807:1}
        v2=run(list(v1),frozen)
        tw=(v2[9770]==v2[18274], v2[3183]==v2[17728])
        print(f"  x_38215={v2[38215]}, x_24026={v2[24026]!=0}, x_3368={v2[3368]!=0}", flush=True)
        print(f"  twist: x9770==x18274 {tw[0]}, x3183==x17728 {tw[1]}", flush=True)
        bad=viol_atoms(A,v2)
        print(f"  violated atoms: {len(bad)}: {sorted(bad)[:16]}", flush=True)
    else:
        print("  target not integer at this bit-setting; try other x_12779=2 settings", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__=='__main__':
    main()
