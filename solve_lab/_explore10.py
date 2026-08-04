import heal_harness as H, json, pickle
p=H.p
def trace(v, depth=0, maxd=4, seen=None):
    if seen is None: seen=set()
    ind="  "*depth
    if v in H.freeinp:
        print(f"{ind}x_{v} = FREE")
        return
    if v in seen or depth>maxd:
        print(f"{ind}x_{v} = ... (gate, #anc={len(H.anc[v])})")
        return
    seen.add(v)
    gi=H.definer[v]; t,rhs,vids=H.gates[gi]
    print(f"{ind}x_{v} = {rhs}   [free_anc={len(H.anc[v])}]")
    for u in vids:
        trace(u, depth+1, maxd, seen)
for v in [9254,29967,1308,19083]:
    print(f"\n########## TRACE x_{v} ##########")
    trace(v, maxd=3)
    print(f"  free ancestors of x_{v}: {sorted(H.anc[v])}")
