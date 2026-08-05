import heal_harness as H, json
p=H.p
d=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()
gates={}
for line in open('atoms/gates.jsonl'):
    g=json.loads(line); gates[g['t']]=(g['rhs'],tuple(g['vids']))
def show_tree(v,depth=0,seen=None):
    if seen is None: seen=set()
    pad="  "*depth
    if v in H.freeinp:
        print(f"{pad}x_{v} = FREE  (val={H.val[v]}, val%p={H.val[v]%p})")
        return
    if v in seen:
        print(f"{pad}x_{v} = (seen)"); return
    seen.add(v)
    rhs,vids=gates.get(v,('?',()))
    print(f"{pad}x_{v} = {rhs}   (val%p={H.val[v]%p})")
    if depth<6:
        for u in vids: show_tree(u,depth+1,seen)
print("========== x_2099 (G1) ==========")
show_tree(2099)
print("\nfree ancestors of x_2099:",sorted(H.anc[2099]&H.freeinp))
print("\n========== x_19964 (G2) ==========")
show_tree(19964)
print("\nfree ancestors of x_19964:",sorted(H.anc[19964]&H.freeinp))
