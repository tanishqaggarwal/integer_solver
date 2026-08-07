"""Is there a control for x_12186 that does NOT need a*b = 1 (which lights the group-2 mirror)?"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, uv10
P = L.P
FREE = [u for u in range(L.NVARS) if L.definer.get(u) is None]

for BITS in [(542,), (542, 47), (438,), (542, 438)]:
    v = uv10.state(BITS, {})
    print(f"=== BITS={BITS}  U={v[7715]} V={v[34554]} ab={v[38170]} cd={v[3896]} "
          f"x15298={v[15298]} x34606={v[34606]} x5647={v[5647]}")
    print(f"    x_10603={str(v[10603])[:40]}  (x_25758 = x_10603 * x_33612, x_33612 free)")
    print(f"    x_10603 % P == 0 : {v[10603] % P == 0}")
    # scan for mod-p controls of x_12186 and x_1308
    hits = {'x12186': [], 'x1308': [], 'x24908': [], 'x19083': []}
    base = (v[12186], v[1308], v[24908], v[19083])
    for u in FREE:
        old = v[u]
        v[u] = old + 1
        fw.forward(v)
        for k, idx in (('x12186', 12186), ('x1308', 1308), ('x24908', 24908), ('x19083', 19083)):
            if (v[idx] - base[('x12186', 'x1308', 'x24908', 'x19083').index(k)]) % P:
                hits[k].append(u)
        v[u] = old
        fw.forward(v)
    for k in hits:
        print(f"    {k}: {len(hits[k])} mod-p controls -> {hits[k][:12]}")
