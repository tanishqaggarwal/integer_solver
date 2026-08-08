#!/usr/bin/env python3
"""demo.py -- end-to-end faithfulness proof of the encoding on scaled instances."""
import sys, time
from ecsmall import curve, find
from ladder import build

def instance(p, B, m, k):
    add, mul = curve(p, B)
    G, order = find(p, B)
    assert order >= (1 << m), f"order {order} too small for {m} bits"
    T = mul(k, G)
    pts = [mul(1 << i, G) for i in range(m)]
    Tpp = add(mul(2, T), mul((1 << m) - 1, mul(order - 1, G)))   # 2T - (2^m-1)G
    return G, order, T, pts, Tpp, add, mul

def run(p, B, m, k, chunk=16, mode='binary', exhaustive=True):
    G, order, T, pts, Tpp, add, mul = instance(p, B, m, k)
    t0 = time.time()
    L, bv = build(p, B, pts, Tpp, chunk=chunk, mode=mode)
    Q = L.qb
    st = Q.stats(); st['build_s'] = round(time.time() - t0, 2)
    print(f"[{mode}] p={p} ({p.bit_length()}b)  B={B}  G={G} order={order}  k={k}  T={T}")
    print(f"  QUBO: {st['vars']} vars, {st['couplers']} couplers, "
          f"AND={st['and_vars']} words={st['word_bits']} carries={st['carry_bits']}, "
          f"|J| range 2^{st['dynamic_range_bits']}, built in {st['build_s']}s")

    def energy_of(kk):
        bits = [(kk >> i) & 1 for i in range(m)]
        wv0 = {f"_b{i}": bits[i] for i in range(m)}
        try:
            x, _ = Q.witness({bv[i]: bits[i] for i in range(m)}, wv0)
        except Exception as ex:
            return None, str(ex)
        return Q.energy(x), None

    e, err = energy_of(k)
    print(f"  witness for the true k={k}: energy={e} err={err}")
    assert e == 0, "true solution is not a zero-energy state!"

    if exhaustive:
        zeros, degen = [], 0
        for kk in range(1 << m):
            ee, err = energy_of(kk)
            if ee is None: degen += 1; continue
            if ee == 0: zeros.append(kk)
            elif kk == k: assert False
        sol = [kk for kk in range(1 << m) if mul(kk, G) == T]
        print(f"  zero-energy k values: {zeros}   true solutions: {sol}   "
              f"degenerate (no witness): {degen}")
        assert set(zeros) <= set(sol), "SPURIOUS zero-energy state!"
        assert k in zeros
    return Q, bv, st

if __name__ == '__main__':
    run(97, 3, 5, 11)
    run(97, 3, 5, 11, mode='wallace')
