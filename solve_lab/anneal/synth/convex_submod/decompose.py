#!/usr/bin/env python3
"""decompose.py -- H = H_sub (submodular, min-cut) + H_super (the AND core).

Two questions:
  (1) Which variables live ONLY in submodular couplers?  Conditioned on the rest
      those form a submodular subproblem and are determined by one min-cut -- they
      never need independent search.  Count them.
  (2) What is the supermodular 'core' (vars in >=1 positive coupler), and how does
      it compare to the true free dimension 2s (the operand bits)?

Also: an actual min-cut minimisation of H_sub on a moderate instance, to show it
is genuinely poly-time solvable, and how its minimiser relates to a joint ground
state (it does NOT recover one -- the coupling to H_super is what carries the
hardness)."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common
from collections import Counter, defaultdict


def split(Q):
    lin = {}
    sub = []   # (i,j,c<=0)
    sup = []   # (i,j,c>0)
    for m, c in Q.Q.items():
        if len(m) == 0:
            continue
        if len(m) == 1:
            lin[m[0]] = c
            continue
        (sub if c <= 0 else sup).append((m[0], m[1], c))
    return lin, sub, sup


def analyze(Q):
    lin, sub, sup = split(Q)
    in_super = set()
    for i, j, c in sup:
        in_super.add(i); in_super.add(j)
    pure_sub = [v for v in range(Q.n) if v not in in_super]
    return dict(n=Q.n, nsub=len(sub), nsup=len(sup),
                vars_in_super=len(in_super),
                vars_pure_sub=len(pure_sub),
                pure_sub_kinds=dict(Counter(Q.kind[v] for v in pure_sub)),
                super_kinds=dict(Counter(Q.kind[v] for v in in_super)))


def mincut_min(Q):
    """exact minimum of H_sub (linear + submodular couplers only) by max-flow
    min-cut.  Submodular quadratic c_ij<=0 x_i x_j = c_ij min via the standard
    two-terminal construction.  Returns the minimum energy of H_sub alone."""
    lin, sub, sup = split(Q)
    n = Q.n
    # node ids 0..n-1 vars, S=n, T=n+1
    S, T = n, n + 1
    cap = defaultdict(float)

    def add(u, v, w):
        if w:
            cap[(u, v)] += w

    const = Q.Q.get((), 0.0)
    # E = const + sum lin[i] x_i + sum c_ij x_i x_j  (c_ij<=0)
    # x_i=1 <=> i on T side. Encode min via source/sink capacities.
    # For a submodular pairwise -w x_i x_j (w>=0): rewrite as
    #   -w x_i x_j = -w x_i + w x_i (1 - x_j) ... use standard reduction:
    # -w x_i x_j (w>=0) is supermodular-negative? c_ij<=0 so term = c_ij x_i x_j,
    # let w=-c_ij>=0: term = -w x_i x_j. This is SUBMODULAR (favours x_i=x_j=1
    # both or splitting?).  -w x_i x_j is minimised (=-w) at x_i=x_j=1: it is a
    # submodular function; graph edge S->i or i->j.  Use the reduction:
    #   -w x_i x_j = -w + w(1-x_i) + w(1-x_j) - w(1-x_i)(1-x_j)
    # simpler known gadget: theta(x_i,x_j) with theta(1,1)=-w else 0 is
    # submodular (theta00+theta11=-w <= theta01+theta10=0).  Add via:
    #   const += theta00
    #   S->i : theta10-theta00 ; i->T handled by lin ; edge i->j: theta01+theta10-theta00-theta11
    # We accumulate lin-like terms then build.
    theta_const = 0.0
    src = defaultdict(float)  # capacity S->i  (cost of x_i=0? ) we follow Kolmogorov
    snk = defaultdict(float)  # capacity i->T
    edge = defaultdict(float)
    # unary: cost a_i for x_i=1
    for i, a in lin.items():
        if a >= 0:
            snk[i] += a       # pay a when x_i=1 (i on T side): edge i->T cap a
        else:
            src[i] += -a      # pay |a| when x_i=0
            theta_const += a
    for i, j, c in sub:
        w = -c  # w>=0, term = -w x_i x_j, theta11=-w
        # theta00=0,theta01=0,theta10=0,theta11=-w
        # submodular check: theta00+theta11 (=-w) <= theta01+theta10 (=0) OK
        # Kolmogorov reduction:
        # const += theta00 (0)
        # snk[i] += theta10-theta00 = 0
        # snk[j] += theta11-theta10 = -w  -> negative cap, handle by moving const
        # edge[i,j] += theta10+theta01-theta00-theta11 = w
        # Because theta11-theta10=-w<0, add to src side instead:
        A = 0.0            # theta10-theta00
        Bc = -w            # theta11-theta10
        # split B: if B<0 add -B to src[j] and const+=B
        snk[i] += A
        if Bc >= 0:
            snk[j] += Bc
        else:
            src[j] += -Bc
            theta_const += Bc
        edge[(i, j)] += w
    # build residual graph max-flow (simple BFS Edmonds-Karp) on nodes 0..n+1
    from collections import deque
    graph = defaultdict(lambda: defaultdict(float))
    for i, w in src.items():
        graph[S][i] += w
    for i, w in snk.items():
        graph[i][T] += w
    for (i, j), w in edge.items():
        graph[i][j] += w
    # Edmonds-Karp
    def bfs():
        par = {S: None}
        q = deque([S])
        while q:
            u = q.popleft()
            for v, w in graph[u].items():
                if w > 1e-12 and v not in par:
                    par[v] = u
                    if v == T:
                        return par
                    q.append(v)
        return None
    flow = 0.0
    while True:
        par = bfs()
        if not par:
            break
        # bottleneck
        b = float('inf'); v = T
        while par[v] is not None:
            u = par[v]; b = min(b, graph[u][v]); v = u
        v = T
        while par[v] is not None:
            u = par[v]; graph[u][v] -= b; graph[v][u] += b; v = u
        flow += b
    mincut_energy = const + theta_const + flow
    return mincut_energy, flow


def main():
    sizes = [8, 16, 32, 64, 128, 256]
    print("SUBMODULAR-CORE DECOMPOSITION  H = H_sub + H_super")
    print(f"{'s':>4} {'n':>8} {'2s(free)':>9} {'vars_in_super':>14} "
          f"{'vars_pure_sub':>14} {'pure/n':>7}")
    rows = []
    for s in sizes:
        p = common.real_p256() if s == 256 else common.prime_of_bits(s)
        Q = common.build_mm(s, p=p)['Q']
        a = analyze(Q)
        a['s'] = s
        rows.append(a)
        print(f"{s:>4} {a['n']:>8} {2*s:>9} {a['vars_in_super']:>14} "
              f"{a['vars_pure_sub']:>14} {100*a['vars_pure_sub']/a['n']:>6.1f}%")

    print("\npure-submodular variable kinds (min-cut-determined given the core):")
    for a in rows:
        print(f"  s={a['s']:>4}: {a['pure_sub_kinds']}")
    print("\nsupermodular-core variable kinds:")
    for a in rows:
        print(f"  s={a['s']:>4}: {a['super_kinds']}")

    # actual min-cut of H_sub on a moderate instance
    print("\nMIN-CUT solve of H_sub alone (H_super dropped) on s=32:")
    Q = common.build_mm(32, p=common.prime_of_bits(32))['Q']
    e, flow = mincut_min(Q)
    print(f"  H_sub min energy = {e:.1f} (poly-time via one max-flow, |flow|={flow:.1f})")
    print("  -- this minimises the ferromagnetic part in isolation; it does NOT")
    print("     yield a joint ground state, because every operand bit also sits in")
    print("     a positive (supermodular) AND coupler that this relaxation ignores.")

    with open(os.path.join(os.path.dirname(__file__), 'decompose.json'), 'w') as f:
        json.dump(rows, f, indent=2)


if __name__ == '__main__':
    main()
