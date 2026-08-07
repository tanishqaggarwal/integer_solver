"""Agent B model v3: fully hash-consed expression DAG over the whole file.

Node encodings (all hash-consed into self.nodes):
  ('C', k)              integer constant k
  ('V', i)              variable x_i
  ('S', ((c1,n1),...))  sum of coef*node, canonical sorted, constants folded
  ('P', (n1,n2,...))    product of nodes (sorted), no constants inside

Roots: one node per equation (the LHS).
"""
import sys, time, pickle, collections
from bparse import parse_line

class DAG:
    def __init__(self):
        self.nodes = []          # list of node tuples
        self.idx = {}            # node tuple -> id
        self.C0 = self.mk(('C', 0))
        self.C1 = self.mk(('C', 1))

    def mk(self, t):
        i = self.idx.get(t)
        if i is None:
            i = len(self.nodes); self.idx[t] = i; self.nodes.append(t)
        return i

    def const(self, k):
        return self.mk(('C', k))

    def var(self, i):
        return self.mk(('V', i))

    def summ(self, terms):
        """terms: list of (coef, node_id).  Returns node id."""
        acc = collections.defaultdict(int)
        k = 0
        for c, n in terms:
            if c == 0: continue
            t = self.nodes[n]
            if t[0] == 'C':
                k += c * t[1]
            elif t[0] == 'S':
                # inline nested sums
                for c2, n2 in t[1]:
                    acc[n2] += c * c2
                # nested sums carry their constant as a ('C',k) member already
            else:
                acc[n] += c
        items = tuple(sorted((n, c) for n, c in acc.items() if c))
        items = tuple((c, n) for n, c in items)
        if k:
            items = items + ((k, self.C1),)
        if not items:
            return self.C0
        if len(items) == 1 and items[0][0] == 1:
            return items[0][1]
        return self.mk(('S', items))

    def prod(self, facs):
        k = 1
        ns = []
        for n in facs:
            t = self.nodes[n]
            if t[0] == 'C':
                k *= t[1]
                if k == 0: return self.C0
            elif t[0] == 'P':
                ns.extend(t[1])
            else:
                ns.append(n)
        if not ns:
            return self.const(k)
        ns = tuple(sorted(ns))
        pid = ns[0] if len(ns) == 1 else self.mk(('P', ns))
        if k == 1:
            return pid
        return self.summ([(k, pid)])

    def build(self, ast):
        k = ast[0]
        if k == 'n': return self.const(ast[1])
        if k == 'v': return self.var(ast[1])
        if k == '+':
            return self.summ([(1, self.build(c)) for c in ast[1]])
        if k == '*':
            return self.prod([self.build(c) for c in ast[1]])
        raise ValueError(k)


def main():
    D = DAG(); roots = []
    t0 = time.time()
    for i, line in enumerate(open('/home/user/integer_solver/EQUATIONS.txt')):
        if not line.strip(): continue
        roots.append(D.build(parse_line(line)))
        if (i+1) % 10000 == 0:
            print("  %d nodes=%d %.1fs" % (i+1, len(D.nodes), time.time()-t0), flush=True)
    print("eqs=%d nodes=%d %.1fs" % (len(roots), len(D.nodes), time.time()-t0))
    kc = collections.Counter(t[0] for t in D.nodes)
    print("node kinds", dict(kc))
    pickle.dump({'nodes': D.nodes, 'roots': roots}, open('dag.pkl','wb'), -1)
    print('wrote dag.pkl')

if __name__ == '__main__':
    main()
