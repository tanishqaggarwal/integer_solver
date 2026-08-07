"""W: symbolically recover ONE law block verbatim from EQUATIONS.txt.  My own expansion.
No import of any other agent's file."""
import sys, os, re, json, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import model, sympy as sp

d = model.get(); A = d['atom_src']; AV = d['atom_vars']; EQ = d['eq_terms']
fw = pickle.load(open('fwd2.pkl','rb'))
definer = fw['definer']            # var -> atom index that defines it (or -1)

# which atoms mention which var
byvar = {}
for i, vs in enumerate(AV):
    for v in vs: byvar.setdefault(v, []).append(i)

P_VAR, Q_VAR = 26064, 24453
PVAL = 115792089237316195423570985008687907853269984665640564039457584007908834671663
QVAL = 97553848499418123410591666447050222001188385549510401465815187079080512838891

Ecands = [(i, s) for i, s in enumerate(A) if re.fullmatch(r'x_\d+ - \(x_\d+ \+ x_24453\)', s)]
print('E atoms:', len(Ecands))
i0, s0 = Ecands[0]
print('block seed atom', i0, s0)

SYM = {}
def sym(v):
    if v not in SYM: SYM[v] = sp.Symbol('x%d' % v)
    return SYM[v]
def parse(s):
    ns = {}
    for m in set(re.findall(r'x_(\d+)', s)): ns['x_%s' % m] = sym(int(m))
    return sp.sympify(s, locals=ns, evaluate=True)

# expand a var backwards through definers, stopping at `stop` set
def expand(expr, stop, maxdepth=40):
    for _ in range(maxdepth):
        fv = [v for v in expr.free_symbols]
        subs = {}
        for s_ in fv:
            v = int(str(s_)[1:])
            if v in stop: continue
            a = definer[v]
            if a < 0: continue
            src = A[a]
            m = re.fullmatch(r'x_%d - (.*)' % v, src)
            if not m: continue
            subs[s_] = parse('(' + m.group(1) + ')')
        if not subs: break
        expr = expr.subs(subs)
    return sp.expand(expr)
if __name__ == '__main__':
    # the block:  E = D + Q.  Find the atoms downstream of E.
    Evar = int(re.match(r'x_(\d+)', s0).group(1))
    Dvar = int(re.findall(r'x_(\d+)', s0)[1])
    print('E var x_%d, D var x_%d' % (Evar, Dvar))
    print('D def:', A[definer[Dvar]] if definer[Dvar] >= 0 else None)
    print('consumers of E:', [(a, A[a][:200]) for a in byvar[Evar] if a != i0])
