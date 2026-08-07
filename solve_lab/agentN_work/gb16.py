"""Dimension and degree of the exact ideal, via Singular; then the exact integer solve.

Reports, IN THIS ORDER:
  1. dim / degree / a standard basis of the collateral ideal I_C (16 generators, 16 unknowns);
  2. dim / degree of I_C + all 12 region rows;
  3. dim / degree of I_C + each 6-subset of the region rows (the subsets that would beat 39,026);
  4. only then, the integer solve on the variety.
"""
import os, sys, json, time, pickle, subprocess, itertools
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(2000000)
D = pickle.load(open(os.path.join(HERE, 'runs', 'kerquad2.pkl'), 'rb'))
d, gens, K0, Rl = D['d'], D['gens'], D['K0'], D['R']
V = ['s(%d)' % i for i in range(d)]


def poly_str(c0, L, Q):
    ts = []
    if c0:
        ts.append('(%d)' % c0)
    for a, c in enumerate(L):
        if c:
            ts.append('(%d)*%s' % (c, V[a]))
    for (a, b), c in sorted(Q.items()):
        ts.append('(%d)*%s*%s' % (c, V[a], V[b]))
    return '+'.join(ts) if ts else '0'


COL = [g for g in gens if not g[4]]
REG = [g for g in gens if g[4]]
Cstr = [poly_str(g[1], g[2], g[3]) for g in COL]
Rstr = [poly_str(g[1], g[2], g[3]) for g in REG]
Rlabels = [g[0] for g in REG]


def singular(script, tag, timeout=1800):
    p = os.path.join(HERE, 'runs', 'sing_%s.sing' % tag)
    open(p, 'w').write(script)
    t0 = time.time()
    try:
        r = subprocess.run(['Singular', '-q', '--no-warn', p], capture_output=True,
                           text=True, timeout=timeout)
        return r.stdout.strip(), time.time() - t0
    except subprocess.TimeoutExpired:
        return 'TIMEOUT after %ds' % timeout, time.time() - t0


RING = 'ring r = 0,(%s),dp;\n' % ','.join(V)
LIBS = 'LIB "primdec.lib";\n'


def dimdeg(ideal_lines, tag, extra=''):
    s = RING + LIBS + ideal_lines + '''
ideal G = std(I);
"dim="; dim(G);
"degree=";
if (dim(G) >= 0) { "  "; degree(G); }
"vdim="; vdim(G);
"ngens_gb="; size(G);
''' + extra + '\nquit;\n'
    return singular(s, tag)


def main():
    print('unknowns %d ; collateral generators %d ; region generators %d'
          % (d, len(Cstr), len(Rstr)), flush=True)

    print('\n### 1. collateral ideal I_C', flush=True)
    I_C = 'ideal I = ' + ',\n  '.join(Cstr) + ';\n'
    out, el = dimdeg(I_C, 'IC', extra='''
ideal Rad = radical(I);
"radical_gens="; size(std(Rad));
"radical_dim="; dim(std(Rad));
list pd = minAssGTZ(I);
"components="; size(pd);
for (int i=1; i<=size(pd); i++) { "  comp"; i; "  dim"; dim(std(pd[i])); }
''')
    print(out, flush=True)
    print('(%.1fs)' % el, flush=True)

    print('\n### 2. I_C + all 12 region rows', flush=True)
    I_all = 'ideal I = ' + ',\n  '.join(Cstr + Rstr) + ';\n'
    out2, el2 = dimdeg(I_all, 'IALL')
    print(out2, flush=True)
    print('(%.1fs)' % el2, flush=True)

    print('\n### 3. I_C + each 6-subset of the 12 region rows (924 subsets)', flush=True)
    lines = [RING, 'int nz=0; int ns=0;\n']
    lines.append('ideal C = ' + ',\n  '.join(Cstr) + ';\n')
    for i, s in enumerate(Rstr):
        lines.append('poly p%d = %s;\n' % (i, s))
    lines.append('list SS;\n')
    body = []
    for cnt, comb in enumerate(itertools.combinations(range(12), 6)):
        body.append('ideal J%d = C' % cnt + ''.join(',p%d' % j for j in comb) + ';\n')
        body.append('ideal G%d = std(J%d); if (dim(G%d)==-1) { nz=nz+1; } '
                    'else { ns=ns+1; "SOLVABLE_OVER_Q subset"; %d; dim(G%d); }\n'
                    % (cnt, cnt, cnt, cnt, cnt))
        body.append('kill J%d; kill G%d;\n' % (cnt, cnt))
    lines += body
    lines.append('"empty_variety="; nz; "nonempty="; ns;\nquit;\n')
    out3, el3 = singular(''.join(lines), 'SUB6', timeout=3600)
    print(out3[-4000:], flush=True)
    print('(%.1fs)' % el3, flush=True)

    json.dump(dict(IC=out, IALL=out2, SUB6=out3[-8000:]),
              open(os.path.join(HERE, 'runs', 'gb16.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
