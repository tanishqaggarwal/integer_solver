#!/usr/bin/env python3
"""Agent Z: independent identification of the selector variables, by pure regex
on the raw file (no other agent's code, no pickles)."""
import re, os, json, collections

EQ = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'EQUATIONS.txt')
txt = open(EQ).read()
print("chars:", len(txt))

# every integer literal in the file, with its bit length
lits = collections.Counter()
for m in re.finditer(r'(?<![\dx_])(\d{20,})', txt):
    lits[m.group(1)] += 1
print("distinct >=20-digit literals:", len(lits))
bl = collections.Counter(int(k).bit_length() for k in lits)
print("bit-length histogram of big literals:", sorted(bl.items()))

# pattern A:  (x_S)*((x_C)-(BIGK))          -> s*(x - K)
pA = re.compile(r'\(x_(\d+)\)\*\(\(x_(\d+)\)-\((\d{60,})\)\)')
A = pA.findall(txt)
print("pattern s*(x-K):", len(A), "distinct:", len(set(A)))

selA = collections.Counter(a[0] for a in A)
coordA = collections.Counter(a[1] for a in A)
kA = collections.Counter(a[2] for a in A)
print("  distinct S:", len(selA), " distinct coord:", len(coordA), " distinct K:", len(kA))

# pattern B:  ((1)-(x_S))*(x_C)  in some parenthesisation
for pat in [r'\(\(\(1\)\)-\(x_(\d+)\)\)\*\(x_(\d+)\)',
            r'\(\(1\)-\(x_(\d+)\)\)\*\(x_(\d+)\)',
            r'\(x_(\d+)\)\*\(\(\(1\)\)-\(x_(\d+)\)\)',
            r'\(x_(\d+)\)\*\(\(1\)-\(x_(\d+)\)\)']:
    B = re.findall(pat, txt)
    print("pattern", pat, "->", len(B), "distinct", len(set(B)))

sels = sorted(int(s) for s in selA)
print("selector count:", len(sels))
print("occurrences per selector in pattern A:", sorted(collections.Counter(selA.values()).items()))
json.dump({'selectors': sels,
           'triples': sorted(set((int(a), int(b), c) for a, b, c in A))},
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'zsel.json'), 'w'))
print("wrote zsel.json")
