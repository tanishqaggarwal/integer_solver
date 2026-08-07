"""Is the co-occurrence in those 48 atoms a genuine COUPLING or just coarse BUNDLING?

Sharp test: a genuine cross-selector coupling requires two selectors to meet inside a single
nonlinear term (a product), or to sit in a cardinality/one-hot sum.  Bundling -- several
independent per-selector load constraints packed into one additive atom because my parse is
coarser than P's -- puts each selector in its own additive term and couples nothing.

So for each atom with >=2 selectors:
  - does it contain a PRODUCT of two distinct selectors?          -> genuine coupling
  - is it a booleanity certificate (sum of x*(1-x) / x*x - x)?    -> not a coupling
  - otherwise: are the selectors in distinct top-level additive terms?  -> bundling
"""
import sys, re, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentS_work')
import common as C
import harness as H

SEL = set(f for f in C.cluster_cone() if C.isbool(f))
multi = []
for a, t in enumerate(H.atoms):
    if not t:
        continue
    vs = {int(m.group(1)) for m in re.finditer(r'x_(\d+)', t)}
    s = vs & SEL
    if len(s) >= 2:
        multi.append((a, s, t))
print("atoms with >=2 selectors: %d" % len(multi), flush=True)

prodpat = re.compile(r'x_(\d+)\s*\*\s*x_(\d+)')
kinds = collections.Counter()
coupled = []
for a, s, t in multi:
    # genuine product of two DISTINCT selectors?
    prods = [(int(x), int(y)) for x, y in prodpat.findall(t)]
    real = [(x, y) for x, y in prods if x in SEL and y in SEL and x != y]
    # booleanity certificate: every selector occurs only as x*(1-x) or x*x - x
    boolish = all(re.search(r'2 \* x_%d \* \(1 - x_%d\)|x_%d \* x_%d - x_%d' % (v, v, v, v, v), t) for v in s)
    if real:
        kinds['PRODUCT-COUPLED'] += 1
        coupled.append((a, real, t[:160]))
    elif boolish:
        kinds['booleanity certificate'] += 1
    else:
        kinds['bundled (separate additive terms)'] += 1
print("\nclassification:", dict(kinds), flush=True)
if coupled:
    print("\n*** GENUINE PRODUCT COUPLINGS BETWEEN DISTINCT SELECTORS ***", flush=True)
    for a, real, t in coupled[:12]:
        print("   a%-6d pairs=%s\n      %s" % (a, real, t), flush=True)
else:
    print("\nNo atom multiplies two distinct selectors: every co-occurrence is either a\n"
          "booleanity certificate or independent per-selector terms bundled into one atom.", flush=True)

# also: does any selector-pair share a single additive term at all?
print("\nselector self-products (booleanity, expected):",
      sum(1 for a, s, t in multi for x, y in prodpat.findall(t) if int(x) == int(y) and int(x) in SEL), flush=True)
