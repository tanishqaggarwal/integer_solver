# AE / item 3 — does U's partition theorem kill §8 for all `|S|` at once?

**Agent AE**, answering the coordinator's mid-task question. Script: `ae_support.py`;
outputs `res_support.json`, `res_927_supports.json`.

---

## 0. FIRST, A SCOPE FLAG YOU SHOULD READ BEFORE THE RESULT

The message says *"the crux I handed you stands"* and *"your object — the 927 lift conditions and
766 off-pins"*. **No such crux was ever handed to me and that is not my object.** My brief was
the cheap structural sweeps on `k₀` (magnitude / quotient / orbit families on the curve), which
is what `STRUCTURED_KEYS.md` reports. I think this part of the message was routed to the wrong
agent.

I answered it anyway, because the question turned out to be short and because the answer is
worth having. But **the measurement below leans on three other agents' artefacts** (each
sanity-checked here first), which is exactly the kind of dependency an independent derivation is
supposed to avoid — so weigh it accordingly, and if you want it re-derived from `EQUATIONS.txt`
by a parser that shares nothing, that is a different and much longer job than this was.

---

## 1. THE ANSWER, UP FRONT

> **U's partition theorem does not kill the growth mechanism, and the reason it does not is the
> opposite of the one proposed. The conjectured chain was "no wraparound ⇒ the fold is exact
> integer addition ⇒ the lift conditions are `|S|`-blind". The first implication is right. The
> second is backwards: exact integer addition makes the fold value at a node an *injective*
> encoding of `S` restricted to that node, so the conditions there see `S` *perfectly*, not not
> at all.**
>
> **What actually limits §8 is a different fact, which I measured: `925` of the `927` conditions
> have a selector support that is a proper subset of `{0..255}`, and `745` of them see `4`
> selectors or fewer. Only `2` of the `927` see all `256`. But support locality alone does not
> close §8 either, because the supports are laminar and a *conjunction* of subtree-local
> conditions can still bound `|S|` globally.**
>
> **So: the crux resolves as "merely bounds it" — but it bounds it to two named conditions.**

## 2. WHY THE PROPOSED IMPLICATION RUNS BACKWARDS

Let `v` be a node of the fold tree with selector support `σ_v`, and let
`m_v = Σ_{i ∈ S ∩ σ_v} 2^i` be the scalar it carries.

* If the fold **wrapped** mod `N`, the map `S ∩ σ_v ↦ m_v mod N` would be many-to-one and the
  node would genuinely lose information about which leaves are on.
* U's theorem says it never wraps: `maskval(σ_v) < N` for every proper support (max
  `0.798718631·N`). Therefore `m_v` is the **exact integer** `Σ_{i∈S∩σ_v} 2^i`, and because the
  summands are distinct powers of two, **binary representation is unique**, so
  `S ∩ σ_v ↦ m_v` is **injective**.

A condition sitting at `v` is a predicate on the wire values there, which are functions of `m_v`
(and of the free variables). Injectivity means no information about `S ∩ σ_v` has been destroyed
before the condition sees it. **No-wraparound therefore preserves the conditions' sensitivity to
`S`; it is the *enabling* hypothesis for §8's mechanism, not its refutation.** U's theorem is
still doing exactly the work U claimed for it — it closes the degenerate-gadget escape hatch —
but it is not a `|S|`-blindness result and should not be cited as one.

## 3. WHAT I MEASURED INSTEAD — the selector support of each of the 927

### 3.1 Dependencies, each sanity-checked before use

| artefact | source | check I ran | result |
|---|---|---|---|
| per-wire selector support | `agentU_work/v_supp2.pkl` | distinct non-empty supports; laminarity; root; singletons | **511** distinct, **0** laminarity violations, exactly **1** support of size 256, its two maximal proper subsupports are **178 / 78**, exactly **256** singletons — reproduces U's table |
| 39,033-atom parse | `agentT_work/mirror/F/circ4.pkl` | 38,748 distinct wires; U's table covers 34,999 of them, and **0** U keys are absent from F | consistent |
| handle / cofactor list | `agentT_work/mirror/L/handles.pkl` | 3,681 cofactors, all parsed | `c == 1`: **2,754**, `c > 1`: **927** — **exactly** L's, P's and T's counts, rebuilt here |

### 3.2 The one thing that could have invalidated the measurement, closed

Of the wires in a condition's atoms, `910` lookups had no entry in U's support table. **All of
them are the `P` wire** in the definition `h = P·u` — and there are only **55 distinct** such
wires across all 927 conditions. Following copy-chains in F's parse, **all 55 resolve to the
literal 256-bit constant `p`** (`= 2^256 − 2^32 − 977`, checked by value). A constant has no
selector dependence, so **the histogram below is exact, not a lower bound.** Had even one `P`
carried the root support, every number in §3.3 would have collapsed to 256.

### 3.3 The histogram

| `\|σ\|` | conditions | | `\|σ\|` | conditions |
|---|---|---|---|---|
| **0** | **48** | | 21–50 | 20 |
| **1** | **420** | | **78** | **2** |
| **2** | **189** | | **88** | **2** |
| **3** | **76** | | **90** | **1** |
| **4** | **60** | | **178** | **2** |
| 5–20 | 105 | | **256** | **2** |

Cumulative: `|σ| ≤ 1`: **468 (50.5 %)** · `≤ 2`: 657 (70.9 %) · `≤ 4`: **793 (85.5 %)** ·
`≤ 16`: 887 (95.7 %) · `≤ 178`: **925 (99.8 %)** · `= 256`: **2**.

Two further structural facts, both exact:

* **48 conditions have empty selector support** — they do not depend on the ON-set at all. They
  are constant predicates and can contribute nothing to any `|S|` constraint, in either
  direction.
* The `879` conditions with non-empty support carry **exactly 511 distinct supports — set-equal
  to U's entire 511-node tree.** Every node of the fold tree carries at least one `c > 1`
  condition. The conditions are not clustered at the leaves or at the root; they are spread
  across the whole tree, one or more per node.

## 4. THE LOCALITY LEMMA, AND EXACTLY WHAT IT BUYS

**Lemma.** A condition `C` whose atoms' wires have combined selector support `σ` satisfies
`C(S) = C(S′)` whenever `S ∩ σ = S′ ∩ σ`.
*Proof.* Immediate from the definition of selector support: every wire feeding `C` is a function
of `S` only through `S ∩ σ`. ∎

**Corollary (per-condition).** If `C` is satisfiable at all, then for any witness `S₀` and any
`S` with `S ∩ σ = S₀ ∩ σ`, `C(S)` also holds. Hence a single condition with `|σ| = s` is
consistent with every `|S|` in a window of width `256 − s`. Applied to the histogram:

* **925 of 927 conditions are individually consistent with `|S|` ranging over ≥ 78 values.**
* **793 of 927 are individually consistent with `|S|` ranging over ≥ 252 values** — i.e. with
  essentially the whole range.
* **Only the 2 root-support conditions can, on their own, see `|S|` at all.**

**And what it does NOT buy — stated plainly, because this is where the stronger claim dies.**
The corollary is per-condition. §8 asks about the **conjunction**. A conjunction of subtree-local
predicates over a laminar family *can* bound cardinality globally: 256 singleton-support
conditions each forcing its own selector off would force `|S| = 0`, and every one of those
conditions is "local". So **locality is necessary but not sufficient**, and the clean result the
coordinator hoped for — "every condition has small support, therefore §8 dies for all `|S|` at
once" — **does not follow from this measurement.** I am not going to dress it up as though it
does.

What the measurement *does* do is collapse the search: any instance-side `|S|` constraint must
either (a) be carried by the two root-support conditions, or (b) arise from the joint action of
many local conditions, which is a statement about their satisfying sets `A_C ⊆ 2^{σ_C}`, not
about their supports.

## 5. THE TWO CONDITIONS THAT COULD CARRY §8 — named, for routing

These are the only two of the 927 whose selector support is all of `{0..255}`:

| cofactor `u` | handle `h` | multiplier `c` | guard atom (F's parse) |
|---|---|---|---|
| `x5146` | `x29804` | `6672769` | `((x15298*x25739)-(6672769*x29804))` |
| `x14393` | `x34243` | `12354891` | `((x24468-x13682)-(12354891*x34243))` |

The next tier down, for completeness (support 178 — the larger root half):

| `x14485` | `x9216` | `2264251` | `((x38170*x15286)-(2264251*x9216))` |
| `x26020` | `x9756` | `11103619` | `((x28258*x21589)-(11103619*x9756))` |

**The obvious next question — do these two conditions actually vary with `|S|`? — is a probe:
it requires evaluating the guard at chosen ON-sets. Per your coordination note that is agent T's
and agent AH's ground, so I did not run it. Routing it to them is the single highest-value
follow-up from this file**, and it is now a two-condition question rather than a 927-condition
one. If both turn out to be `|S|`-insensitive, §8 loses its only individually-capable carriers
and what remains is the harder joint question in §4.

## 6. WHAT IS NOT CLAIMED

* No infeasibility claim, and no claim about `w`.
* The histogram is a fact about **F's parse + L's handle list + U's support closure**, jointly.
  All three reproduce their published counts here, and the `P`-wire audit closes the one gap I
  could find, but this is a **join of three artefacts, not a fourth independent derivation.**
* §8 is **not** settled by this file, in either direction.
