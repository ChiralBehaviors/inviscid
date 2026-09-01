# SU(2) boundary conditions for the jitterbug medium — brainstorm record

2026-08-31. Status: **resolved 2026-09-01** — the §8.6 gate plan ran in full
(§9); the closure question is answered and decided (OWNER DECISION 22: the
double). Origin: owner question "how do we deal with boundary conditions — not
a torus, but a double covering of 3D space," clarified to mean **the
properties of SU(2)**.

> **§1–7 are the brainstorm as first written; §8 (five critique lenses) revises
> several of their claims — read §8 before citing anything above it. §9 records
> what the gates measured and the closure decision.**

## 1. The question

Finite realizations of the medium need a closure. The record already shows the
stakes: free boundaries tumble and wind (the three-cell animation,
inviscid-qvf.26); fixed-centre rigs manufacture a lock the medium does not have
(DECISION 20, the retired rig_lock line). Periodic identification (the 3-torus)
is the unexamined default. The proposal on the table: the correct closure
involves the SU(2) → SO(3) double cover.

## 2. The constraint that shapes any answer

R³ is simply connected, so it has **no nontrivial connected unbranched covers**.
"The medium is a double cover of 3-space" can only mean:

1. a **branched** double cover (two sheets glued along a branch locus), or
2. a cover of a **quotient**: the medium closes on a flat 3-manifold/orbifold
   that is not the torus, and the faithful medium is that quotient's connected
   double cover. The identification is then a translation composed with a
   nontrivial deck transformation; the torus is the case where the deck map is
   the identity.

Either way the proposal is only sharp once the deck transformation is named.

## 3. Three deck-map candidates already in the model

**(a) The exchange screw.** Measured (ring.py, gated): a half-lattice
translation carries a cell site to a void site and the void runs at
`b = a + 60°`. So τ = (translate by (1,1,1)) ∘ (fold +60°) ∘ (cell↔void) is a
symmetry of the family. An even-period torus represents it; an odd period
breaks it unless the screw is built into the identification — and then the
faithful cell/void medium is a connected double cover of the screw quotient.
This lands in the classified list of flat closed 3-manifolds.

**(b) The inside-out cover.** `L(a+180°) = −L(a)` (the spacing law, measured);
the far half of the overdrive cycle looks like the array inverted through its
centre. Conjectured deck map: central inversion ∘ fold+180°. Shape space would
have period 180° with the identity-tracking medium double-covering it; the two
collapses (+60°, +240°) are where the sheets touch in ambient space.
**Open: `X(a+180°) = −X(a) + const` has not been verified on configurations.**

**(c) The two-sided-plate reading.** The voids are the backs of the cells' own
plates; emptying them worked because they carry no independent DOF ("empty" =
"not independent — the other sheet"). Formally: one plate complex, two ambient
sheets glued along it — a branched double cover of R³ with branch locus the
plate complex and deck map the cell↔void exchange. Finite version: close a
block not into a torus but into its **double** — glue a second copy along the
boundary plate-backs, closing every half-weld, leaving no boundary and no rig.

## 4. The SU(2) sharpening

What SU(2) adds beyond §3(b): π₁(SO(3)) = ℤ₂ — a 2π rotation is a
noncontractible loop, 4π contracts. The overdrive measured that **positions**
close at 360° of fold (`X(300°) = X(−60°)` to 1e-15, gate R5). The SU(2)
question is whether the **rotation history** closes too.

Hypothesis H1: it does not — the cycle is the nontrivial loop.

- Each plate spins about its face normal as the fold advances (the classic
  jitterbug rotates faces through the VE→octa leg). If plate spin tracks the
  fold roughly linearly, one 360° cycle gives each plate a full 2π about its
  normal: positions restored, quaternion lift returns as **−q**.
  **Open: the per-cycle plate rotation has NOT been measured; the claim's sign
  depends on it.** Note a caution from the record: the coherent motion has
  `turn = a − b + 60° = 0` across a face weld (face_to_face law with
  `b = a + 60`), i.e. **relative** rotation across shared faces vanishes — if
  twist accumulates anywhere it is at the **vertex joints** (corners of plates
  of different bodies) or between a plate and its body frame, not at the face
  welds.
- Micro-scale ℤ₂, measured: the two joints on one vertex point swap by passing
  head-on (Two Joints, One Point; separation directions anti-parallel through
  the VE). Two swaps per cycle is the **full twist**: trivial as a permutation,
  nontrivial as a braid.

**Where the ℤ₂ physically lives — the actual modeling decision.** As coded,
welds and ties are point identifications, and a point remembers no winding.
That model *quotients out* the ℤ₂: it is the SO(3) medium by construction. A
rubber joint with physical extent remembers twist: 2π leaves it twisted, 4π
relaxes (belt trick; the anti-twister mechanism is the statement that a body
tethered to its surroundings can rotate indefinitely at 4π periodicity). A
medium that is all tether is a candidate mechanical realization. So:

- **SO(3) medium** (current): joints are points; one overdrive cycle is the
  identity; no spin sector.
- **SU(2) medium** (proposal): joints carry twist mod 2; one cycle restores
  every position but leaves every joint 2π-twisted — not the identity; 720°
  closes. Consequences if true: a ℤ₂ quantum number; an antiperiodic (spin)
  sector in the Bloch analysis alongside the periodic one; twist defects (an
  unpaired 2π joint) as transportable excitations.

Under this reading the "double cover" is not over position space at all: it is
over the orientation/configuration bundle, and "torus vs double cover" becomes
"does the state space include the spin lift or not."

## 5. Proposed gates (none run)

- **G1 (decisive, cheapest):** continuous quaternion lift of one plate's
  orientation over one 360° overdrive cycle; read the sign.
  `q(300°) = −q(−60°)` ⇒ H1; `+q` kills it. Plate frames are rigid throughout
  (the collapse degenerates positions, not plate frames), so the lift is
  well-defined.
- **G2:** accumulated relative rotation across one **vertex joint** per cycle
  (the two plates whose corners it ties), and its lift sign. This is where
  twist would physically sit, given the face-weld turn law is 0.
- **G3:** braid class of the joint-pair swap per cycle (full twist vs trivial).
- **G4:** `X(a+180°) = −X(a) + const` on a patch — pins deck map (b).
- **G5:** spectrum of a doubled block (§3c) vs free block vs torus — does the
  double lose the seven free-boundary modes and converge to the Bloch bands
  fastest?
- **G6:** twisted Bloch sector — dispersion with a 60°-screw identification
  along (1,1,1) (deck map (a)) vs periodic.

## 6. Open questions

- Composition: is the exchange (+60°) a "sixth root" of the spinor flip — six
  exchanges = 2π? The collapses sit at +60° and +240° = 60°+180°; suggestive,
  unverified, and possibly numerology.
- Does the physical 60° arc (all real struts allow) ever exercise the ℤ₂, or
  is it overdrive-only — constraining the mathematics of the medium but not
  its physical motion?
- Which closure should actual finite simulations (Necronomata, the honeycomb
  animations) use: free + contact (current), torus, screw quotient, or the
  double of the block?

## 7. Terminology hygiene

"Double covering" collides with `model/double_covering/` (cells on both
sublattices — the retired line). If this direction proceeds it needs its own
name: **spin cover** (§4), **exchange cover** (§3a), **plate double** (§3c).

## 8. Critiques (five lenses, 2026-08-31)

**Read this section before citing anything in §1–7.** Sections above are kept
verbatim as the record of the brainstorm; several of their claims are revised
here. (All five critics also note: T2 write-back was unavailable this session —
this file is the durable record.)

### 8.1 Mathematical rigor

- **§3c is the wrong object as named.** A branched cover needs a
  codimension-≥2 branch locus; the plate complex is codimension 1 — plate
  interiors double trivially. The correct object: a **chamber-swap double
  cover branched over the 1-skeleton** of the plate arrangement (its edges and
  vertices).
- **§4 overclaims "the medium's cycle is THE nontrivial loop."** Nothing
  computes π₁ of the configuration space; the −q inference is sound only via
  the pullback along one plate's map to SO(3), and it presumes the cycle is a
  loop in the *full* configuration space (orientations included) — which is
  part of what G1 must show, so the prose was circular.
- **§3a's τ is not a deck transformation of a spatial quotient** — fold+60°
  moves between fibers of the family. Correct home: a **mapping torus over the
  fold circle**, valid only if six applications reduce to a pure lattice
  translation — asserted, never checked.
- The 60+180=240 composition (§6) is confirmed content-free: no shared group
  has been constructed for the two deck maps to compose in.

### 8.2 Mechanics and prior art

- **The central mechanism fails: a bonded point-joint cannot remember twist
  mod 2.** ℤ₂ memory needs tether topology — slack and room to route around
  the rotating body. The model's joint is a torsion bushing: strain energy
  smooth and monotonic in twist, direct untwist always available. And a fully
  face-welded, vertex-tied lattice has no slack anywhere for the 4π maneuver —
  the anti-twister citation is decorative.
- Deflate the vocabulary: the real classical content of the "spin sector" is
  **antiperiodic Bloch boundary conditions** — no exchange antisymmetry, no
  quantum number.
- Prior art exists and should be cited: Kane–Lubensky topological mechanics
  (different invariant); the **spinor linkage** (Gallagher & Weiss 2022 — what
  a genuinely tether-capable joint requires); **Möbius kaleidocycles**
  (Schönke & Fried, PNAS 2019 — real linkages with twisted closure, via
  engineered hinge angles and slack).
- What survives untouched: the purely kinematic claims — lift signs and
  spectra are facts about paths, not about material remembering them.

### 8.3 Consistency with the measured record

Nothing in §1–7 contradicts the record; confidence was miscalibrated in
**both directions**:

- **G1's answer is already predicted.** `jitterbug.py`'s closed form is exact
  in the fold angle; `one_cell.py` R2 gates spin(a) = a to 1e-9; body
  quaternions sit at identity throughout the drive. One 360° cycle is
  therefore 2π of plate spin: **q(300°) = −q(−60°) is a near-corollary of
  gated facts** — G1 is a cheap confirmation, not a discovery.
- **G4 is near-formality**: L(a+180)=−L(a) is an algebraic identity of the
  gated spacing law, and the R3 front census flips 100%→0% under a+180.
- **§3a overstates**: the Bieberbach/flat-manifold claim needs a free,
  properly-discontinuous action and injective tiling by iterated τ — nothing
  gated checks either.
- **"Two swaps per cycle" overstates**: only the cell-type VE crossing swap is
  gated (vertex_point R4, ±15° window); the void-type swap is inferred by
  symmetry, never measured.
- Bonus fact for §6: the gated crossing counts (R4: 108–1620 mid-passage)
  are **hard evidence that real struts cannot traverse the loop** — see 8.5.

### 8.4 Experiment design

- G1 needs **null controls on the lift code itself** (trivial loop → +q, bare
  2π → −q, bare 4π → +q) and far finer than 30° sampling through the
  inside-out region and the collapses.
- **−q is necessary, not sufficient**, for "the medium is SU(2)" — the
  point-joint model as coded structurally registers no winding regardless.
- **G2 must measure both co-located joints** ({O,A} and {B,D}) — "the vertex
  joint" is not singular.
- **G3 is ill-posed**: the exchange is collinear (separation directions
  anti-parallel), no transverse framing to braid in. Drop it or fold into G2.
- Cost order: validate lift code → G1 (multiple plates, patches) → G4
  (free — data already in `overdrive.drive()`) → G2 → then, only if wanted,
  G5/G6 (both need assembly constructors that don't exist).
- **Missing gate**: one restricted to the physical 60° arc — §6 bullet 2 made
  executable.

### 8.5 Epistemics and method

- **Scope drift is the strongest objection overall**: the originating question
  was closure for finite simulations; the document answers a different,
  grander question. Four of six gates serve the spin question; only G5/G6
  serve closure selection — and those are exactly the two needing unbuilt
  machinery. A document answering the actual question would lead with them.
- **Even a successful G1 may be physically inert**: by §6's own admission —
  now hardened by 8.3's crossing counts — the noncontractible loop lives
  entirely in the analytic continuation; the physical 60° arc is a
  contractible path. The ℤ₂ is (pending G1) a true fact about the
  parametrization, not established as a fact about the physical medium.
- §4 asserts its reframe declaratively two paragraphs after flagging its
  premise unmeasured (an inversion 8.3 partly repairs: the premise was in the
  record all along). This project's history (rig_lock, the retracted
  exchange) shows unhedged prose outliving its basis.
- Deflate "spin"-family vocabulary throughout to **classical ℤ₂ holonomy /
  antiperiodic sector**; §7's name for the §3c object should follow 8.1
  ("chamber-swap cover", not "plate double").

### 8.6 Synthesis

What the five lenses agree on:

1. **The SU(2) intuition is (very likely) correct kinematics** — the record
   already implies each plate's per-cycle path is the nontrivial SO(3) loop.
   Run the controls + fine-grained G1 to close it, expecting −q.
2. **It is not correct mechanics as proposed** — point joints cannot store
   the ℤ₂; storing it would require redesigned joints (framed
   tether/ribbon, cf. the spinor linkage), a modeling *decision* to be argued
   on its own merits, not a property the current medium has.
3. **It does not answer the boundary-condition question.** The physical arc
   never traverses the loop (crossings forbid it), so closure selection for
   finite simulations rests on the G5/G6 lane — doubled block vs screw
   quotient vs torus — which is where the machinery investment belongs if the
   BC question is the goal. The τ-quotient additionally owes the g⁶ and
   properly-discontinuous checks (8.1, 8.3) before it is even a candidate.
4. Corrected revised gate plan: **G0** lift-code controls; **G1** fine-grained,
   multi-plate (confirmation); **G4** free-ride; **G2** both joints; G3
   dropped; **G7** physical-arc gate (new); **G5/G6** as the actual BC lane,
   gated on constructors and the τ checks.

## 9. Resolution (2026-09-01)

The §8.6 gate plan ran in full: eight beads, eight gated modules in
`analysis/model/su2/` (PRs #45–#52), 31 rows, 0 failures. Condensed record:
T2 `su2-gate-lane-results-2026-09-01` [23914]. What each gate returned:

- **G0** (`lift.py`): the lift instrument calibrated — trivial loops +q, bare
  2π −q, 4π +q; coarse or unclosed paths refused by the instrument itself.
- **G1** (`plate_holonomy.py`): `q(300°) = −q(−60°)` for **432/432 plates**
  on all three patches, cells and voids alike — §8.3's near-corollary is now
  a measurement. Bonus: the spin law `R ≡ rot(u_f, σ_f(a+60°))` holds to
  2e-15 through the whole cycle.
- **G4** (`inside_out_cover.py`): `X(a+180°) = −X(a)` exact, the conjectured
  const **zero**, inversion centre the lattice origin. Deck map (b) pinned.
  New: sheet-touching in ambient space is patch geometry — a site-symmetric
  patch's two sheets are ambiently indistinguishable at *every* angle.
- **G2** (`joint_twist.py`, subsuming G3): the relative-rotation lift across
  **both** co-located vertex joints is **+q** for every plate pair — the ℤ₂
  cancels pairwise (sign(rel) = product of absolute signs, unjoined controls
  included). Stronger than §8.2: there is no ℤ₂ in any pairwise relative
  history *to* store; even a tether-capable joint would return untwisted.
  The −q exists only against an external frame.
- **G7** (`physical_arc.py`): the physical sixty degrees never exercises the
  holonomy — crossing-free exactly on [−60°, 0°], walled within one degree
  beyond both ends; oscillation of any amplitude lifts +q, only monotone
  traversal reaches −q. §8.5's strongest objection is a gated fact.
- **τ prerequisites** (`screw_prerequisites.py`): single applications of τ
  leave the coherent family (the parity flip demands +60° and −60° at once);
  τ⁶ is a pure lattice translation exactly; τ³ = translate ∘ per-body central
  inversion; the action is free and properly discontinuous on the physical
  arc only. The Bieberbach naming of §3a stays withdrawn.
- **G5** (`doubled_block.py`): the §3c double **built** — every half-weld
  closed onto the twin sheet, all cells at weld-degree six. It does *not*
  lose the seven zero modes: its spectrum is exactly the free block's
  (Neumann sector) plus a pinned-boundary Dirichlet sector with none. In
  Kolmogorov distance to the Bloch density of states it runs at less than
  half the free block's distance at every size and at the exact periodic
  reference's level from side 3.
- **G6** (`screw_sector.py`): resolved in the **negative** — odd multiples of
  (1,1,1) put zero solid sites on solid sites (under DECISION 21 the screw
  has nothing to act on); no isometry pairs the cell and void bodies at any
  generic fold (0/48); the glide pairing at a = −30° (12+12 of 48) belongs
  to the retired both-kinds covering.

**Verdict of the lane.** The SU(2) intuition is correct kinematics (G1),
physically inert (G7), with no mechanical carrier even under joint redesign
(G2) and no spectral sector (G6). The ℤ₂ is a property of the fold circle's
analytic continuation, visible only against an external frame. The original
question of §1 is answered by G5's race: free boundaries are dominated; the
double and the torus are spectrally equivalent closures at the sizes tried.

**OWNER DECISION 22 (2026-09-01): the closure of record is the DOUBLE** —
§3c's two-sided-plate closure, chosen over the torus ("double for sure. that
is the most elegant"). It keeps all seven zero modes where the torus keeps
four, it is a real assembly rather than operator surgery, and it realises the
branched-cover picture this note began from as a buildable constraint system.
Accepted cost: twice the cells of the box it closes. Constructor:
`analysis/model/su2/doubled_block.py::double(side)`. Record: T2
`DECISION-22-2026-09-01-closure-is-the-double` [23932].
