"""Step G: does the 60-120 interference survive removing the faces?

Raised by the project owner (2026-08-11): "this is a model, not physical
reality. Maybe the triangles are patterns."

That is a serious objection to every collision number recorded so far, because
`verify_critic.interpenetrations` counts EDGE-THROUGH-FACE events, which assumes
the triangles are filled. Fuller's own jitterbug has no faces at all -- it is 24
struts cohered at 12 flexible joints (Synergetics 460.00) -- and an edge crossing
the open middle of a strut triangle would hit nothing. If the boundary at a=60
were an artefact of filling the faces, it would be an artefact of an assumption
nobody had stated, which is this project's recurring failure mode.

So: drop the faces entirely and measure strut against strut.

RESULT: the boundary survives. The struts themselves cross on exactly the same
interval. Minimum clearance falls 1.224744871 (a=0) -> 0.040296585 (a=59) ->
EXACTLY 0 throughout 60-120 -> 0.394346678 (a=130) -> 1.224744871 (a=180).
Two independent tests -- edge-through-filled-face, and edge-to-edge distance with
no faces at all -- agree on (60,120) U (240,300). The interference is a property
of the EDGE geometry.

What that leaves open is one step further in, and it is a modelling decision
rather than a measurement: whether interference DISQUALIFIES a configuration.
Tracked as bead inviscid-qvf.8.
"""
import numpy as np

from analysis.model.jitterbug import corners, L_EDGE
from analysis.model.linkage_variety import PAIRS


def segment_distance(p1, q1, p2, q2):
    """Shortest distance between two 3D segments (clamped parametric solve)."""
    d1, d2, r = q1 - p1, q2 - p2, p1 - p2
    a, e = d1 @ d1, d2 @ d2
    if a < 1e-14 or e < 1e-14:
        return np.linalg.norm(r)
    b, c, f = d1 @ d2, d1 @ r, d2 @ r
    den = a * e - b * b
    s = np.clip((b * f - c * e) / den, 0.0, 1.0) if den > 1e-12 else 0.0
    t = np.clip((b * s + f) / e, 0.0, 1.0)
    s = np.clip((b * t - c) / a, 0.0, 1.0)
    return np.linalg.norm((p1 + d1 * s) - (p2 + d2 * t))


def struts(a):
    """The 24 struts as ((face,corner),(face,corner), p, q). No faces involved."""
    X = corners(a)
    return [((i, j), (i, (j + 1) % 3), X[i, j], X[i, (j + 1) % 3])
            for i in range(8) for j in range(3)]


def _hinged(ends_a, ends_b):
    """Do these two struts legitimately meet at a shared hinge vertex?"""
    if ends_a & ends_b:
        return True
    for (pa, pb) in PAIRS:
        if (pa in ends_a and pb in ends_b) or (pb in ends_a and pa in ends_b):
            return True
    return False


def min_strut_gap(a):
    S = struts(a)
    best = np.inf
    for m in range(len(S)):
        for n in range(m + 1, len(S)):
            ea, eb = {S[m][0], S[m][1]}, {S[n][0], S[n][1]}
            if _hinged(ea, eb):
                continue
            best = min(best, segment_distance(S[m][2], S[m][3], S[n][2], S[n][3]))
    return best


if __name__ == "__main__":
    print(f"strut length = {L_EDGE:.9f}   (24 struts, 12 flexible joints -- no faces)")
    print("\n     a      min strut-strut gap    state")
    for a in (0, 15, 30, 45, 55, 59, 59.9, 60, 61, 70, 80, 90,
              100, 110, 119, 120, 121, 130, 150, 180, 240, 270, 300, 330):
        g = min_strut_gap(a)
        state = "STRUTS CROSS" if g < 1e-9 else ""
        print(f"  {a:7.2f}    {g:.9f}          {state}")

    print("\ncross-check against the face-based count.")
    print("The two tests agree everywhere EXCEPT the merge angles 60/120/240/300,")
    print("and that exception is a real distinction rather than a discrepancy: there")
    print("the struts COINCIDE in pairs (the octahedron's doubled edges, multiplicity")
    print("2) instead of crossing, so the gap is 0 while nothing pierces an interior.")
    print("Touching is not penetrating -- which is exactly why a=60 is the BOUNDARY.\n")
    from analysis.model.interpenetration_check import interpenetrations
    MERGE = {60.0, 120.0, 240.0, 300.0}
    for a in (30, 59, 60, 61, 90, 119, 120, 121, 180, 250, 310):
        g, n = min_strut_gap(a), interpenetrations(corners(a))
        if float(a) in MERGE:
            ok = g < 1e-9 and n == 0
            note = "boundary: struts coincide, nothing penetrates"
        else:
            ok = (g < 1e-9) == (n > 0)
            note = "agree"
        print(f"  a={a:6.1f}  strut gap {g:.9f}  face-pierce count {n:3d}   "
              f"{note if ok else 'UNEXPECTED'}")
