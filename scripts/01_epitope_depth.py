#!/usr/bin/env python3
"""
Split the class A hotspots by burial depth, and flag the lipid-facing surface.

Why this is needed
------------------
Two failure modes sit between the computed hotspot list and a usable
RFdiffusion input, and both are specific to a GPCR target:

1. DEPTH. The 27 class A hotspots include residues lining the deep 7TM core
   (E387 and friends). A folded 60-90 aa scaffold physically cannot reach them;
   only an extended N-terminal "message module" can. If those residues are
   passed to RFdiffusion as ppi.hotspot_res it will try to bury a scaffold in a
   pocket it cannot enter and produce garbage. The scaffold should be aimed at
   the extracellular vestibule RIM; the deep residues are the message module's
   job (brief section 5).

2. LIPID. RFdiffusion has no concept of a membrane. The TMD's lipid-facing
   surface is a big hydrophobic convex patch -- exactly what the model likes to
   design against. Left unconstrained it will happily generate binders stuck to
   the side of the transmembrane bundle, which cannot bind in a detergent
   micelle assay. We compute the membrane normal so these residues can be
   identified and excluded.

Method
------
* Membrane normal: first principal axis of the TMD Ca coordinates. For a 7TM
  bundle this is the membrane normal to within a few degrees.
* Extracellular direction: the sign that puts the ECD centroid on the positive
  side. Derived from the structure, not assumed.
* Depth of a residue: projection onto that axis, measured from the most
  extracellular TMD Ca. Larger = deeper into the membrane.
* Burial: count of receptor heavy atoms within 10 A of the residue's side-chain
  centroid. High burial + accessible depth = pocket rim; high burial + deep =
  core.
* Lipid-facing: side-chain centroid points away from the bundle axis (radially
  outward) AND sits inside the membrane slab.

Outputs target/epitope_classes.json and prints RFdiffusion-ready hotspot strings.
"""

import json
import os
import sys
from collections import defaultdict

import numpy as np

TMD_PDB = "target/glp1r_TMD.pdb"
ECD_PDB = "target/glp1r_ECD.pdb"
HOTSPOTS = "target/hotspots.json"

MESSAGE_MAX_RES = 15     # GLP-1 res 15 marks the top of the message module
RIM_DEPTH_MAX = 12.0     # retained for reference; superseded by the GLP-1 ruler
BACKBONE = {"N", "CA", "C", "O"}


def read_res(pdb):
    res = defaultdict(lambda: {"name": None, "atoms": []})
    for l in open(pdb):
        if l[:4] != "ATOM":
            continue
        n = int(l[22:26])
        res[n]["name"] = l[17:20].strip()
        res[n]["atoms"].append((l[12:16].strip(),
                                np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return dict(res)


def sc_centroid(r):
    sc = [xyz for an, xyz in r["atoms"] if an not in BACKBONE]
    if not sc:
        sc = [xyz for an, xyz in r["atoms"] if an == "CA"]
    return np.mean(sc, axis=0)


def main():
    for f in (TMD_PDB, ECD_PDB, HOTSPOTS):
        if not os.path.isfile(f):
            sys.exit(f"missing {f}; run scripts/00_prep_target.py first")

    tmd = read_res(TMD_PDB)
    ecd = read_res(ECD_PDB)
    hs = json.load(open(HOTSPOTS))
    classA = hs["targets"]["A_orthosteric_vestibule"]["hotspots"]

    ca = np.array([next(xyz for an, xyz in tmd[n]["atoms"] if an == "CA")
                   for n in sorted(tmd)])
    ca_nums = sorted(tmd)
    centre = ca.mean(axis=0)

    # membrane normal = first principal axis of the TM bundle
    u, s, vt = np.linalg.svd(ca - centre)
    axis = vt[0] / np.linalg.norm(vt[0])
    print(f"TM bundle principal axis explains "
          f"{s[0]**2 / (s**2).sum() * 100:.0f}% of Ca variance")

    # orient: extracellular = side the ECD sits on
    ecd_centroid = np.mean([sc_centroid(r) for r in ecd.values()], axis=0)
    if np.dot(ecd_centroid - centre, axis) < 0:
        axis = -axis
    print(f"extracellular direction set from ECD centroid "
          f"({np.dot(ecd_centroid - centre, axis):.1f} A on the + side)")

    proj = {n: float(np.dot(c - centre, axis)) for n, c in zip(ca_nums, ca)}
    top = max(proj.values())
    depth = {n: top - p for n, p in proj.items()}   # 0 = most extracellular

    # radial direction for lipid-facing test
    all_sc = {n: sc_centroid(tmd[n]) for n in tmd}
    radial = {}
    for n, c in all_sc.items():
        v = c - centre
        v_perp = v - np.dot(v, axis) * axis
        r = np.linalg.norm(v_perp)
        ca_n = next(xyz for an, xyz in tmd[n]["atoms"] if an == "CA")
        v2 = ca_n - centre
        v2_perp = v2 - np.dot(v2, axis) * axis
        # does the side chain point outward relative to its own backbone?
        radial[n] = (r, float(np.dot(v_perp - v2_perp, v2_perp /
                                     (np.linalg.norm(v2_perp) + 1e-9))))

    coords_all = np.vstack([xyz for r in tmd.values() for _, xyz in r["atoms"]])
    burial = {}
    for n in tmd:
        c = all_sc[n]
        burial[n] = int((np.linalg.norm(coords_all - c, axis=1) < 10.0).sum())

    med_r = np.median([radial[n][0] for n in tmd])

    # ---- use GLP-1 itself as the depth ruler -------------------------------
    # Depth measured from the most extracellular TMD Ca is misleading: that Ca
    # is an extracellular loop tip, not the vestibule mouth. The native ligand
    # is the honest ruler -- whatever depth GLP-1 reaches, a ligand can reach.
    glp1 = read_res("target/6X18_P.pdb") if os.path.isfile("target/6X18_P.pdb") else None
    if glp1 is None:
        glp1 = {}
        for l in open("target/6X18.pdb"):
            if l[:4] == "ATOM" and l[21] == "P":
                n = int(l[22:26])
                glp1.setdefault(n, {"name": l[17:20].strip(), "atoms": []})
                glp1[n]["atoms"].append(
                    (l[12:16].strip(),
                     np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))

    glp1_depth = {}
    for n, r in glp1.items():
        c = next((xyz for an, xyz in r["atoms"] if an == "CA"), None)
        if c is not None:
            glp1_depth[n] = top - float(np.dot(c - centre, axis))
    print("\nGLP-1 depth along the same axis (the ruler):")
    for n in sorted(glp1_depth):
        if n <= 20 or n % 4 == 0:
            print(f"  GLP-1 {n:>3} {glp1[n]['name']:>4}  depth {glp1_depth[n]:>6.1f} A")

    # The vestibule mouth is where GLP-1 crosses out of the bundle. Take it as
    # the depth of the shallowest message residue (res 15).
    mouth = glp1_depth[min(glp1_depth, key=lambda k: abs(k - MESSAGE_MAX_RES))]
    print(f"\nvestibule mouth taken as depth of GLP-1 {MESSAGE_MAX_RES}: {mouth:.1f} A")
    print(f"deepest GLP-1 residue (His7): {glp1_depth[min(glp1_depth)]:.1f} A")

    pc = hs["per_glp1_residue_contacts"]
    shallowest_contact = {}
    for p, partners in pc.items():
        for rnum in partners:
            rn = int(rnum)
            shallowest_contact[rn] = min(shallowest_contact.get(rn, 1e9),
                                         glp1_depth.get(int(p), 1e9))

    rows = []
    for n in classA:
        if n not in tmd:
            continue
        d = depth[n]
        r, outward = radial[n]
        reach = shallowest_contact.get(n, 1e9)   # shallowest GLP-1 depth touching it
        lipid = (outward > 0.5) and (r > med_r) and (d > mouth)
        if lipid:
            cls = "LIPID"
        elif reach <= mouth:
            cls = "rim"
        else:
            cls = "deep"
        rows.append((n, tmd[n]["name"], d, burial[n], r, cls))

    rows.sort(key=lambda x: x[2])
    print(f"\n{'res':>5} {'aa':>4} {'depth':>7} {'burial':>7} {'radius':>7}  class")
    for n, aa, d, b, r, cls in rows:
        print(f"{n:>5} {aa:>4} {d:>7.1f} {b:>7} {r:>7.1f}  {cls}")

    rim = [n for n, _, _, _, _, c in rows if c == "rim"]
    deep = [n for n, _, _, _, _, c in rows if c == "deep"]
    lip = [n for n, _, _, _, _, c in rows if c == "LIPID"]

    print(f"\nrim  (scaffold epitope, depth <= {RIM_DEPTH_MAX:.0f} A): {rim}")
    print(f"deep (message-module territory)              : {deep}")
    print(f"lipid-facing (EXCLUDE)                       : {lip}")

    # all lipid-facing TMD residues, for later filtering of design interfaces
    lipid_all = sorted(n for n in tmd
                       if radial[n][1] > 0.5 and radial[n][0] > med_r
                       and depth[n] > RIM_DEPTH_MAX)
    print(f"\n{len(lipid_all)} of {len(tmd)} TMD residues are lipid-facing overall")

    out = {
        "method": "PCA membrane normal on TMD Ca; depth from most extracellular Ca; "
                  "burial = receptor heavy atoms within 10 A of side-chain centroid",
        "depth_ruler": "GLP-1 Ca depth along the TM axis; vestibule mouth = depth of GLP-1 res 15",
        "classA_rim_hotspots": rim,
        "classA_deep_residues": deep,
        "classA_lipid_facing": lip,
        "tmd_lipid_facing_all": lipid_all,
        "depth_A": {str(n): round(depth[n], 2) for n in sorted(tmd)},
    }
    json.dump(out, open("target/epitope_classes.json", "w"), indent=2)
    print("\nwrote target/epitope_classes.json")

    print("\nRFdiffusion hotspot string (class A scaffold):")
    print("  ppi.hotspot_res=[" + ",".join(f"A{n}" for n in rim) + "]")
    print("RFdiffusion hotspot string (class C, ECD groove):")
    print("  ppi.hotspot_res=["
          + ",".join(f"A{n}" for n in hs["targets"]["C_ecd_groove"]["hotspots"]) + "]")


if __name__ == "__main__":
    main()
