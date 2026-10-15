#!/usr/bin/env python3
"""Etapa 2 (analise): metricas de pose, nao so escore.

Para cada ligante encaixado (poses/<id>_out.pdbqt) calcula, por modo do Vina:
  escore            : escore do Vina (ordena poses; NAO e energia livre)
  core_rmsd         : RMSD (A) do esqueleto bicíclico de 10 carbonos contra a pose de referencia do
                      lapachol (reference/best_pose.pdb, a que ficou empilhada na MD). Sem alinhamento;
                      minimo sobre as simetrias do esqueleto.
  d_centroide       : distancia (A) entre o centroide do anel quinonico (anel do esqueleto com mais O
                      terminais) e o centroide do anel de 14 atomos da isoaloxazina do FAD (cadeia A).
  angulo_planos     : angulo agudo (graus) entre os planos desses dois aneis (0 = paralelos).
  d_N5              : menor distancia (A) entre o N5 do FAD e os 6 atomos do anel quinonico.
  d_min_FAD         : menor distancia (A) entre qualquer atomo pesado do ligante e a isoaloxazina
                      (incl. O2, O4 e metilas).
Resumo por ligante: modo 1 + quantos modos tem core_rmsd <= 2.0 A ("tipo-referencia") e o melhor
escore entre eles. O limiar de 2.0 A e o mesmo criterio de aceitacao do redocking no README;
aqui ele mede semelhanca com a pose de referencia, nao correcao.

d_N5 e d_centroide sao descritores geometricos, nao previsao de reatividade.

Tambem confere o redocking da duroquinona (controle) contra reference/dqn_clean.pdb.
"""
import csv
import math
import os
import sys

import numpy as np
from rdkit import Chem

HERE = os.path.dirname(os.path.abspath(__file__))
REC = os.path.join(HERE, "receptor", "receptor.pdbqt")
REF_LAP = os.path.join(HERE, "reference", "best_pose.pdb")
REF_DQN = os.path.join(HERE, "reference", "dqn_clean.pdb")
LIMIAR_RMSD = 2.0

ISO_RING = ["N1", "C2", "N3", "C4", "C4X", "N5", "C5X", "C6", "C7", "C8", "C9", "C9A", "N10", "C10"]
ISO_ALL = ISO_RING + ["O2", "O4", "C7M", "C8M"]
RADII = {"C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "S": 1.05, "Cl": 1.02, "Br": 1.20, "P": 1.07}
ADTYPE = {"C": "C", "A": "C", "N": "N", "NA": "N", "OA": "O", "O": "O", "SA": "S", "S": "S",
          "F": "F", "Cl": "Cl", "CL": "Cl", "Br": "Br", "BR": "Br", "P": "P", "HD": "H", "H": "H"}

# esqueleto naftaleno (10 C, 2 aneis fundidos); ligacoes '~' = qualquer ordem
CORE = Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#6]2~[#6](~[#6]~1)~[#6]~[#6]~[#6]~[#6]~2")
RING_A, RING_B = (0, 1, 2, 3, 4, 5), (3, 4, 6, 7, 8, 9)


# ---------------------------------------------------------------- leitura
def read_pdbqt_models(path):
    models, cur, score = [], None, None
    for ln in open(path):
        if ln.startswith("MODEL"):
            cur, score = [], None
        elif ln.startswith("REMARK VINA RESULT"):
            score = float(ln.split()[3])
        elif ln.startswith(("ATOM", "HETATM")) and cur is not None:
            t = ln.split()[-1]
            el = ADTYPE.get(t, t[:1])
            if el == "H":
                continue
            cur.append((ln[12:16].strip(), el,
                        np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])))
        elif ln.startswith("ENDMDL") and cur is not None:
            models.append(dict(score=score, atoms=cur))
            cur = None
    return models


def read_pdb_heavy(path):
    atoms = []
    for ln in open(path):
        if ln.startswith("ENDMDL"):
            break
        if ln.startswith(("ATOM", "HETATM")):
            el = ln[76:78].strip().capitalize() or "".join(c for c in ln[12:16] if c.isalpha())[:1].upper()
            if el == "H":
                continue
            atoms.append((ln[12:16].strip(), el,
                          np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])))
    return atoms


def read_fad(path):
    fad = {}
    for ln in open(path):
        if ln.startswith(("ATOM", "HETATM")) and ln[17:20].strip() == "FAD" and ln[21] == "A":
            fad[ln[12:16].strip()] = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    missing = [n for n in ISO_ALL if n not in fad]
    if missing:
        raise SystemExit(f"atomos do FAD ausentes no receptor: {missing}")
    return fad


# ---------------------------------------------------------------- geometria
def skeleton(atoms):
    m = Chem.RWMol()
    for _, el, _ in atoms:
        m.AddAtom(Chem.Atom(el))
    xyz = np.array([a[2] for a in atoms])
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            lim = RADII[atoms[i][1]] + RADII[atoms[j][1]] + 0.45
            if np.linalg.norm(xyz[i] - xyz[j]) < lim:
                m.AddBond(i, j, Chem.BondType.SINGLE)
    m = m.GetMol()
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    return m


def quinone_ring(skel, match):
    def n_term_o(ring):
        n = 0
        for idx in ring:
            for nb in skel.GetAtomWithIdx(match[idx]).GetNeighbors():
                if nb.GetSymbol() == "O" and nb.GetDegree() == 1:
                    n += 1
        return n
    a, b = n_term_o(RING_A), n_term_o(RING_B)
    ring = RING_A if a >= b else RING_B
    return [match[i] for i in ring], max(a, b)


def plane(P):
    c = P.mean(axis=0)
    return c, np.linalg.svd(P - c)[2][2]


def pose_metrics(atoms, fad, ref_core):
    skel = skeleton(atoms)
    matches = skel.GetSubstructMatches(CORE, uniquify=False, maxMatches=2000)
    nan = dict(core_rmsd=np.nan, d_centroide=np.nan, angulo_planos=np.nan,
               d_N5=np.nan, d_min_FAD=np.nan)
    if not matches:
        return nan
    xyz = np.array([a[2] for a in atoms])
    qring, n_o = quinone_ring(skel, matches[0])
    if n_o < 2:
        return nan
    Q = xyz[qring]
    cq, nq = plane(Q)
    iso = np.array([fad[n] for n in ISO_RING])
    cf, nf = plane(iso)
    iso_all = np.array([fad[n] for n in ISO_ALL])
    out = dict(
        d_centroide=float(np.linalg.norm(cq - cf)),
        angulo_planos=math.degrees(math.acos(min(1.0, abs(float(np.dot(nq, nf)))))),
        d_N5=float(np.linalg.norm(Q - fad["N5"], axis=1).min()),
        d_min_FAD=float(np.linalg.norm(xyz[:, None, :] - iso_all[None, :, :], axis=2).min()),
    )
    out["core_rmsd"] = (min(math.sqrt(float(np.mean(np.sum((xyz[list(m)] - ref_core) ** 2, axis=1))))
                            for m in matches) if ref_core is not None else np.nan)
    return out


def sym_rmsd(pose_atoms, ref_atoms):
    """RMSD sem alinhamento, atomos pesados, minimo sobre automorfismos (esqueleto so com ligacoes simples)."""
    ref_sk, pose_sk = skeleton(ref_atoms), skeleton(pose_atoms)
    ref_xyz = np.array([a[2] for a in ref_atoms])
    pose_xyz = np.array([a[2] for a in pose_atoms])
    ms = pose_sk.GetSubstructMatches(ref_sk, uniquify=False, maxMatches=5000)
    if not ms:
        return np.nan
    return min(math.sqrt(float(np.mean(np.sum((pose_xyz[list(m)] - ref_xyz) ** 2, axis=1)))) for m in ms)


def fmt(x, nd=2):
    return "" if x is None or (isinstance(x, float) and math.isnan(x)) else round(x, nd)


# ---------------------------------------------------------------- principal
def main():
    fad = read_fad(REC)
    lap_ref = read_pdb_heavy(REF_LAP)
    skel_ref = skeleton(lap_ref)
    m0 = skel_ref.GetSubstructMatches(CORE, uniquify=False, maxMatches=2000)
    if not m0:
        raise SystemExit("esqueleto bicíclico nao encontrado na pose de referencia")
    ref_xyz = np.array([a[2] for a in lap_ref])
    ref_core = ref_xyz[list(m0[0])]
    ref = pose_metrics(lap_ref, fad, ref_core)
    print("Referencia (best_pose.pdb, lapachol): " +
          ", ".join(f"{k}={v:.2f}" for k, v in ref.items()))

    with open(os.path.join(HERE, "ligands", "manifest.csv"), encoding="utf-8") as fh:
        manifest = list(csv.DictReader(fh))
    props = {}
    with open(os.path.join(HERE, "analogos_lapachol.csv"), encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            props[r["id"]] = r

    # controle de redocking
    dq = os.path.join(HERE, "poses", "CTRL_DQN_neutro_out.pdbqt")
    if os.path.exists(dq):
        ref_dq = read_pdb_heavy(REF_DQN)
        mods = read_pdbqt_models(dq)
        rm = [sym_rmsd(m["atoms"], ref_dq) for m in mods]
        k = int(np.nanargmin(rm)) + 1
        print(f"\nRedocking duroquinona (controle): modo 1 escore {mods[0]['score']:.2f}, RMSD {rm[0]:.2f} A; "
              f"menor RMSD entre {len(mods)} modos = {np.nanmin(rm):.2f} A (modo {k}, escore {mods[k - 1]['score']:.2f}). "
              f"Criterio < 2.0 A no modo 1: {'sim' if rm[0] < 2.0 else 'NAO'}; em algum modo: "
              f"{'sim' if np.nanmin(rm) < 2.0 else 'NAO'}. (README, preparo com obabel: 1.11 A no modo 1.)")

    rows = []
    for m in manifest:
        if m["id"] == "CTRL_DQN":
            continue
        path = os.path.join(HERE, "poses", f"{m['ligand_id']}_out.pdbqt")
        if not os.path.exists(path):
            continue
        models = read_pdbqt_models(path)
        per = [pose_metrics(md["atoms"], fad, ref_core) for md in models]
        scores = [md["score"] for md in models]
        like = [i for i, p in enumerate(per) if not math.isnan(p["core_rmsd"]) and p["core_rmsd"] <= LIMIAR_RMSD]
        # modos vem ordenados por escore (mais negativo = melhor): o primeiro 'tipo-ref' e o melhor deles
        best_like = like[0] if like else None
        p1 = per[0]
        p = props[m["id"]]
        rows.append(dict(
            ligand_id=m["ligand_id"], id=m["id"], nome=m["nome"], estado=m["estado"],
            logP_Crippen=p["logP_Crippen"], logD7_4_aprox=p["logD7_4_aprox"],
            escore_modo1=fmt(scores[0]),
            core_rmsd_modo1=fmt(p1["core_rmsd"]), d_centroide_modo1=fmt(p1["d_centroide"]),
            angulo_planos_modo1=fmt(p1["angulo_planos"], 1), d_N5_modo1=fmt(p1["d_N5"]),
            d_min_FAD_modo1=fmt(p1["d_min_FAD"]),
            n_modos=len(models), n_modos_tipo_ref=len(like),
            melhor_escore_tipo_ref=fmt(scores[best_like]) if best_like is not None else "",
            modo_tipo_ref=best_like + 1 if best_like is not None else "",
        ))

    base = {(r["estado"]): float(r["escore_modo1"]) for r in rows if r["id"] == "L00"}
    for r in rows:
        r["delta_escore_vs_L00_mesmo_estado"] = (
            fmt(float(r["escore_modo1"]) - base[r["estado"]]) if r["estado"] in base else "")

    if not rows:
        print("nenhuma pose de biblioteca encontrada em poses/")
        return
    out = os.path.join(HERE, "resultados_docking.csv")
    with open(out, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"\n{'ligante':14} {'escore':>6} {'rmsd':>5} {'d_cent':>6} {'ang':>5} {'dN5':>5} {'dmin':>5} "
          f"{'tipo-ref':>8} {'melhor':>7}")
    for r in rows:
        print(f"{r['ligand_id']:14} {r['escore_modo1']:>6} {r['core_rmsd_modo1']:>5} "
              f"{r['d_centroide_modo1']:>6} {r['angulo_planos_modo1']:>5} {r['d_N5_modo1']:>5} "
              f"{r['d_min_FAD_modo1']:>5} {str(r['n_modos_tipo_ref']) + '/' + str(r['n_modos']):>8} "
              f"{r['melhor_escore_tipo_ref']:>7}")
    print(f"\n{len(rows)} ligantes -> {out}")


if __name__ == "__main__":
    sys.exit(main())
