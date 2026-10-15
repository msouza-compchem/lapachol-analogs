#!/usr/bin/env python3
"""Etapa 3 do funil: descritores eletronicos (GFN2-xTB, vertical IP/EA) dos 17 compostos.

Replica o protocolo de quinone-ml/02_xtb_descriptors.py (modo 'sp'):
  3D com RDKit ETKDGv3 (seed 42, useSmallRingTorsions) + MMFF (500 passos), SEM otimizar no xTB,
  xtb --gfn 2 --alpb water --vipea ; mu = -(IP+EA)/2, eta = (IP-EA)/2, omega = mu^2/(2 eta).

Por que direto no xTB e nao pelo modelo: para 17 moleculas o xTB (10-20 s cada) e a referencia
com a qual o GNN foi treinado; o modelo so compensa em bibliotecas de milhares de compostos.

Estados: neutro (todos) e anion (os 10 acidos; xtb --chrg -1). O dataset do quinone-ml so tem
especies NEUTRAS; os numeros do anion sao EXPLORATORIOS (solvatacao implicita de anions e menos
confiavel e mu/eta/omega do anion nao sao comparaveis com os dos neutros).

Uso: python 04_descritores_xtb.py --xtb /caminho/xtb [--quinone-ml /caminho/quinone-ml]
"""
import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
HERE = Path(__file__).resolve().parent
ACID = Chem.MolFromSmarts("[OX2H1][C;!a]=[C;!a][C;!a]=O")


def desprotonar(mol):
    hits = mol.GetSubstructMatches(ACID)
    if len(hits) != 1:
        raise ValueError(f"esperava 1 OH acido, achei {len(hits)}")
    rw = Chem.RWMol(mol)
    rw.GetAtomWithIdx(hits[0][0]).SetFormalCharge(-1)
    Chem.SanitizeMol(rw)
    return rw.GetMol()


def to_xyz(mol, path):
    """Igual a quinone-ml: ETKDGv3, seed 42, MMFF 500 passos."""
    mol = Chem.AddHs(mol)
    p = AllChem.ETKDGv3()
    p.randomSeed = 42
    p.useSmallRingTorsions = True
    if AllChem.EmbedMolecule(mol, p) != 0:
        p.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, p) != 0:
            return False
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        pass
    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), "mol"]
    for a in mol.GetAtoms():
        q = conf.GetAtomPosition(a.GetIdx())
        lines.append(f"{a.GetSymbol()} {q.x:.6f} {q.y:.6f} {q.z:.6f}")
    Path(path).write_text("\n".join(lines) + "\n")
    return True


def parse_xtb(text):
    out = {}
    m = re.search(r"delta SCC IP \(eV\):\s+(-?[\d.]+)", text)
    if m:
        out["ip"] = float(m.group(1))
    m = re.search(r"delta SCC EA \(eV\):\s+(-?[\d.]+)", text)
    if m:
        out["ea"] = float(m.group(1))
    m = re.search(r"TOTAL ENERGY\s+(-?[\d.]+)\s+Eh", text)
    if m:
        out["e_total"] = float(m.group(1))
    return out


def run_one(job):
    lig_id, nome, estado, smiles, xtb_bin, xtbhome = job
    rec = dict(id=lig_id, nome=nome, estado=estado, smiles=smiles, status="ok")
    mol = Chem.MolFromSmiles(smiles)
    charge = Chem.GetFormalCharge(mol)
    wd = Path(tempfile.mkdtemp(prefix="xtb_"))
    try:
        if not to_xyz(mol, wd / "mol.xyz"):
            rec["status"] = "embed_failed"
            return rec
        env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
        if xtbhome:
            env["XTBHOME"] = xtbhome
        cmd = [xtb_bin, "mol.xyz", "--gfn", "2", "--alpb", "water", "--vipea"]
        if charge:
            cmd += ["--chrg", str(charge)]
        proc = subprocess.run(cmd, cwd=wd, capture_output=True, text=True, timeout=600, env=env)
        parsed = parse_xtb(proc.stdout)
        if "ip" not in parsed or "ea" not in parsed:
            rec["status"] = "scf_failed" if proc.returncode else "parse_failed"
            return rec
        rec.update(parsed)
        ip, ea = parsed["ip"], parsed["ea"]
        mu, eta = -(ip + ea) / 2.0, (ip - ea) / 2.0
        rec.update(mu=round(mu, 4), eta=round(eta, 4),
                   omega=round(mu * mu / (2 * eta), 4) if eta > 0.05 else None)
        return rec
    except subprocess.TimeoutExpired:
        rec["status"] = "timeout"
        return rec
    finally:
        shutil.rmtree(wd, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xtb", default=shutil.which("xtb") or "xtb")
    ap.add_argument("--quinone-ml", default="", help="pasta do repo quinone-ml, para cruzar com o dataset")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--only", default="", help="ids (ex.: L00,L07)")
    a = ap.parse_args()
    xtbhome = str(Path(a.xtb).resolve().parent.parent) if os.path.exists(a.xtb) else ""

    jobs = []
    with open(HERE / "analogos_lapachol.csv", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            if a.only and r["id"] not in a.only.split(","):
                continue
            mol = Chem.MolFromSmiles(r["smiles"])
            jobs.append((r["id"], r["nome"], "neutro", Chem.MolToSmiles(mol), a.xtb, xtbhome))
            if r["estado_ionizacao_pH7_4"].startswith("anion"):
                jobs.append((r["id"], r["nome"], "anion", Chem.MolToSmiles(desprotonar(mol)), a.xtb, xtbhome))

    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        rows = list(pool.map(run_one, jobs))

    # cruza com o dataset do quinone-ml (so especies neutras)
    ds = {}
    if a.quinone_ml:
        import pandas as pd
        d = pd.read_csv(Path(a.quinone_ml) / "data" / "descriptors_final.csv")
        for r in d.itertuples():
            m = Chem.MolFromSmiles(r.smiles)
            if m is not None:
                ds[Chem.MolToSmiles(m)] = r
    for r in rows:
        hit = ds.get(r["smiles"]) if r["estado"] == "neutro" else None
        r["no_dataset_quinone_ml"] = "sim" if hit is not None else "nao"
        r["ea_dataset"] = round(float(hit.ea), 4) if hit is not None else ""
        r["ip_dataset"] = round(float(hit.ip), 4) if hit is not None else ""

    fields = ["id", "nome", "estado", "smiles", "status", "ip", "ea", "mu", "eta", "omega",
              "e_total", "no_dataset_quinone_ml", "ip_dataset", "ea_dataset"]
    out = HERE / "descritores_xtb.csv"
    with open(out, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    print(f"{'id':4} {'estado':7} {'IP':>6} {'EA':>6} {'mu':>7} {'eta':>6} {'omega':>6}  dataset(EA)  {'nome'}")
    for r in rows:
        f = lambda k: f"{r[k]:6.3f}" if isinstance(r.get(k), float) else f"{'-':>6}"
        print(f"{r['id']:4} {r['estado']:7} {f('ip')} {f('ea')} {r.get('mu', float('nan')):7.3f} "
              f"{f('eta')} {f('omega')}  {str(r['ea_dataset']):>11}  {r['nome']}"
              + ("" if r["status"] == "ok" else f"  [{r['status']}]")
              if r["status"] == "ok" else f"{r['id']:4} {r['estado']:7} FALHOU: {r['status']}")
    print(f"\n{len(rows)} calculos -> {out}")


if __name__ == "__main__":
    sys.exit(main())
