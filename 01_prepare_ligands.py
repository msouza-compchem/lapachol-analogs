#!/usr/bin/env python3
"""Etapa 2 (preparo): ligantes em PDBQT com estado de protonacao EXPLICITO.

Entrada : analogos_lapachol.csv (saida da etapa 1)
Saida   : ligands/<id>_<estado>.pdbqt  e  ligands/manifest.csv

Estados:
  neutro : como escrito no SMILES.
  anion  : 2-OH enolico (vinilogo-acido da naftoquinona) desprotonado, carga -1.
           Gerado para todo composto cuja classe tem pKa estimado < ~6
           (coluna estado_ionizacao_pH7_4 do CSV comeca com 'anion').

Por que explicito: o protocolo original (obabel -p 7.4 a partir de .xyz) nao desprotona
esse tipo de enol; o lapachol foi encaixado e simulado na forma NEUTRA.
O Vina nao tem termo eletrostatico na funcao de escore: o anion difere do neutro
apenas pela tipagem de atomos (doador/aceptor). Tratar a diferenca como fraca.

Geometria 3D: RDKit ETKDG + MMFF (UFF se MMFF nao tipar). O Vina trata o ligante como
flexivel, entao o confomero inicial importa pouco.
"""
import csv
import os
import sys

from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "analogos_lapachol.csv")
OUTD = os.path.join(HERE, "ligands")
os.makedirs(OUTD, exist_ok=True)

# C2-OH conjugado a carbonila vinilogamente (nao casa fenol aromatico nem alcool alifatico).
ACID = Chem.MolFromSmarts("[OX2H1][C;!a]=[C;!a][C;!a]=O")

CONTROLE = ("CTRL_DQN", "Duroquinona (controle de redocking)", "CC1=C(C)C(=O)C(C)=C(C)C1=O", False)


def desprotonar(mol):
    hits = mol.GetSubstructMatches(ACID)
    if len(hits) != 1:
        raise ValueError(f"esperava 1 OH acido, achei {len(hits)}")
    rw = Chem.RWMol(mol)
    o = rw.GetAtomWithIdx(hits[0][0])
    o.SetFormalCharge(-1)
    Chem.SanitizeMol(rw)
    return rw.GetMol()


def para_pdbqt(mol, seed=42):
    mh = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mh, randomSeed=seed) != 0:
        raise RuntimeError("falha no embed 3D")
    if AllChem.MMFFOptimizeMolecule(mh, maxIters=2000) == -1:
        AllChem.UFFOptimizeMolecule(mh, maxIters=2000)
    setups = MoleculePreparation().prepare(mh)
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise RuntimeError(f"Meeko: {err}")
    ntors = pdbqt.count("BRANCH") // 2
    return pdbqt, ntors


def main():
    itens = []
    with open(SRC, encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            tem_anion = r["estado_ionizacao_pH7_4"].startswith("anion")
            itens.append((r["id"], r["nome"], r["smiles"], tem_anion))
    itens.append(CONTROLE)

    manifest = []
    for cid, nome, smi, tem_anion in itens:
        mol = Chem.MolFromSmiles(smi)
        estados = [("neutro", mol)]
        if tem_anion:
            estados.append(("anion", desprotonar(mol)))
        for estado, m in estados:
            lig_id = f"{cid}_{estado}"
            pdbqt, ntors = para_pdbqt(m)
            path = os.path.join(OUTD, f"{lig_id}.pdbqt")
            with open(path, "w") as fh:
                fh.write(pdbqt)
            manifest.append(dict(
                ligand_id=lig_id, id=cid, nome=nome, estado=estado,
                carga_formal=Chem.GetFormalCharge(m), torsoes=ntors,
                smiles_docado=Chem.MolToSmiles(m), pdbqt=os.path.relpath(path, HERE)))
            print(f"{lig_id:18} carga {Chem.GetFormalCharge(m):+d}  torsoes {ntors}  {nome}")

    with open(os.path.join(OUTD, "manifest.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(manifest[0].keys()))
        w.writeheader()
        w.writerows(manifest)
    print(f"\n{len(manifest)} ligantes (neutro + anion) -> {OUTD}")


if __name__ == "__main__":
    sys.exit(main())
