#!/usr/bin/env python3
"""Biblioteca de analogos do lapachol: propriedades calculadas (RDKit).

Etapa 1 do funil: propriedades calculadas -> docking -> descritores eletronicos -> MD / DG(z).
Uso: python gerar_analogos.py saida.csv

Limites conhecidos (ver colunas e notas):
- logP e Crippen (aditivo por atomos), forma NEUTRA. Nao e medido.
- pKa e heuristico por classe (intervalo), NAO calculado. Confirmar com ferramenta de pKa.
- logD7.4 = logP + log10(fracao neutra), aproximacao grosseira (ignora particao de par ionico).
"""
import csv
import sys

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

PH = 7.4
OUT = sys.argv[1]

# Intervalos de pKa por classe (heuristicos).
ACIDO = (5.0, 6.0)   # 2-hidroxi-3-alquil-1,4-naftoquinona (faixa informada para o lapachol)
LAWSONE = (3.5, 4.5)  # sem substituinte alquila no C3 (mais acido)
NEUTRO = None

# (id, nome, grupo, SMILES, pKa, racional/hipotese)
LIB = [
    ("L00", "Lapachol (referencia)", "controle",
     "CC(C)=CCC1=C(O)C(=O)c2ccccc2C1=O", ACIDO,
     "Referencia do funil."),
    ("L01", "Lawsone (sem prenil)", "cadeia lateral",
     "OC1=CC(=O)c2ccccc2C1=O", LAWSONE,
     "Remove a cadeia lipofilica; mostra quanto do perfil de membrana vem do prenil."),
    ("L02", "Norlapachol", "cadeia lateral",
     "CC(C)=CC1=C(O)C(=O)c2ccccc2C1=O", ACIDO,
     "Cadeia um carbono mais curta, alceno conjugado ao anel."),
    ("L03", "Dihidrolapachol", "cadeia lateral",
     "CC(C)CCC1=C(O)C(=O)c2ccccc2C1=O", ACIDO,
     "Prenil saturado; testa o papel da insaturacao (sitio de oxidacao)."),
    ("L04", "3-butil-2-hidroxi-NQ", "cadeia lateral",
     "CCCCC1=C(O)C(=O)c2ccccc2C1=O", ACIDO,
     "Serie de alquilas lineares; separa efeito de comprimento do efeito de ramificacao."),
    ("L05", "Geranil-analogo (C10)", "cadeia lateral",
     "CC(C)=CCCC(C)=CCC1=C(O)C(=O)c2ccccc2C1=O", ACIDO,
     "Cadeia longa; testa o limite superior de lipofilia e a ancoragem na bicamada."),
    ("L06", "3-(3-hidroxi-3-metilbutil)", "cadeia lateral",
     "CC(C)(O)CCC1=C(O)C(=O)c2ccccc2C1=O", ACIDO,
     "Hidroxilacao da cadeia; reduz logP e adiciona um doador de H."),
    ("L07", "beta-Lapachona", "ciclizado",
     "CC1(C)CCC2=C(O1)c1ccccc1C(=O)C2=O", NEUTRO,
     "o-Quinona ciclizada, neutra; substrato conhecido da NQO1 (controle positivo plausivel)."),
    ("L08", "alfa-Lapachona", "ciclizado",
     "CC1(C)CCC2=C(O1)C(=O)c1ccccc1C2=O", NEUTRO,
     "p-Quinona ciclizada, neutra; isomero da beta."),
    ("L09", "2-Metoxi-lapachol", "bloqueio de acidez",
     "CC(C)=CCC1=C(OC)C(=O)c2ccccc2C1=O", NEUTRO,
     "Elimina a ionizacao e o doador de H; testa permeabilidade da forma neutra vs encaixe no sitio."),
    ("L10", "Acetato de lapachol", "bloqueio de acidez",
     "CC(C)=CCC1=C(OC(C)=O)C(=O)c2ccccc2C1=O", NEUTRO,
     "Ester neutro (pro-farmaco em hipotese); hidrolise esperada em meio biologico."),
    ("L11", "2-Amino-3-prenil-NQ", "bioisostero de OH",
     "CC(C)=CCC1=C(N)C(=O)c2ccccc2C1=O", NEUTRO,
     "OH -> NH2: mantem doador de H, perde a acidez."),
    ("L12", "2-(Metilamino)-3-prenil-NQ", "bioisostero de OH",
     "CC(C)=CCC1=C(NC)C(=O)c2ccccc2C1=O", NEUTRO,
     "Variante N-metilada; mais lipofilica."),
    ("L13", "2-Metil-3-prenil-NQ", "bioisostero de OH",
     "CC(C)=CCC1=C(C)C(=O)c2ccccc2C1=O", NEUTRO,
     "Sem grupo acido nem doador de H; p-quinona do tipo vitamina K."),
    ("L14", "7-Fluoro-lapachol", "anel aromatico",
     "CC(C)=CCC1=C(O)C(=O)c2cc(F)ccc2C1=O", ACIDO,
     "F no anel benzenico; tende a reduzir o pKa e alterar o potencial redox (magnitude nao estimada)."),
    ("L15", "6-Cloro-lapachol", "anel aromatico",
     "CC(C)=CCC1=C(O)C(=O)c2ccc(Cl)cc2C1=O", ACIDO,
     "Cl no anel benzenico; aumenta logP; mesma tendencia de pKa que o F."),
    ("L16", "5-Hidroxi-lapachol", "anel aromatico",
     "CC(C)=CCC1=C(O)C(=O)c2cccc(O)c2C1=O", ACIDO,
     "Tipo juglona; OH peri ao carbonila pode formar ligacao de H intramolecular."),
]


def frac_ion(pka, ph=PH):
    return 1.0 / (1.0 + 10 ** (pka - ph))


params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
pains = FilterCatalog(params)

rows = []
seen = {}
ref_logp = None

for cid, nome, grupo, smi, pka, why in LIB:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise SystemExit(f"SMILES invalido: {cid} {smi}")
    can = Chem.MolToSmiles(mol)
    if can in seen:
        raise SystemExit(f"Duplicata: {cid} == {seen[can]}")
    seen[can] = cid

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    heavy = mol.GetNumHeavyAtoms()

    if pka is None:
        pka_txt, ion_txt, logd, estado = "n/a", "~0", logp, "neutro"
    else:
        lo, hi = pka
        mid = (lo + hi) / 2
        f_lo, f_hi = frac_ion(hi), frac_ion(lo)  # pKa maior -> menos ionizado
        pka_txt = f"{lo:.1f}-{hi:.1f}"
        ion_txt = f"{100*f_lo:.1f}-{100*f_hi:.1f}"
        import math
        logd = logp + math.log10(1.0 - frac_ion(mid))
        estado = "anion (>95%)" if f_lo > 0.95 else "misto"

    lip_viol = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    veber_ok = tpsa <= 140 and rotb <= 10
    pains_hits = ";".join(sorted({e.GetDescription() for e in pains.GetMatches(mol)})) or "-"

    if cid == "L00":
        ref_logp = logp
        assert formula == "C15H14O3", formula
        assert abs(mw - 242.27) < 0.05, mw

    rows.append(dict(
        id=cid, nome=nome, grupo=grupo, smiles=smi, formula=formula,
        massa=round(mw, 2), heavy_atoms=heavy,
        logP_Crippen=round(logp, 2), TPSA=round(tpsa, 1),
        HBD=hbd, HBA=hba, rot_bonds=rotb,
        pKa_estimado_heuristico=pka_txt, pct_ionizado_pH7_4=ion_txt,
        estado_ionizacao_pH7_4=estado, logD7_4_aprox=round(logd, 2),
        Lipinski_violacoes=lip_viol, Veber_ok="sim" if veber_ok else "nao",
        PAINS=pains_hits, racional_hipotese=why,
    ))

for r in rows:
    r["delta_logP_vs_lapachol"] = round(r["logP_Crippen"] - ref_logp, 2)

fields = ["id", "nome", "grupo", "smiles", "formula", "massa", "heavy_atoms",
          "logP_Crippen", "delta_logP_vs_lapachol", "TPSA", "HBD", "HBA", "rot_bonds",
          "pKa_estimado_heuristico", "pct_ionizado_pH7_4", "estado_ionizacao_pH7_4",
          "logD7_4_aprox", "Lipinski_violacoes", "Veber_ok", "PAINS", "racional_hipotese"]

with open(OUT, "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)

hdr = f"{'id':4} {'nome':30} {'MW':>7} {'logP':>5} {'TPSA':>5} {'HBD':>3} {'HBA':>3} {'RB':>2} {'%ion':>9} {'logD':>5}  PAINS"
print(hdr)
for r in rows:
    print(f"{r['id']:4} {r['nome'][:30]:30} {r['massa']:7.2f} {r['logP_Crippen']:5.2f} "
          f"{r['TPSA']:5.1f} {r['HBD']:3d} {r['HBA']:3d} {r['rot_bonds']:2d} "
          f"{r['pct_ionizado_pH7_4']:>9} {r['logD7_4_aprox']:5.2f}  {r['PAINS']}")
print(f"\n{len(rows)} compostos; CSV: {OUT}")
