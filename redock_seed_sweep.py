#!/usr/bin/env python3
"""Controle de redocking da duroquinona em várias sementes do Vina, com dois preparos do ligante.

Mostra o quanto o modo 1 do redocking depende da semente e do preparo, e se a pose tipo-cristal
(RMSD < 2,0 Å contra 1DXO) aparece entre os 9 modos. Mesmo protocolo de 02_dock_library.py.

Preparos:
  meeko  : ligands/CTRL_DQN_neutro.pdbqt (gerado por 01_prepare_ligands.py)
  obabel : `obabel reference/dqn_clean.pdb -p 7.4 --partialcharge gasteiger` (como em lapachol-nqo1);
           só roda se o executável `obabel` existir.

Uso: python redock_seed_sweep.py --vina /caminho/vina [--seeds 42 1 2 3 4 5 6 7 8 9] [--cpu 2]
Saídas temporárias em scratch/ (ignorado pelo git).
"""
import argparse
import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("an", HERE / "03_analyze_poses.py")
an = importlib.util.module_from_spec(spec)
spec.loader.exec_module(an)

BOX = ["--center_x", "19.754", "--center_y", "-6.423", "--center_z", "5.313",
       "--size_x", "22", "--size_y", "22", "--size_z", "22",
       "--exhaustiveness", "32", "--num_modes", "9", "--energy_range", "4"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vina", default="vina")
    ap.add_argument("--seeds", type=int, nargs="+", default=[42, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    ap.add_argument("--cpu", type=int, default=2)
    a = ap.parse_args()

    scratch = HERE / "scratch"
    scratch.mkdir(exist_ok=True)
    ref = an.read_pdb_heavy(HERE / "reference" / "dqn_clean.pdb")
    preps = {"meeko": HERE / "ligands" / "CTRL_DQN_neutro.pdbqt"}
    if shutil.which("obabel"):
        ob = scratch / "dqn_obabel.pdbqt"
        subprocess.run(["obabel", str(HERE / "reference" / "dqn_clean.pdb"), "-O", str(ob),
                        "-p", "7.4", "--partialcharge", "gasteiger"],
                       check=True, capture_output=True)
        preps["obabel"] = ob
    else:
        print("obabel não encontrado: só o preparo meeko será testado\n")

    for nome, lig in preps.items():
        print(f"== preparo {nome} ==")
        print("seed | escore modo 1 | RMSD modo 1 | menor RMSD (modo) | escore desse modo")
        ok = tipo_cristal = 0
        for sd in a.seeds:
            out = scratch / f"dq_{nome}_{sd}.pdbqt"
            subprocess.run([a.vina, "--receptor", str(HERE / "receptor" / "receptor.pdbqt"),
                            "--ligand", str(lig), "--out", str(out), "--cpu", str(a.cpu),
                            "--seed", str(sd)] + BOX,
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            mods = an.read_pdbqt_models(out)
            rm = [an.sym_rmsd(m["atoms"], ref) for m in mods]
            sc = [m["score"] for m in mods]
            ib = int(np.nanargmin(rm))
            ok += rm[0] < 2.0
            tipo_cristal += rm[ib] < 2.0
            print(f"{sd:4d} | {sc[0]:13.2f} | {rm[0]:11.2f} | {rm[ib]:7.2f} ({ib + 1})       | {sc[ib]:6.2f}")
        n = len(a.seeds)
        print(f"modo 1 com RMSD < 2,0 Å: {ok}/{n} sementes; pose tipo-cristal entre os modos: {tipo_cristal}/{n}\n")


if __name__ == "__main__":
    sys.exit(main())
