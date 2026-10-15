#!/usr/bin/env python3
"""Etapa 2: docking grosso da biblioteca no NQO1, com o MESMO protocolo do repositorio lapachol-nqo1.

Parametros copiados de 02_docking/vina_config.txt:
  caixa 22 x 22 x 22 A centrada em (19.754, -6.423, 5.313) [sitio da cadeia A, FAD retido],
  exhaustiveness 32, num_modes 9, energy_range 4, seed 42. Receptor rigido (receptor.pdbqt).

Uso:
  python 02_dock_library.py                       # tudo o que esta em ligands/manifest.csv
  python 02_dock_library.py --only L00_neutro,CTRL_DQN_neutro --cpu 2
  python 02_dock_library.py --vina /caminho/vina  # se o executavel nao estiver no PATH

Retoma de onde parou: pula ligantes cujo poses/<id>_out.pdbqt ja existe (use --force para refazer).
O escore do Vina ordena poses; NAO e energia livre de ligacao.
"""
import argparse
import csv
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))

BOX = dict(center_x=19.754, center_y=-6.423, center_z=5.313, size_x=22, size_y=22, size_z=22)
SEARCH = dict(exhaustiveness=32, num_modes=9, energy_range=4, seed=42)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vina", default="vina")
    ap.add_argument("--receptor", default=os.path.join(HERE, "receptor", "receptor.pdbqt"))
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--only", default="", help="ids separados por virgula (ex.: L00_neutro,L07_neutro)")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--seed", type=int, default=SEARCH["seed"],
                    help="semente do Vina; saidas com semente != 42 ganham sufixo _s<seed>")
    ap.add_argument("--budget", type=float, default=0,
                    help="segundos; nao inicia novo ligante depois disso (0 = sem limite)")
    a = ap.parse_args()
    t_start = time.time()
    suf = "" if a.seed == SEARCH["seed"] else f"_s{a.seed}"

    with open(os.path.join(HERE, "ligands", "manifest.csv"), encoding="utf-8") as fh:
        manifest = list(csv.DictReader(fh))
    if a.only:
        want = {x.strip() for x in a.only.split(",")}
        manifest = [m for m in manifest if m["ligand_id"] in want]
        faltam = want - {m["ligand_id"] for m in manifest}
        if faltam:
            sys.exit(f"ids nao encontrados no manifest: {sorted(faltam)}")

    os.makedirs(os.path.join(HERE, "poses"), exist_ok=True)
    for i, m in enumerate(manifest, 1):
        lid = m["ligand_id"]
        out = os.path.join(HERE, "poses", f"{lid}{suf}_out.pdbqt")
        log = os.path.join(HERE, "poses", f"{lid}{suf}.log")
        if os.path.exists(out) and not a.force:
            print(f"[{i}/{len(manifest)}] {lid}: ja existe, pulando")
            continue
        if a.budget and time.time() - t_start > a.budget:
            print(f"[{i}/{len(manifest)}] orcamento de tempo esgotado; rode de novo para retomar")
            break
        cmd = [a.vina, "--receptor", a.receptor,
               "--ligand", os.path.join(HERE, m["pdbqt"]), "--out", out, "--cpu", str(a.cpu)]
        for k, v in {**BOX, **{**SEARCH, "seed": a.seed}}.items():
            cmd += [f"--{k}", str(v)]
        t0 = time.time()
        with open(log, "w") as fh:
            rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode
        if rc != 0:
            if os.path.exists(out):
                os.remove(out)
            print(f"[{i}/{len(manifest)}] {lid}: FALHOU (ver {log})")
            continue
        print(f"[{i}/{len(manifest)}] {lid}: ok em {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
