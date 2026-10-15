# lapachol-analogs

**Funnel for the rational design of lapachol analogs: calculated properties → docking → electronic descriptors → MD and ΔG(z) for a few candidates.**
Part of my computational chemistry portfolio ([msouza-compchem](https://github.com/msouza-compchem)). Work in progress: every step is documented with its parameters, scripts, results and limitations.

*Resumo em português: funil em etapas para desenhar análogos do lapachol com boa permeabilidade em membrana e bom encaixe na NQO1 humana. Este repositório reúne a biblioteca de 17 compostos, o docking grosso (protocolo do `lapachol-nqo1`), os descritores eletrônicos (protocolo do `quinone-ml`) e a seleção de poucos candidatos para dinâmica molecular e ΔG(z) (`lapachol-membrane`). Os resultados são hipóteses a testar, não conclusões. Os ensaios experimentais são previstos em colaboração com um grupo de bancada (ver `PREDICTIONS.md`).*

## Context

Three earlier studies feed this one:

| Repository | What it provides here |
|---|---|
| [`lapachol-nqo1`](https://github.com/msouza-compchem/lapachol-nqo1) | DFT, docking validated by redocking, 10.7 ns of MD of the complex (protocol reused in stage 2) |
| [`quinone-ml`](https://github.com/msouza-compchem/quinone-ml) | GFN2-xTB vertical IP/EA protocol for quinones (reused in stage 3) |
| [`lapachol-membrane`](https://github.com/msouza-compchem/lapachol-membrane) | Umbrella sampling of lapachol in POPC (stage 4; round 1 not converged) |

## Question

> Which lapachol analogs could both cross a membrane and fit the catalytic site of human NQO1, and how do we choose the few worth the cost of MD and ΔG(z)?

**What this project does not answer:** docking scores are not binding free energies; ΔG(z) alone is not permeability (that also needs the diffusion profile D(z)); a classical force field with one lipid composition is not a cell membrane; nothing here says a compound is an NQO1 substrate (that requires measurement).

## Funnel and status

| Stage | Method | Script | Status |
|---|---|---|---|
| 1. Properties | RDKit: mass, Crippen logP, TPSA, HBD/HBA, heuristic pKa per class, approximate logD at pH 7.4 | `00_propriedades.py` | Done, 17 compounds |
| 2. Docking | AutoDock Vina 1.2.7, protocol of `lapachol-nqo1`; neutral and anion for the 10 acidic compounds | `01_prepare_ligands.py`, `02_dock_library.py`, `03_analyze_poses.py` | Done; coarse filter only |
| 3. Electronic descriptors | GFN2-xTB (ALPB water, `--vipea`), protocol of `quinone-ml` | `04_descritores_xtb.py` | Done; anion values exploratory |
| 4. Membrane MD, ΔG(z) | Umbrella sampling in POPC (GROMACS) | in `lapachol-membrane` | Not started for analogs; the lapachol reference is not converged |
| 5. Experiments | pKa, logD, POPC liposome partition, NQO1 kinetics | `PREDICTIONS.md` | Proposed, to be done with a collaborating group |

`05_consolidar_selecao.py` merges stages 1–3 into `funil_consolidado.csv` and writes `selecao_candidatos.md`. `compare_rounds.py` compares two umbrella-sampling rounds (for example conformations A and B of the lapachol OH) while taking the drift of each round into account.

## Library

17 compounds in six groups: the reference (lapachol), side-chain variants (lawsone, norlapachol, dihydrolapachol, butyl, geranyl, hydroxylated chain), cyclised isomers (α- and β-lapachone), acidity-blocked variants (2-methoxy, acetate), OH bioisosteres (NH₂, NHMe, Me) and ring-substituted variants (7-F, 6-Cl, 5-OH). Full table with SMILES and the hypothesis behind each one: `analogos_lapachol.csv`.

## Results so far

- **Protonation state matters more than any substituent.** Ionising a 2-hydroxy-3-alkyl compound lowers the calculated vertical electron affinity by 0.97 eV on average (0.88 to 1.02 eV), about twice the whole range among the neutral compounds (3.42 to 3.92 eV, a 0.50 eV span). The earlier lapachol docking and MD used the neutral form; the anion is expected to dominate at pH 7.4 if the pKa is 5 to 6, which is a class-level estimate and has not been measured.
- **Docking does not discriminate.** All compounds stack on the FAD isoalloxazine ring (ring-centroid distance 3.64 to 4.28 Å) and scores span −6.89 to −5.30, with seed-to-seed variation of about 0.03. Pose and distance to FAD are reported together with the score.
- **xTB reproduces the DFT reference for lapachol:** vertical EA 3.677 eV against 3.770 eV (adiabatic DFT). This is the only DFT point of comparison.
- **Selection by contrast, not by ranking.** None of the three stages separates the analogs well enough to rank them, so `selecao_candidatos.md` proposes compounds that each answer a different question (ionisation, hydrogen-bond donor, lipophilicity, negative and positive controls).

## Limitations

- **Redocking is fragile.** The duroquinone control is reproduced at 1.11 Å in mode 1 with the original preparation (Open Babel, seed 42), but mode 1 is below 2.0 Å in only 3 of 10 seeds with that preparation and in 0 of 10 with Meeko. A crystal-like pose (0.87 to 1.18 Å) is among the 9 modes in all 20 runs, with top scores within about 0.1 of each other. The analysis therefore looks at all modes, not only the first.
- **Neutral and anion differ only in atom typing in Vina**, which has no electrostatic term. Do not read the small score differences as a preference.
- **pKa and logD are heuristics** (Crippen logP plus a pKa interval per compound class). They give an order of magnitude only.
- **Single DFT reference** for the xTB descriptors; differences of a few tenths of an eV between analogs are within the method error.
- **Rigid receptor**, one site (chain A), one conformer of 1D4A.
- **Anion descriptors are exploratory:** implicit solvation of anions is less reliable and the `quinone-ml` dataset contains only neutral species.
- **Membrane stage:** the lapachol reference run has not converged (round 1, 5 ns per window; drift documented in `lapachol-membrane`). No contrast between analogs is interpretable before it does. Fixed protonation also prevents the anion from neutralising inside the bilayer.

## Repository layout

```
00_propriedades.py ... 05_consolidar_selecao.py   pipeline, one script per stage
compare_rounds.py        comparison of two umbrella-sampling rounds
redock_seed_sweep.py     duroquinone redocking over seeds and ligand preparations
                         (output in reference/redock_seed_sweep.txt)
analogos_lapachol.csv    library and calculated properties (stage 1)
resultados_docking.csv   docking metrics per ligand (stage 2)
descritores_xtb.csv      xTB descriptors per species (stage 3)
funil_consolidado.csv    the three stages merged, with the role in the selection
selecao_candidatos.md    proposed selection, sequence and cost
PREDICTIONS.md           predictions to be frozen before any measurement
ligands/                 PDBQT of every ligand and manifest.csv
poses/                   Vina output and log for every ligand
receptor/                receptor.pdbqt (from lapachol-nqo1)
reference/               reference poses and round 1 convergence tables
```

## Reproducing

```bash
python 00_propriedades.py analogos_lapachol.csv      # stage 1
python 01_prepare_ligands.py                         # SMILES -> PDBQT (neutral and anion)
python 02_dock_library.py --cpu 4 --vina /path/to/vina   # ~10 min on 2 cores
python 03_analyze_poses.py                           # -> resultados_docking.csv
python 04_descritores_xtb.py --xtb /path/to/xtb --quinone-ml /path/to/quinone-ml
python 05_consolidar_selecao.py                      # -> funil_consolidado.csv, selecao_candidatos.md
```

`02_dock_library.py` resumes where it stopped and accepts `--seed` for replicas. Versions used: RDKit 2026.9.1, Meeko 0.8.0, gemmi 0.7.5, AutoDock Vina 1.2.7, xTB 6.7.1, Python 3.13. `environment.yml` lists the packages but has not been tested with conda; the runs above used pip and the standalone Vina and xTB binaries.

## Provenance of copied files

- `receptor/receptor.pdbqt`, `reference/best_pose.pdb`, `reference/dqn_clean.pdb`: from `lapachol-nqo1` (commit `a52d367`).
- `reference/round1_B_discard.csv`, `reference/round1_B_blocks.csv`: from `lapachol-membrane` (commit `7f22839`, files `results/convergence_stage1_discard.csv` and `convergence_stage1_blocks.csv`).

## Collaboration

The experimental arm (synthesis of a few compounds, pKa, logD, liposome partition, NQO1 kinetics) is meant to be done with a bench group. The plan is to freeze the predictions in `PREDICTIONS.md` with a dated tag before any measurement, so that results in either direction are informative. Contact through GitHub issues.

## Licence

MIT, see [`LICENSE`](LICENSE).
