#!/usr/bin/env python3
"""Compara duas rodadas (ex.: A e B) de umbrella sampling a partir da varredura de descarte.

Entrada: dois CSV no formato que scripts/analysis_convergence.py (lapachol-membrane) escreve em
results/convergence_discard.csv, com as colunas
  data_from_ps, data_to_ps, min_kJ_mol, z_min_nm, G_center_kJ_mol [, barrier_kJ_mol]
Para a rodada B existe results/convergence_stage1_discard.csv no repositorio. Para A, rode
analysis_convergence.py nas janelas de A com os MESMOS descartes e o mesmo tmax.

Metricas (mesma definicao do README do lapachol-membrane, zero do perfil na agua):
  poco    = minimo do perfil (kJ/mol; mais negativo = mais fundo)
  centro  = G(z = 0)
  barreira = centro - poco
Delta = A - B; poco com delta > 0 significa "A mais raso que B".

Leitura (heuristica de leitura, nao teste estatistico):
  deriva = (max - min) de cada metrica ao longo da varredura de descarte, dentro de cada rodada.
  O erro de bootstrap do WHAM nao enxerga essa deriva (README: bootstrap <= 2.5 kJ/mol).
  - diferenca resolvida : delta com o mesmo sinal em todos os descartes e |delta| menor > deriva maior
  - sugestiva           : mesmo sinal em todos os descartes, mas dentro da deriva
  - nao resolvida       : sinal muda entre descartes
  Uma rodada so e dada como estavel se a metrica varia <= --tol nos ultimos 3 descartes.

Uso: python compare_rounds.py --A convergence_discard_A.csv --B convergence_stage1_discard.csv
"""
import argparse
import sys

import pandas as pd

NOMES = {"min_kJ_mol": "poco", "G_center_kJ_mol": "centro", "barrier_kJ_mol": "barreira"}


def carregar(path):
    d = pd.read_csv(path)
    d["barrier_kJ_mol"] = d["G_center_kJ_mol"] - d["min_kJ_mol"]  # sempre recalculada, mesma definicao nas duas rodadas
    return d.set_index("data_from_ps")[list(NOMES)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--A", required=True)
    ap.add_argument("--B", required=True)
    ap.add_argument("--tol", type=float, default=2.5, help="kJ/mol; criterio de estabilidade (default: bootstrap max do README)")
    a = ap.parse_args()

    A, B = carregar(a.A), carregar(a.B)
    comum = sorted(set(A.index) & set(B.index))
    if len(comum) < 3:
        sys.exit(f"descartes em comum insuficientes ({comum}); use os mesmos --discards nas duas rodadas")
    A, B = A.loc[comum], B.loc[comum]
    delta = A - B

    print(f"descartes em comum (ps): {comum}\n")
    print("Delta = A - B (kJ/mol)")
    print(delta.rename(columns=NOMES).round(1).to_string())

    print("\nDeriva dentro de cada rodada (max - min ao longo dos descartes) e estabilidade nos 3 ultimos:")
    print(f"{'metrica':9} {'deriva A':>9} {'deriva B':>9} {'ult3 A':>7} {'ult3 B':>7}")
    veredito = {}
    for col, nome in NOMES.items():
        dA, dB = A[col].max() - A[col].min(), B[col].max() - B[col].min()
        uA, uB = A[col].iloc[-3:].max() - A[col].iloc[-3:].min(), B[col].iloc[-3:].max() - B[col].iloc[-3:].min()
        print(f"{nome:9} {dA:9.1f} {dB:9.1f} {uA:7.1f} {uB:7.1f}")
        s = delta[col]
        mesmo_sinal = bool((s > 0).all() or (s < 0).all())
        deriva = max(dA, dB)
        if mesmo_sinal and s.abs().min() > deriva:
            v = "diferenca resolvida"
        elif mesmo_sinal:
            v = "sugestiva (mesmo sinal, mas dentro da deriva)"
        else:
            v = "nao resolvida (o sinal muda entre descartes)"
        est = []
        if uA > a.tol:
            est.append("A")
        if uB > a.tol:
            est.append("B")
        veredito[nome] = (v, est, s.mean())

    print()
    for nome, (v, est, m) in veredito.items():
        extra = (f"; rodada(s) ainda instavel(is) nos ultimos 3 descartes (> {a.tol} kJ/mol): {', '.join(est)}"
                 if est else "")
        print(f"{nome:9}: {v}; delta medio {m:+.1f} kJ/mol{extra}")
    if any(e for _, e, _ in veredito.values()):
        print("\nAviso: com rodada(s) instavel(is), mesmo uma diferenca 'resolvida' descreve o dado atual, "
              "nao o perfil convergido.")


if __name__ == "__main__":
    main()
