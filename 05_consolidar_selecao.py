#!/usr/bin/env python3
"""Consolida as três etapas (propriedades, docking, descritores xTB) e gera a seleção de candidatos.

Saídas:
  funil_consolidado.csv     uma linha por espécie (composto x estado), com o papel na seleção
  selecao_candidatos.md     regra de seleção, tabela de contrastes e pendências antes de gastar GPU

A seleção é uma PROPOSTA por contraste de hipótese, não um ranking de "melhores": nenhuma das
três etapas separa os análogos o bastante para isso (ver números no .md).
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
rd = lambda n: pd.read_csv(HERE / n, encoding="utf-8-sig")
c = lambda x, nd=2: f"{x:.{nd}f}".replace(".", ",")  # vírgula decimal nos textos

props = rd("analogos_lapachol.csv")
dock = rd("resultados_docking.csv")
xtb = rd("descritores_xtb.csv")

df = xtb[["id", "estado", "ip", "ea", "mu", "eta", "omega"]].merge(
    props[["id", "nome", "grupo", "massa", "logP_Crippen", "logD7_4_aprox", "TPSA", "HBD",
           "pct_ionizado_pH7_4"]], on="id", how="left").merge(
    dock[["id", "estado", "escore_modo1", "core_rmsd_modo1", "d_centroide_modo1", "d_N5_modo1",
          "n_modos_tipo_ref", "n_modos"]], on=["id", "estado"], how="left")

# (id, estado, ordem, papel, pergunta que responde, o que contaria contra / ressalva)
SEL = [
    ("L00", "neutro", "0", "pré-requisito: lapachol neutro (conformação A) convergido",
     "Ponto de comparação de todos os outros. A rodada 1 não convergiu com 5 ns por janela; sem o lapachol "
     "neutro convergido, nenhum contraste abaixo é interpretável.",
     "Se o perfil continuar derivando com 15 a 20 ns por janela, o protocolo (campo de força, coordenada ou "
     "tempo) precisa de revisão antes de qualquer análogo."),
    ("L00", "anion", "1", "efeito da ionização (mesma molécula)",
     "O que muda no perfil quando o lapachol entra como ânion, a espécie dominante em água (pKa 5 a 6)?",
     "Se o perfil do ânion for igual ao do neutro dentro do ruído, o estado de protonação não é o fator "
     "limitante neste modelo."),
    ("L09", "neutro", "1", "par pareado sem acidez",
     "A troca de OH por OMe preserva forma, pose no sítio e eletrônica (docking a 0,07 Å da referência, EA "
     "quase igual): isola o papel do doador de H e da acidez no perfil de membrana.",
     "Se L09 ficar igual ao lapachol neutro, o doador de H não pesa; se ficar igual só ao ânion, a carga "
     "é o que importa."),
    ("L15", "anion", "2", "candidato de projeto (mais lipofílico)",
     "O Cl no anel aumenta logP e logD (a mudança no pKa não foi estimada): melhora a partição do ânion? "
     "A EA também sobe.",
     "Se o poço ficar mais fundo sem a barreira central cair, o ganho é de sequestro na membrana, "
     "não de permeação."),
    ("L01", "anion", "2", "controle negativo (polar)",
     "O lawsone é a molécula mais polar da lista: o perfil deve sair claramente pior que o do lapachol.",
     "Se não sair, o umbrella sampling não discrimina neste protocolo e os outros contrastes não valem."),
    ("L07", "neutro", "3", "controle positivo NQO1, isômero neutro",
     "Isômero do lapachol (mesma fórmula), neutro e sem torções livres: o perfil depende da forma ou só da "
     "carga? É também substrato conhecido da NQO1.",
     "Ressalva: no docking o modo 1 não cai na pose de referência (RMSD do esqueleto 2,7 Å, comparando "
     "topologias diferentes), embora empilhe no FAD."),
    ("L05", "anion", "opc", "opcional: extremo lipofílico",
     "A cadeia geranil (logP 4,96, perto do limite de Lipinski) mostra o regime de sequestro, em que poço "
     "fundo não significa permeabilidade.",
     "Custo alto de amostragem (6 torções livres) para uma pergunta que L15 já sonda."),
]
papel = {(i, e): f"{o} - {p}" for i, e, o, p, _, _ in SEL}
df["papel_na_selecao"] = [papel.get((r.id, r.estado), "") for r in df.itertuples()]
cols = ["id", "nome", "grupo", "estado", "papel_na_selecao", "massa", "logP_Crippen", "logD7_4_aprox",
        "TPSA", "HBD", "pct_ionizado_pH7_4", "escore_modo1", "core_rmsd_modo1", "d_centroide_modo1",
        "d_N5_modo1", "n_modos_tipo_ref", "n_modos", "ip", "ea", "mu", "eta", "omega"]
df = df.sort_values(["id", "estado"]).reset_index(drop=True)[cols]
df.to_csv(HERE / "funil_consolidado.csv", index=False, encoding="utf-8-sig")

# ------------------------------------------------------------------ números para o texto
neu = df[df.estado == "neutro"]
piv = df.pivot(index="id", columns="estado", values="ea").dropna()
queda = (piv["neutro"] - piv["anion"])          # queda de EA ao ionizar (positiva)
faixa_neu = neu.ea.max() - neu.ea.min()
g = lambda i, e, k: df[(df.id == i) & (df.estado == e)].iloc[0][k]
acidos = set(piv.index)

# ruído atual do protocolo, medido na rodada 1 (B) de lapachol-membrane (varredura de descarte e blocos)
b1 = rd("reference/round1_B_discard.csv")
bl = rd("reference/round1_B_blocks.csv")
b1["barreira"] = b1.G_center_kJ_mol - b1.min_kJ_mol
deriva_poco = b1.min_kJ_mol.max() - b1.min_kJ_mol.min()
deriva_centro = b1.G_center_kJ_mol.max() - b1.G_center_kJ_mol.min()
deriva_barreira = b1.barreira.max() - b1.barreira.min()
bl3 = bl[bl.scheme == "3x1500ps"].G_center_kJ_mol
dif_blocos = bl3.max() - bl3.min()
NS_DIA = 344.0                       # ns/dia na L4 com o lapachol (README de lapachol-membrane)
dias_15ns = 37 * 15 / NS_DIA         # 37 janelas x 15 ns
dias_5ns = 37 * 5 / NS_DIA


def dados(i, e):
    r = df[(df.id == i) & (df.estado == e)].iloc[0]
    if i in acidos and e == "anion":
        lip = f"logD {c(r.logD7_4_aprox)} (logP do neutro {c(r.logP_Crippen)})"
    else:
        lip = f"logP {c(r.logP_Crippen)}"
    return (f"{lip}, TPSA {c(r.TPSA, 1)}, HBD {int(r.HBD)}; docking {c(r.escore_modo1)} "
            f"(modo 1 a {c(r.core_rmsd_modo1)} Å da referência; {int(r.n_modos_tipo_ref)}/{int(r.n_modos)} "
            f"modos tipo-referência); EA {c(r.ea)} eV")


linhas = [f"| {o} | {i} {g(i, e, 'nome')}, {e.replace('anion', 'ânion')} | {perg} | {dados(i, e)} | {contra} |"
          for i, e, o, p, perg, contra in SEL]

md = f"""# Seleção de candidatos para MD e ΔG(z)

Gerado por `05_consolidar_selecao.py` a partir de três passos: propriedades (RDKit), docking (Vina, protocolo do
`lapachol-nqo1`) e descritores eletrônicos (xTB, protocolo do `quinone-ml`). A tabela completa está em
`funil_consolidado.csv`. Tudo aqui é hipótese a testar.

## Por que a seleção é por contraste e não por ranking

- **Docking:** todos os compostos empilham na isoaloxazina (distância entre centroides de
  {c(df.d_centroide_modo1.min())} a {c(df.d_centroide_modo1.max())} Å) e os escores ficam entre
  {c(df.escore_modo1.min())} e {c(df.escore_modo1.max())}. O filtro grosso não elimina ninguém.
- **Descritores (neutros):** EA de {c(neu.ea.min())} a {c(neu.ea.max())} eV e ω de {c(neu.omega.min(), 1)} a
  {c(neu.omega.max(), 1)} eV, uma faixa estreita. A validação do xTB contra o DFT existe só para o lapachol
  (EA 3,677 contra 3,770 eV), então diferenças dessa ordem entre análogos são hipótese.
- **Ionização:** passar de neutro a ânion reduz a EA em {c(queda.mean())} eV em média
  ({c(queda.min())} a {c(queda.max())} eV), cerca do dobro da faixa total entre os neutros
  ({c(faixa_neu)} eV). Sobre a EA, o estado de protonação pesa mais que qualquer substituinte da lista. Os
  números do ânion são exploratórios: a solvatação implícita de ânions é menos confiável e o dataset do
  `quinone-ml` só tem espécies neutras.

Por isso cada composto abaixo entra para responder uma pergunta que o lapachol sozinho não responde.

## Proposta

Ordem: 0 = pré-requisito, 1 = primeiro, 2 = segundo, 3 = terceiro, opc = só se sobrar tempo de GPU. Os
sistemas de ordem 1 em diante só devem rodar depois do pré-requisito (ver a sequência abaixo).

| Ordem | Sistema | Pergunta | Dados que sustentam | O que contaria contra / ressalva |
|---|---|---|---|---|
{chr(10).join(linhas)}

Desenho: lapachol (OH) e L09 (OMe) em neutro, mais o lapachol ânion, formam um 2 × 2 incompleto (o OMe não
ioniza) que separa o efeito do doador de H do efeito da carga. L09 neutro contra lapachol ânion é a comparação de
projeto: o que se ganha bloqueando a acidez.

## Sequência e custo

(a) **Quantificar A contra B** (só CPU): rodar `scripts/analysis_convergence.py` nas janelas de A com os mesmos
    descartes que B e comparar com `compare_rounds.py`. Hoje o "poço mais raso em A" ainda não tem número.
(b) **Convergir o lapachol neutro (conformação A)**: estágio 2, 15 a 20 ns por janela (proposta), descartando o
    início. Critério de parada proposto: o perfil mudar menos que o erro de bootstrap (2,5 kJ/mol) ao descartar
    mais 1 ns.
(c) **Lapachol ânion.** Sistema novo: parâmetros CGenFF do ânion, equilibração e puxada. O ânion não tem o OH, então
    some a armadilha da torção que atrasou a rodada 1. Em compensação, espécies carregadas costumam pedir mais
    amostragem (água, contra-íon, deformação da membrana); isso é esperado, não medido aqui.
(d) **L09 neutro.** Também sistema novo; sem OH, mas a metoxila tem a própria torção, que merece o mesmo
    escaneamento xTB contra CGenFF que o OH recebeu.
(e) Só depois, conforme os resultados: L15 ânion, L01 ânion e L07.

Custo, estimado com os seus 344 ns/dia na L4 e as 37 janelas: cerca de {c(dias_15ns, 1)} dia de GPU por sistema a
15 ns por janela (e {c(dias_5ns, 1)} dia a 5 ns). Cinco sistemas a 15 ns são cerca de {c(5 * dias_15ns, 0)} dias de L4,
fora equilibração e puxada. É esse orçamento que justifica fazer (a) e (b) antes de qualquer análogo.

## Antes de gastar GPU

1. **Estado de protonação (confirmado).** As rodadas A e B são do lapachol neutro. O ânion, a espécie dominante
   em água, ainda não foi simulado.
2. **Ruído atual do protocolo** (rodada 1, 5 ns por janela, lida em `lapachol-membrane`): ao descartar mais do
   início, o poço vai de {c(b1.min_kJ_mol.iloc[0], 1)} a {c(b1.min_kJ_mol.min(), 1)} kJ/mol (deriva de
   {c(deriva_poco, 1)}), o centro varia {c(deriva_centro, 1)} e a barreira {c(deriva_barreira, 1)} kJ/mol; os blocos
   de tempo diferem até {c(dif_blocos, 1)} kJ/mol no centro. O bootstrap (até 2,5 kJ/mol) não enxerga isso.
   Até convergir, contrastes menores que essas derivas não são resolvíveis, e isso vale também para A contra B.
3. **pKa.** O intervalo de 5 a 6 é heurístico por classe. Confirmar com o DFT do ânion que está pendente e, se
   possível, com uma ferramenta de pKa.
4. **Parâmetros do ânion.** O CGenFF do ânion precisa de checagem de penalidade e de carga. No `lapachol-nqo1`
   o lapachol neutro teve penalidade 28,5, ainda sem a validação recomendada.
5. **Carga fixa.** Com protonação fixa, o ânion não pode se neutralizar ao entrar na região apolar da
   bicamada, onde o pKa tende a subir. Neutro e ânion são dois limites do problema, não a resposta.
6. **Permeabilidade.** ΔG(z) sozinho não é permeabilidade sem o perfil de difusão. Poço fundo pode ser
   sequestro, não vantagem.

## Docking e MD do complexo

O docking não discrimina, então não serve de critério para escolher quem vai para MD do complexo. Para fechar a
pendência do estado de protonação na NQO1, rodar o complexo com o lapachol ânion (a pose de docking dele também
é tipo-referência, a {c(g("L00", "anion", "core_rmsd_modo1"))} Å) e comparar com os 10,7 ns do neutro. L09 e L15
só entram no MD do complexo se passarem pela etapa de membrana.
"""
(HERE / "selecao_candidatos.md").write_text(md, encoding="utf-8")
print(f"funil_consolidado.csv: {len(df)} linhas | selecao_candidatos.md: {len(md.splitlines())} linhas")
print(f"EA neutros {neu.ea.min():.3f}-{neu.ea.max():.3f} (faixa {faixa_neu:.3f}); queda de EA ao ionizar "
      f"{queda.mean():.3f} ({queda.min():.3f} a {queda.max():.3f}); razão {queda.mean() / faixa_neu:.2f}")
