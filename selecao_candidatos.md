# Seleção de candidatos para MD e ΔG(z)

Gerado por `05_consolidar_selecao.py` a partir de três passos: propriedades (RDKit), docking (Vina, protocolo do
`lapachol-nqo1`) e descritores eletrônicos (xTB, protocolo do `quinone-ml`). A tabela completa está em
`funil_consolidado.csv`. Tudo aqui é hipótese a testar.

## Por que a seleção é por contraste e não por ranking

- **Docking:** todos os compostos empilham na isoaloxazina (distância entre centroides de
  3,64 a 4,28 Å) e os escores ficam entre
  -6,89 e -5,30. O filtro grosso não elimina ninguém.
- **Descritores (neutros):** EA de 3,42 a 3,92 eV e ω de 8,6 a
  11,2 eV, uma faixa estreita. A validação do xTB contra o DFT existe só para o lapachol
  (EA 3,677 contra 3,770 eV), então diferenças dessa ordem entre análogos são hipótese.
- **Ionização:** passar de neutro a ânion reduz a EA em 0,97 eV em média
  (0,88 a 1,02 eV), cerca do dobro da faixa total entre os neutros
  (0,50 eV). Sobre a EA, o estado de protonação pesa mais que qualquer substituinte da lista. Os
  números do ânion são exploratórios: a solvatação implícita de ânions é menos confiável e o dataset do
  `quinone-ml` só tem espécies neutras.

Por isso cada composto abaixo entra para responder uma pergunta que o lapachol sozinho não responde.

## Proposta

Ordem: 0 = pré-requisito, 1 = primeiro, 2 = segundo, 3 = terceiro, opc = só se sobrar tempo de GPU. Os
sistemas de ordem 1 em diante só devem rodar depois do pré-requisito (ver a sequência abaixo).

| Ordem | Sistema | Pergunta | Dados que sustentam | O que contaria contra / ressalva |
|---|---|---|---|---|
| 0 | L00 Lapachol (referencia), neutro | Ponto de comparação de todos os outros. A rodada 1 não convergiu com 5 ns por janela; sem o lapachol neutro convergido, nenhum contraste abaixo é interpretável. | logP 3,23, TPSA 54,4, HBD 1; docking -6,02 (modo 1 a 0,75 Å da referência; 6/8 modos tipo-referência); EA 3,68 eV | Se o perfil continuar derivando com 15 a 20 ns por janela, o protocolo (campo de força, coordenada ou tempo) precisa de revisão antes de qualquer análogo. |
| 1 | L00 Lapachol (referencia), ânion | O que muda no perfil quando o lapachol entra como ânion, a espécie dominante em água (pKa 5 a 6)? | logD 1,33 (logP do neutro 3,23), TPSA 54,4, HBD 1; docking -6,07 (modo 1 a 0,38 Å da referência; 4/9 modos tipo-referência); EA 2,67 eV | Se o perfil do ânion for igual ao do neutro dentro do ruído, o estado de protonação não é o fator limitante neste modelo. |
| 1 | L09 2-Metoxi-lapachol, neutro | A troca de OH por OMe preserva forma, pose no sítio e eletrônica (docking a 0,07 Å da referência, EA quase igual): isola o papel do doador de H e da acidez no perfil de membrana. | logP 3,32, TPSA 43,4, HBD 0; docking -5,97 (modo 1 a 0,07 Å da referência; 5/9 modos tipo-referência); EA 3,70 eV | Se L09 ficar igual ao lapachol neutro, o doador de H não pesa; se ficar igual só ao ânion, a carga é o que importa. |
| 2 | L15 6-Cloro-lapachol, ânion | O Cl no anel aumenta logP e logD (a mudança no pKa não foi estimada): melhora a partição do ânion? A EA também sobe. | logD 1,98 (logP do neutro 3,89), TPSA 54,4, HBD 1; docking -6,22 (modo 1 a 0,62 Å da referência; 6/9 modos tipo-referência); EA 2,91 eV | Se o poço ficar mais fundo sem a barreira central cair, o ganho é de sequestro na membrana, não de permeação. |
| 2 | L01 Lawsone (sem prenil), ânion | O lawsone é a molécula mais polar da lista: o perfil deve sair claramente pior que o do lapachol. | logD -1,89 (logP do neutro 1,51), TPSA 54,4, HBD 1; docking -5,45 (modo 1 a 1,57 Å da referência; 4/9 modos tipo-referência); EA 2,77 eV | Se não sair, o umbrella sampling não discrimina neste protocolo e os outros contrastes não valem. |
| 3 | L07 beta-Lapachona, neutro | Isômero do lapachol (mesma fórmula), neutro e sem torções livres: o perfil depende da forma ou só da carga? É também substrato conhecido da NQO1. | logP 2,75, TPSA 43,4, HBD 0; docking -6,89 (modo 1 a 2,66 Å da referência; 3/8 modos tipo-referência); EA 3,64 eV | Ressalva: no docking o modo 1 não cai na pose de referência (RMSD do esqueleto 2,7 Å, comparando topologias diferentes), embora empilhe no FAD. |
| opc | L05 Geranil-analogo (C10), ânion | A cadeia geranil (logP 4,96, perto do limite de Lipinski) mostra o regime de sequestro, em que poço fundo não significa permeabilidade. | logD 3,06 (logP do neutro 4,96), TPSA 54,4, HBD 1; docking -5,64 (modo 1 a 0,44 Å da referência; 5/9 modos tipo-referência); EA 2,68 eV | Custo alto de amostragem (6 torções livres) para uma pergunta que L15 já sonda. |

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

Custo, estimado com os seus 344 ns/dia na L4 e as 37 janelas: cerca de 1,6 dia de GPU por sistema a
15 ns por janela (e 0,5 dia a 5 ns). Cinco sistemas a 15 ns são cerca de 8 dias de L4,
fora equilibração e puxada. É esse orçamento que justifica fazer (a) e (b) antes de qualquer análogo.

## Antes de gastar GPU

1. **Estado de protonação (confirmado).** As rodadas A e B são do lapachol neutro. O ânion, a espécie dominante
   em água, ainda não foi simulado.
2. **Ruído atual do protocolo** (rodada 1, 5 ns por janela, lida em `lapachol-membrane`): ao descartar mais do
   início, o poço vai de -7,9 a -13,6 kJ/mol (deriva de
   5,7), o centro varia 12,6 e a barreira 7,0 kJ/mol; os blocos
   de tempo diferem até 25,5 kJ/mol no centro. O bootstrap (até 2,5 kJ/mol) não enxerga isso.
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
é tipo-referência, a 0,38 Å) e comparar com os 10,7 ns do neutro. L09 e L15
só entram no MD do complexo se passarem pela etapa de membrana.
