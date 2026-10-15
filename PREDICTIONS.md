# Previsões a congelar antes dos ensaios

**Estado: rascunho, ainda não congelado (sem tag).** Para congelar: revisar este arquivo, fazer o commit e criar a tag antes de enviar qualquer composto ou resultado a um grupo experimental.

```bash
git tag -a predictions-v1 -m "Previsões congeladas antes dos ensaios"
git push origin predictions-v1
```

Depois da tag, alterações vão em um arquivo novo ou em uma versão nova (`predictions-v2`), nunca por edição silenciosa deste texto.

## Previsões

| Medida | Previsão | O que contaria contra |
|---|---|---|
| pKa do lapachol | Entre 5 e 6; o lawsone, mais ácido que o lapachol | Lapachol fora de 4,5 a 6,5: a classificação por estado de ionização precisa ser refeita |
| logD em pH 7,4 | Lawsone < lapachol < β-lapachona < 2-metoxi-lapachol (valores calculados: −1,9; 1,3; 2,8; 3,3) | Lapachol acima da β-lapachona: o logP de Crippen com pKa por classe não serve para esta família |
| Partição em POPC | Sem previsão por enquanto: o perfil ΔG(z) do lapachol neutro ainda não convergiu | Será acrescentada e datada antes de qualquer medida de partição |
| Potencial redox | Sem previsão de ordem: a diferença calculada de afinidade eletrônica entre os quatro compostos é de até 0,16 eV, do tamanho do erro do xTB (0,09 eV no lapachol) | A voltametria entra como validação do xTB, não como teste |
| Cinética na NQO1 | Sem previsão: o docking não separa os compostos e não diz se algum é substrato | A β-lapachona serve de controle positivo |

Os valores de logD vêm de um logP aditivo (Crippen) corrigido por um pKa heurístico, então são uma ordem de grandeza e não uma previsão fina. A fonte dos números é `funil_consolidado.csv`.

## Conjunto mínimo de compostos e medidas

Compostos: lapachol, lawsone, β-lapachona e 2-metoxi-lapachol (pureza acima de 95% por HPLC). Medidas, em ordem de prioridade: pKa por titulação espectrofotométrica; logD em pH 7,4; partição em lipossomos de POPC em pH ácido e em pH 7,4. Opcional: cinética de consumo de NADH com NQO1 recombinante, com dicumarol e β-lapachona como controles.
