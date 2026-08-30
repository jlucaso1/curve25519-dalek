# Continuação: o caminho de grupo usa o backend vetorial, e faz dezoito inversões por mensagem

**Este lote continua o trabalho do PR já aberto neste fork** (`perf(x25519): four field-arithmetic changes on the hot path, plus ADX/BMI2, AVX2 and wasm32 build findings`). Trabalhe **no mesmo branch e no mesmo PR** — não abra outro. Acrescente ao `docs/perf-x25519-field-arithmetic.md` existente em vez de criar um documento novo, e estenda a seção `## Checked and not changed` do corpo do PR conforme o que você medir.

Leia primeiro o que já está lá. Vários resultados deste texto **confirmam** conclusões suas com um segundo caminho de código, e um deles abre uma direção que aquele lote deferiu por falta de motivo.

## O que um segundo perfil mostrou

O perfil que motivou o lote anterior era de troca de mensagens direta (Signal Double Ratchet, X25519 por mensagem). Um segundo perfil, agora de **mensagens de grupo** (Sender Key) no mesmo cliente e na mesma máquina, dá um quadro diferente. Como antes: **estes números vêm de fora deste repositório e são o motivo da tarefa, não a especificação.**

| | DM | grupo (128 membros) |
| --- | ---: | ---: |
| `vector::avx2` (Ed25519 verify) | ~0% | **40,8%** |
| `serial::u64` | ~64% | 15,4% |
| ciclos por mensagem | 423 778 | 358 906 |
| IPC | 3,01 | 2,31 |

**Duas coisas que interessam a você diretamente:**

1. **O backend vetorial é usado, e domina.** A pilha é `verify_signature` → `vartime_double_scalar_mul_basepoint` → `vartime_double_base_mul` → `_impl_mul`, com 18,4% de custo inclusivo. É a confirmação em produção do que você mediu isoladamente (38 491 ns sob `"simd"` contra 47 949 ns sob `"serial"`): o vetorial entrega **onde está conectado**. A ladder de Montgomery simplesmente não o chama; a verificação de assinatura chama.
2. **A inversão de campo deixou de ser marginal.** Você quantificou a inversão em 5,3% das instruções de `mul_clamped` e ~20% de `mul_base_clamped`, e deferiu safegcd por isso. No caminho de grupo, `pow2k` aparece com **12,13% de custo próprio no perfil de tempo**, e a contagem de chamadas explica por quê.

## A contagem que abre a direção

Callgrind sobre o cenário de grupo, chamadas por mensagem:

| chamador | → alvo | por mensagem |
| --- | --- | ---: |
| `FieldElement51::invert` | `pow2k` | **22,00** |
| `curve_models::ProjectivePoint::as_affine` | `invert` | **16,00** |
| `EdwardsPoint::compress` | `invert` | 2,00 |

São **dezoito inversões de campo por mensagem**, dezesseis delas vindas de conversões projetivo→afim.

Uma inversão custa uma cadeia de exponenciação inteira. `batch_invert` existe nesta crate e resolve exatamente esta forma: inverter \\(n\\) elementos com **uma** inversão e \\(3(n-1)\\) multiplicações, em vez de \\(n\\) inversões.

**A pergunta que decide o lote, e que só você pode responder lendo o código:** essas dezesseis `as_affine` são independentes entre si, ou estão em série, cada uma alimentando a próxima? Se forem independentes, o batch se aplica e o ganho é grande. Se forem sequenciais — o resultado de uma decidindo a entrada da seguinte — **não há batch possível e o lote fecha aqui**, o que é um desfecho perfeitamente válido.

Comece por aí, antes de qualquer outra coisa.

## O que você não tem

Você **não tem** o cliente que produziu o perfil, nem o servidor de fixture. Não tente reproduzir os 40,8% nem as dezoito inversões: elas vêm do consumidor, não desta crate. O que você pode fazer aqui dentro:

- Rastrear, no código, quais chamadas de `as_affine` e `compress` um `verify` de Ed25519 faz, e se elas são independentes.
- Estender o harness que você já construiu com um kernel de **verificação de assinatura** (`ed25519-dalek` está no workspace), que hoje ele não cobre — os kernels existentes são X25519 e keygen.
- Medir `batch_invert` contra inversões repetidas nos tamanhos que interessarem, se o rastreio mostrar independência.

## As perguntas, em ordem

1. **As dezesseis `as_affine` são batchable?** Sim/não com a referência de código. Este é o item 1 e pode encerrar o lote.
2. **Se forem, onde o batch caberia** sem mudar API pública nem semântica — dentro do próprio caminho de verificação, ou exigiria que o chamador (`ed25519-dalek`, ou o consumidor) mudasse? Se exigir mudança no chamador, isso muda a natureza do lote e precisa ser dito.
3. **safegcd, agora com um motivo maior.** Você o deferiu com a inversão em 5,3%. No caminho de verificação ela pesa mais. Reavalie **apenas se** o item 1 der negativo — se o batch resolver, safegcd continua sem justificativa, porque uma inversão por mensagem não sustenta a complexidade que você mesmo apontou (a propriedade constant-time dele é global à iteração, não local).
4. **Há headroom no `vartime_double_scalar_mul_basepoint`?** Ele é 18,4% inclusivo do cliente em grupo e você já mostrou que o vetorial ganha 20% ali sobre o serial. Vale medir se as tabelas ou o wNAF têm folga — mas lembre-se de que você **já mediu e rejeitou** tabelas maiores para o caso base (radix-32 é um empate, radix-64 é pior). Verifique se aquela conclusão transfere para o caminho vartime antes de repetir o experimento.
5. **Os 15,4% em `serial::u64` no grupo — de onde vêm?** Se forem majoritariamente as inversões acima, o item 1 já os cobre. Se houver outra fonte, diga qual.

## Quando não fazer

Recuse por item, com números. **Recusar continua sendo desfecho completo**, e neste lote é o resultado mais provável em pelo menos dois itens.

- **Se as `as_affine` forem sequenciais**, feche o item 1 e diga isso. Não force um batch onde não há independência.
- **Se o batch exigir mudar API pública.** Este fork tem um consumidor esperando por um upstream; ampliar superfície pública torna o merge upstream mais difícil, não menos.
- **Se o ganho ficar abaixo do que este host resolve.** Você já estabeleceu que ele é ruidoso e já usou pareamento com contagem exata de instruções para contornar. Use o mesmo padrão, e a mesma honestidade que usou na mudança 4.
- **Se o caminho for do `ed25519-dalek` e não do `curve25519-dalek`.** Diga em qual crate o problema vive. O lote anterior manteve o escopo em `curve25519-dalek` de propósito.

## Restrições

Valem as mesmas do lote anterior, e por isso não vou repeti-las em detalhe: constant-time é critério de morte, nada de `unsafe` injustificado, sem mudança de representação de limbs, sem mudança de API pública, wasm32 e alvos sem SIMD não podem regredir, sem dependência nova.

Uma a mais: **o que já foi commitado no PR não deve ser tocado** a não ser que este trabalho mostre que algo ali está errado — e nesse caso, corrija e diga no PR o que mudou de leitura.

## O que commitar

O mesmo contrato do lote anterior:

- **Acrescente ao `docs/perf-x25519-field-arithmetic.md`** uma seção sobre o caminho de verificação, com máquina, flags, iterações e dispersão, incluindo os itens recusados e o número que os recusou.
- **Código apenas do que der benefício comprovado.** Se o item 1 der negativo, o documento entra e o código não.
- **Os kernels de bench novos**, mesmo os dos itens recusados.

## PR

Não abra outro. Faça push no branch existente e **atualize o corpo do PR**: acrescente as medições novas em `## Measurements` e os itens recusados em `## Checked and not changed`, mantendo os que já estão lá. Se o título deixar de descrever o conteúdo, ajuste-o.

## Testes mínimos

Os mesmos do lote anterior — suíte completa em cada configuração de backend, testes diferenciais para qualquer aritmética nova, vetores RFC 7748 se `montgomery.rs` mudar, build wasm32, `fmt` e `clippy`. Se você tocar em verificação de assinatura, acrescente os vetores conhecidos de Ed25519 (RFC 8032) ao conjunto.
