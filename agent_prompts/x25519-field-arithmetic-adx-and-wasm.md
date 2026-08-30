# O AVX2 desta crate não alcança o X25519, e o backend que o wasm32 recebe nunca foi medido

Você é um agente autônomo trabalhando num fork de `curve25519-dalek`. Leia o `README.md`, o `CONTRIBUTING.md` e o `build.rs` antes de começar — a seleção de backend acontece lá e é metade deste lote.

Este lote é sobre **aritmética de corpo no caminho do X25519**, em dois alvos que o consumidor deste fork usa: **x86_64 com ADX/BMI2** e **wasm32**. Ele nasce de um perfil externo e a maior parte dele é medição; a edição de código só acontece onde a medição justificar.

## De onde vem a motivação

Perfil de um cliente do protocolo WhatsApp (`oxidezap/whatsapp-rust`) que usa esta crate via `wacore-libsignal`, num cenário de troca de mensagens ping-pong, `perf` em hardware real (não sob valgrind), Zen 4 com `avx2`, `sha_ni`, `adx` e `bmi2`. Números são **o motivo da tarefa, não a especificação**:

| símbolo | self |
| --- | ---: |
| `serial::u64::field::FieldElement51::pow2k` | 30,8% |
| `serial::u64::field::FieldElement51 as Mul>::mul` | 25,3% |
| `RustCryptoProvider::x25519_agreement` | 7,7% |
| `sha2::sha256::x86_sha::compress` | 3,3% |

**Cerca de 64% do tempo do cliente está em aritmética de corpo serial**, com duas operações X25519 por mensagem. O binário **tinha** `backend::vector::avx2` compilado e 531 instruções AVX2, e mesmo assim o caminho quente é serial.

A leitura de código que explica isso, e que você deve confirmar em vez de aceitar:

- `src/field.rs` resolve `FieldElement` sempre para um tipo **serial** (`serial::u64::FieldElement51`, `fiat_u64`, ou `serial::u32::FieldElement2625`). Nenhum ramo do `cfg_if!` resolve para um tipo vetorial.
- `src/montgomery.rs` usa `crate::field::FieldElement`, logo herda o serial.
- `grep -rn montgomery src/backend/vector/` não retorna nada: o backend vetorial implementa `EdwardsPoint`, não a ladder de Montgomery.

Se isso estiver certo, **o AVX2 desta crate acelera ed25519 e não acelera X25519**, e é por isso que um binário com AVX2 gasta 56% em `serial::u64`.

## As duas frentes

Meça as duas. Elas são independentes e podem ter desfechos opostos.

### Frente A — x86_64 com ADX/BMI2

`src/backend/serial/u64/field.rs` faz a multiplicação 64×64→128 com `u128` do Rust (`fn m(x, y) -> u128`), e acumula em `u128`. A CPU alvo tem `adx` (`adcx`/`adox`) e `bmi2` (`mulx`), que existem exatamente para cadeias de multiply-accumulate: duas cadeias de carry independentes, sem tocar nas flags uma da outra.

Implementações em assembly com ADX (BoringSSL, libsodium, fiat-crypto com ADX) tipicamente ganham 20–40% em X25519 sobre aritmética portável. **Isso é referência de terceiros, não uma medição deste repositório.** Se o LLVM já estiver emitindo `mulx`/`adcx` a partir do `u128` quando compilado com `-C target-feature=+adx,+bmi2`, o ganho pode ser zero — e essa é a primeira coisa a verificar.

**Verifique antes de escrever qualquer código**: compile `serial::u64::field` com e sem `-C target-cpu=native`, desmonte `mul` e `pow2k`, e conte `mulx`, `adcx`, `adox`. Se já aparecerem, o compilador já faz o trabalho e a frente A provavelmente morre aqui. Reporte a contagem nos dois casos.

### Frente B — wasm32

`build.rs` escolhe `curve25519_dalek_bits` por `target_pointer_width`. wasm32 tem ponteiro de 32 bits, então recebe `DalekBits::Dalek32`, ou seja **`serial::u32::FieldElement2625`**.

Mas wasm32 **tem `i64` nativo**: `i64.mul`, `i64.add` são instruções de primeira classe, e um `u64` não é emulado como seria num ARM de 32 bits real. O próprio `build.rs` marca isso como pendente:

```
//TODO(Wasm32): Needs tests + benchmarks to back this up
```

**Este lote é a chance de fechar esse TODO com números.** Meça `curve25519_dalek_bits="64"` contra o padrão `"32"` em wasm32. Se o u64 ganhar, o benefício é imediato para qualquer consumidor wasm e não exige uma linha de aritmética nova.

Uma advertência que vem de fora e que você deve levar a sério: **perfis de wasm32 não generalizam a partir de nativo.** Efeitos de capacidade e de layout já foram observados invertendo de sinal entre wasm32 e 64-bit em outros lotes. Meça wasm no wasm; não extrapole do x86_64 em nenhuma direção.

## Você não tem o harness que produziu a motivação

Você **não tem** o cliente WhatsApp, nem o servidor de fixture, nem o benchmark multi-processo. Não tente reproduzir os 64%. Este lote é **somente sobre esta crate**, e toda medição sua tem que nascer aqui dentro.

O que existe: `benches/dalek_benchmarks.rs`, com um grupo `montgomery_benches` que já exercita a ladder (`montgomery.rs`), e `required-features = ["alloc", "rand_core"]`.

O que você precisa criar:

1. **Bench de `FieldElement51::mul` e `pow2k` isolados.** É o alvo direto da frente A, e é a única forma de saber se uma mudança na aritmética fez o que promete antes de olhar a ladder.
2. **Bench de X25519 ponta a ponta na crate**: `MontgomeryPoint::mul_clamped` e `mul_base_clamped`, que é o que o consumidor realmente chama. O ganho aqui é o que importa; o do item 1 é diagnóstico.
3. **Bench comparativo de backends**: o `serial::u64` padrão, `curve25519_dalek_backend="fiat"`, e a sua variante ADX se ela existir. O fiat entra como **controle barato**: é um cfg, zero código novo, e se ele já ganhar, a frente A tem que superar ele, não o padrão.
4. **Bench de wasm32.** Criterion não roda em wasm direto. Use o que couber — `wasm-pack test --node` com um harness de tempo próprio, ou um binário wasm chamado por `wasmtime`/`node` que rode N iterações e reporte tempo. Diga qual escolheu e por quê; o método faz parte do resultado.

Reporte **ciclos ou tempo por operação**, não só razões, e diga o número de iterações e como o ruído foi controlado.

## Meça também estas, mesmo que não vire código

O lote pede o inventário, não só o alvo:

- **Quantas operações de corpo uma `mul_clamped` faz.** A ladder são 255 iterações de `differential_add_and_double` (`montgomery.rs:430`), cada uma com um número fixo de `mul`/`square`. Conte-as e mostre a aritmética que liga "X% em `mul`" a "N `mul` por X25519". Isso transforma o perfil externo em algo verificável aqui.
- **`pow2k` contra `square` repetido.** `pow2k(k)` existe para amortizar; confirme que ele está sendo usado onde deveria e que o ganho dele ainda existe no seu alvo.
- **O custo do backend `fiat`** nos dois alvos. Ele é formalmente verificado, o que tem valor próprio; se empatar em performance, isso é informação para o consumidor.
- **Se o backend vetorial poderia cobrir Montgomery.** Não implemente — apenas diga se é viável e o que custaria, porque essa é a pergunta que a motivação levanta e que ninguém respondeu. Uma frase fundamentada vale mais que um palpite longo.

## Quando não fazer

Recuse por frente, com os números. **Recusar é um desfecho completo e válido, e aqui é provável em pelo menos uma das frentes.**

- **Frente A morre se o LLVM já emite `mulx`/`adcx`/`adox`** com `target-feature=+adx,+bmi2`. Conte as instruções antes de escrever assembly.
- **Frente A morre se o ganho isolado em `mul` não sobreviver na `mul_clamped`.** A ladder tem dependências seriais; um `mul` 30% mais rápido pode render muito menos ponta a ponta. Reporte os dois números lado a lado, e se divergirem, o número que vale é o da `mul_clamped`.
- **Frente B morre se `bits="64"` empatar ou perder em wasm32.** É um resultado publicável e fecha o TODO do `build.rs` com evidência.
- **Qualquer frente morre se quebrar constant-time.** Isto é código criptográfico: nenhum branch, índice de memória ou early-return pode depender de dado secreto. Uma otimização que introduza isso está errada mesmo que seja mais rápida, e não há número que a justifique.
- **Qualquer frente morre se exigir `unsafe` que você não consiga justificar.** `core::arch` intrinsics e `asm!` são `unsafe` por natureza; se usar, isole no menor escopo possível e documente a precondição de cada bloco.

## Restrições

- **wasm32 e alvos sem ADX não podem regredir.** Qualquer caminho novo fica atrás de `cfg(target_feature = ...)` ou equivalente, com o código atual intacto como fallback. Prefira gate em tempo de compilação a detecção em runtime: esta crate é `no_std`-friendly e `is_x86_feature_detected!` exige `std`.
- Não mude a API pública. `FieldElement51` é `pub(crate)`; o que sai da crate são `MontgomeryPoint`, `Scalar`, `EdwardsPoint`.
- Não mude a representação de limbs (5×51 bits) nem as pré-condições documentadas de `bit excess`. Elas são invariantes de correção, não detalhes.
- Não mude o formato serializado nem nada observável de fora.
- Sem dependência nova. `fiat` já está no grafo por trás de um cfg.
- Não toque em `ed25519-dalek` nem em `x25519-dalek` (os outros membros do workspace) a não ser que a medição peça, e nesse caso diga por quê.

## O que commitar

Este é o ponto onde este lote difere dos habituais.

1. **Um `.md` com o progresso e os resultados de tudo que rodou**, incluindo as frentes que não renderam. Coloque em `agent_prompts/results/` ou em `docs/`, escolha e seja consistente. Ele precisa conter:
   - a máquina: CPU, features (`lscpu`), versão do rustc, flags de compilação de cada medição;
   - cada bench com número, unidade, número de iterações e dispersão;
   - as frentes recusadas, com o número que as recusou;
   - a contagem de `mulx`/`adcx`/`adox` antes e depois, se a frente A foi adiante;
   - o resultado wasm32 com o método usado.
2. **Código apenas do que deu benefício comprovado.** Se a frente A não render, o `.md` entra e o código não. Se só a frente B render, commite só o ajuste de `build.rs` (ou a documentação do cfg) e o bench que provou.
3. **Os benches sintéticos que você criou**, mesmo os das frentes recusadas — eles são o que permite refazer a medição depois.

Não commite código de uma frente recusada "para referência". Ele vira dívida que alguém vai reativar sem os números.

## PR

- Branch: `perf/x25519-field-arithmetic` a partir do default deste fork.
- Título: `perf(field): ADX/BMI2 path and wasm32 backend measurement for X25519`
  - Se só a frente B render: `perf(wasm32): use the 64-bit field backend`
  - Se nenhuma render: `docs(perf): measure X25519 field arithmetic on x86_64 and wasm32`
- Corpo: `## Summary` (as duas frentes e o que a motivação externa dizia) · `## Changes` (vazio é válido) · `## Measurements` (a tabela principal, com a máquina) · `## Checked and not changed` (as frentes recusadas com seus números; a viabilidade do vetorial em Montgomery) · `## Validation` (suíte, constant-time, build wasm32).
- Aponte o `.md` completo a partir do corpo do PR em vez de colar tudo nele.

## Testes mínimos

- A suíte da crate inteira, nos dois caminhos: com e sem o cfg do caminho novo. Um caminho que só é testado quando ativado não está testado.
- Um teste diferencial: para entradas aleatórias, a `mul` e a `square` novas produzem bit a bit o mesmo que as atuais, incluindo limbs no limite superior do `bit excess` documentado.
- Build wasm32 (`cargo build --target wasm32-unknown-unknown`) provando que o gate segura.
- `cargo fmt --all`, `cargo clippy --all-targets -- -D warnings`.

## Pós-PR (obrigatório — não encerre antes)

1. Rode a suíte e o clippy uma última vez no HEAD do PR.
2. Se o CI do fork estiver ativo, `gh pr checks <numero> --watch` sem `--fail-fast`; confirme antes se a falha já existe no branch base.
3. Não faça merge e não adicione labels.
