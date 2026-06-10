# DRAFT — "How arrays are called" concepts text (for README)

Status: draft for review. Three alternative write-ups of how TideCluster
describes each tandem-repeat array (founder / HOR / SSR), written for an
**informed user** (knows what a tandem repeat / satellite / microsatellite is,
but not the internal pipeline). Pick one to fold into the README "Concepts"
section; they can also be mixed (e.g. Option A intro + Option B table).

A worked example used by all three (for consistency if we add one):
*founder 53 bp, strongest 159 bp → HOR ×3 (three 53-bp monomers per
higher-order unit).*

---

## Option A — Plain-language narrative (analogy-led)

A tandem repeat is a stretch of DNA where one short unit is copied over and
over, like a word repeated down a line. TideCluster finds these arrays and then
asks a simple question of each one: *what is the basic repeating unit, and is
there a larger pattern built on top of it?*

- **Founder** — the basic building block: the shortest unit the array is
  genuinely made of. Think of it as the single "word" being repeated.
- **Strongest** — the unit that the sequence matches *most cleanly*. Often this
  is the founder itself, but sometimes the cleanest match is a longer block made
  of several founders in a row.
- **HOR (higher-order repeat)** — when that cleanest block is a whole number of
  founders (2, 3, … copies), the array has a higher-order pattern: the words are
  grouped into repeating "sentences". The **multiplicity** (×k) is how many
  founder copies make up one sentence.
- **SSR (simple sequence repeat / microsatellite)** — some arrays are just a
  very short, simple motif repeated (e.g. `AT`, `ATC`). These are recognised
  earlier, by their consensus, and reported with the motif and its founder
  length directly — no higher-order analysis needed.

For every HOR call we also say *how confident* we are in the ×k pattern: a
clean, low-multiple repeat is reported as a solid HOR; a very long or irregular
one is reported as "founder recovered" — we trust the building block, but treat
the exact ×k as approximate.

---

## Option B — Glossary / definitional (crisp, table-friendly)  ← chosen

**Terms used to describe a tandem-repeat array (TRA)**

| Term | In simple words |
|---|---|
| **Founder** | The basic repeat unit — the shortest building block the array is actually made of (its monomer). |
| **Strongest** | The period the sequence matches most precisely. Equals the founder for a simple array; a multiple of it when there's higher-order structure. |
| **Multiplicity (×k)** | How many founder copies make up the strongest unit. ×1 = no higher-order structure. |
| **HOR (higher-order repeat)** | An array where the strongest unit is k≥2 founders joined together — whole monomers organised into a larger repeating block. |
| **Subrepeat** | A short tandem repeat *inside* the founder monomer (period much shorter than the founder — only part of one monomer). The opposite direction to HOR: HOR groups whole monomers into a bigger unit, a subrepeat is a smaller pattern within a single monomer. |
| **SSR (simple sequence repeat)** | Just a tandem repeat with a very short monomer — a simple short motif repeated (e.g. `AT`, `ATC`). Recognised from its consensus; the motif and its length are reported directly. |
| **HOR-order confidence** | How trustworthy the ×k call is: **strict/supported** = a confident higher-order pattern; **weak** = the founder is reliable but the exact ×k is only approximate (very long, irregular, or rescued). |

**In short:** every array gets a *founder* (its monomer). If the cleanest match
is several whole founders long, it's an *HOR* (×k). If there's a shorter pattern
*within* the monomer, that's a *subrepeat*. An array whose monomer is itself a
very short simple motif is an *SSR*.

---

## Evidence behind the calls — three companions to pair with Option B

A short paragraph to sit under the Option B table, explaining (simplified, for a
genomics/biology audience) *what information the founder/HOR calls are based on*.
Pick one.

### Companion 1 — prose

These calls come from two complementary signals that TideCluster (via `kitehor`)
measures along each array. The first is **periodicity**: the array is scanned for
the spacing at which its sequence repeats, yielding a set of candidate monomer
sizes, each with a strength score (how dominant that period is). The second is an
**identity scan**: for each candidate period, copies one period apart are
compared and their median sequence identity is recorded — a genuine repeat unit
has well-conserved copies (high identity), whereas a spurious period does not.
The **founder** is the shortest period that is both a real divisor of the array
and well-conserved; the **strongest** is the best-conserved period overall. When
the strongest is a whole-number multiple of the founder, the array is an **HOR**.

### Companion 2 — "two signals" bullets

Founder/HOR calls combine two measurements taken along each array:

- **Periodicity (period + strength)** — the spacing at which the sequence
  repeats. A scan returns candidate monomer sizes, each scored by how strong its
  periodic signal is. This proposes *which* repeat units might be present.
- **Identity scan (conservation)** — for a candidate period, how similar the
  repeated copies are to each other (median pairwise identity). This tells us
  *which* of those candidates are real, well-conserved units rather than
  by-chance signal (short, low-complexity periods can look periodic but score low
  identity).

The **founder** is the shortest well-conserved unit; the **strongest** is the
most conserved; a clean integer ratio between them is an **HOR**.

### Companion 3 — mini-table (matches B's table style)

**What each call is based on**

| Signal | What it measures | What it decides |
|---|---|---|
| **Periodicity** | The spacing(s) at which the array repeats, with a strength score per candidate monomer size. | Proposes the candidate repeat units. |
| **Sequence identity** | How well copies one period apart match each other (median pairwise identity). | Confirms which candidates are genuine, well-conserved units. |
| **Founder ÷ Strongest ratio** | Whether the best-conserved period is a whole-number multiple of the basic unit. | An integer ratio (k≥2) → **HOR ×k**; how clean/regular it is → the confidence tier. |

---

## Option C — Decision walkthrough (how the pipeline decides)

For each array, TideCluster works through three questions:

1. **Is it a simple sequence (SSR)?** If the array's consensus is overwhelmingly
   one short, simple motif (e.g. `AT`, `ATC`), it's called an **SSR**
   (microsatellite). We report the motif and its length as the founder, and stop
   here — no higher-order analysis.

2. **What is the basic unit (founder)?** Otherwise the array is treated as a
   satellite-type repeat, and we find its **founder** — the shortest unit it is
   genuinely built from. This is the number to trust most.

3. **Is there a higher-order pattern (HOR)?** We then check the unit the sequence
   matches most cleanly (the **strongest**). If that unit is a whole number of
   founders — k≥2 copies — the array is a **higher-order repeat (HOR)** with
   **multiplicity ×k**: the monomers are organised into a larger repeating block.
   If the strongest unit *is* the founder (×1), there's no higher-order structure.

Finally, every HOR is tagged with a **confidence**: a clean, low ×k is a
confident HOR; a very long, irregular, or rescued one is reported as **"founder
recovered"** — we keep the building block but flag the exact ×k as approximate.

---

## Open questions for review

1. Which option (or which combination) goes into the README?
2. Include the worked example (`founder 53 → strongest 159 → HOR ×3`) inline, or
   leave the text purely conceptual?
3. Keep the `strict / supported / weak` tier names visible here, or describe the
   confidence in words only and let the report legend carry the tier names?
