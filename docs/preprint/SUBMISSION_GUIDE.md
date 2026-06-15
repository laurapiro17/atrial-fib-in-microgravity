# How to publish this — the free path

Two steps, both €0. Do the preprint first (it timestamps your priority), then Cureus.

## Step 0 (strongly recommended first): get a faculty co-author

This is the single highest-leverage move and it is free. A hypothesis paper in
space cardiology is far more credible — and far more likely to be accepted — with
a senior co-author (a UB cardiology or physiology professor, or anyone working in
cardiac electrophysiology / space medicine). It also helps your MD-PhD path
(mentorship, a real reference). You can still be first and corresponding author.
Send them `cureus_manuscript.md` + the three figures and ask if they'll co-sign.

Not strictly required, but do not skip it lightly.

## Step 1 — Preprint (free, instant DOI, fixes priority)

Pick ONE:

- **medRxiv** (recommended — health/clinical home). Account at medrxiv.org →
  submit. No publication fee. You must affirm there are no patient data (true —
  this is in-silico). medRxiv screens before posting (a few days). Category:
  "Cardiovascular Medicine".
- **arXiv** (`q-bio.TO` Tissues and Organs, or `physics.med-ph`). Free, but a
  first-time submitter needs an **endorsement** from an existing arXiv author in
  that category — a real hurdle if you don't know one. If a co-author already uses
  arXiv, this is trivial; otherwise prefer medRxiv.

Use `scaling_law_manuscript.md` (the generic version) or `cureus_manuscript.md`.
Convert to PDF first (see "Producing the file" below).

## Step 2 — Cureus (free, peer-reviewed, PubMed/PMC-indexed)

1. Account at cureus.com. **There is no submission or publication fee** for a
   standard article. (Cureus offers optional paid editing/expedited review — you
   do NOT need it. If anything ever asks for payment to publish, stop and check.)
2. Article type: **Original Article**.
3. Paste from `cureus_manuscript.md`: it already has Cureus's required pieces —
   structured abstract, Introduction / Materials And Methods / Results /
   Discussion / Conclusions, the **Additional Information** disclosures block, and
   numbered references with DOIs.
4. Upload the three figures with captions:
   - Figure 1 → `figures/gravity_law.png` (𝒩_g vs g, with g*)
   - Figure 2 → `figures/gravity_corroboration.png` (6-seed sheet burden + 𝒩_g)
   - Figure 3 → `figures/gravity_sensitivity.png` (g* vs L₀ and map shape)
5. Disclosures are pre-filled: no human/animal subjects, no conflicts, no funding,
   data openly on GitHub.
6. Cureus runs editorial + community (SIQ) review. Expect revision requests —
   that's normal and good.

## Producing the file to upload

The manuscripts are Markdown. To get a PDF/Word for submission, ask me to convert
`cureus_manuscript.md` (or the generic one) — your machine has md→PDF (pandoc +
Chrome) and a docx route. Figures are already PNGs in `figures/`.

## Honest expectations

- This is a **hypothesis paper**: it proposes and tests a falsifiable scaling
  argument in silico. Reviewers will (rightly) push on the `g→stretch` assumption
  and on `g*` being conditional on atrial size. The manuscript already concedes
  these in Limitations — keep that honesty; it is a strength, not a weakness.
- A rejection or a "major revision" is not failure; it is the normal path. A
  co-author and one revision cycle usually get a sound hypothesis paper in.
