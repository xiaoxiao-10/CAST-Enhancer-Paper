# Stage 12 citation integrity report

Status: **PASS — citation order, identifiers, metadata and claim alignment verified**

Audit date: 2026-08-27
Scope: `paper/manuscript_v1.md`, `paper/references.bib`, `paper/literature_search_report.md`, `paper/manuscript_validation.json` and `paper/reference_verification.json`.

## Mechanical integrity

- The English manuscript contains 40 citation markers covering all 23 retained references.
- First appearances follow Vancouver order 1–23 without gaps; rendered reference-list numbering is also exactly 1–23.
- Every numeric citation maps to the expected BibTeX key, and no cited key or bibliography entry is orphaned.
- All 23 BibTeX entries contain a persistent identifier: 19 DOIs and 4 official archival URLs.
- Current online verification resolved all 23 identifiers. Crossref returned exact DOI and title matches for all 19 DOI records; the four official arXiv, OpenReview and PMLR URLs also resolved. No metadata mismatch or network warning remained in the frozen verification report.
- The online audit sent only bibliography identifiers and titles to public metadata endpoints; unpublished manuscript text was not transmitted.

## Citation-to-claim alignment

| Evidence role | References | Manuscript use | Audit result |
|---|---|---|---|
| CAGE and enhancer transcription | 1–6 | Operational CAGE/eRNA task and independent nascent-transcription context | Aligned; no claim that transcription is the only enhancer definition |
| Accessibility and orthogonal chromatin assays | 7–9 | ATAC-seq interpretation, ENCODE provenance and H3K27ac context | Aligned; no causal assay claim |
| Prior multimodal work | 10–13 | DeepCAPE, CoRE-ATAC, EPCOT and EpiBERT define the novelty boundary | Aligned; the manuscript explicitly rejects a generic first-fusion claim |
| Sequence backbones and baselines | 14–19 | Nucleotide Transformer, DNABERT-2, sequence models, gapped k-mers and LoRA | Aligned with methods and benchmark context |
| Attribution and motif resources | 20–23 | Explains the gated non-use of model attribution and the exploratory JASPAR/MEME motif test | Aligned; zero motif discoveries and no mechanism claim |

No direct quotation is used. Hidden citation anchors are therefore recorded as `none:` rather than inventing page or section locators; DOI and official-source links provide article-level traceability. The strict search report documents 36 screened candidates, 23 retained primary or official sources and the search-bounded contribution boundary. No retraction flag was identified in the frozen screen or current metadata audit, although this is not a substitute for a journal’s final post-publication check.

## Residual pre-submission actions

- Convert the reference list to the exact style and journal abbreviations required by the selected venue.
- Re-run the identifier audit immediately before submission to detect late corrections or retractions.
- Add any new reference only through the same search, metadata and citation-to-claim checks; do not renumber manually without rerunning the validator.

No citation change is required to pass the Stage 12 draft integrity gate.
