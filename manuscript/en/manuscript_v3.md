# Nucleotide composition accounts for the shared component of enhancer activity across unseen cell types

## Abstract

Enhancer activity is cell-type specific, yet cross-cell generalisation of
sequence models is usually evaluated by holding out cell types while sharing all
candidate loci. We combined a cell-type holdout with a locus holdout: 20% of
FANTOM5 candidate enhancers were excluded from all development-cell training,
three model families were retrained on the remaining 80%, and performance was
evaluated in 11 primary cell types never used in training. A 21-parameter
logistic nucleotide-composition model reached a 5.66-fold lift over prevalence
on unseen loci; a 307,265-parameter from-scratch convolutional network did not
exceed it (5.06-fold). A 50-million-parameter pretrained Nucleotide Transformer
with low-rank adaptation reached 8.35-fold and exceeded the logistic
composition model. This comparison depended on the composition readout: an
XGBoost model trained on the same 20 composition features reached 0.243 AUPRC
on unseen development-active loci, compared with 0.161 for logistic composition
and 0.241 for the transformer; the transformer was indistinguishable from
composition-XGBoost (−0.0011; 95% CI −0.0358 to +0.0365; better in 6 of 11 cell
types). Adding 879 JASPAR motif counts did not improve on composition-XGBoost.
The transformer retained an advantage only on retained development-active loci
(+0.0389; 95% CI 0.0155 to 0.0614), consistent with a locus-specific rather than
a transferable signal, and no model improved on composition-XGBoost for loci
inactive in all development cells. Nucleotide composition therefore accounts
for the shared component of enhancer activity across unseen cell types;
comparisons with linear composition alone overstate the transferable advantage
of pretrained sequence models.

**Keywords:** enhancer RNA; CAGE; cross-cell generalisation; nucleotide
composition; pretrained sequence models; locus holdout

## 1. Introduction

Enhancers are regulatory DNA elements whose activity depends on the sequence of
the locus and on the state of the cell [1]<!--ref:Yang2024--><!--anchor:none:-->;
super-enhancers are composite domains that combine classical enhancers and
facilitators [2]<!--ref:Blayney2023--><!--anchor:none:-->.
Bidirectional transcription initiation
is a recurrent feature of active enhancers, and cap analysis of gene expression
(CAGE) has made enhancer transcription measurable across large collections of
human cell types and tissues [3]<!--ref:Andersson2014--><!--anchor:none:--> [4]<!--ref:Kawaji2017--><!--anchor:none:--> [5]<!--ref:Forrest2014--><!--anchor:none:-->.
Nascent-transcription assays, including PRO-seq [8]<!--ref:Mimoso2023--><!--anchor:none:-->,
show that initiation at enhancers and promoters shares architectural features
[6]<!--ref:Narita2021--><!--anchor:none:--> [7]<!--ref:Fong2022--><!--anchor:none:-->.
CAGE-defined enhancer activity is therefore an operational readout of regulatory
activity, without implying that transcription is the only definition of enhancer
function [9]<!--ref:Han2022--><!--anchor:none:-->.

Chromatin accessibility supplies a complementary view; ATAC-seq protocols
[10]<!--ref:Grandi2022--><!--anchor:none:--> and the ENCODE catalogue
[11]<!--ref:ENCODE2020--><!--anchor:none:--> provide the standard measurements
that are paired with sequence models;
H3K27 acetylation marks a further active state [12]<!--ref:Kang2021--><!--anchor:none:-->.
Sequence-based prediction has been pursued with convolutional models
[13]<!--ref:Penzar2023--><!--anchor:none:--> [14]<!--ref:Liu2022--><!--anchor:none:-->,
feature-based and k-mer approaches [15]<!--ref:Smith2023--><!--anchor:none:-->, long-context
architectures [16]<!--ref:Avsec2021--><!--anchor:none:--> [17]<!--ref:Avsec2021BPNet--><!--anchor:none:-->,
and sequence-to-function models that predict regulatory activity genome-wide
[18]<!--ref:Chen2022Sei--><!--anchor:none:--> [19]<!--ref:Yuan2022--><!--anchor:none:-->.
Reviews of the field [20]<!--ref:Yue2023--><!--anchor:none:--> [21]<!--ref:Routhier2022--><!--anchor:none:-->
and of genomic language models [22]<!--ref:Benegas2025--><!--anchor:none:-->
summarise these developments. Multimodal studies
have combined sequence with accessibility or related epigenomic measurements
[23]<!--ref:Chen2021--><!--anchor:none:--> [24]<!--ref:Thibodeau2021--><!--anchor:none:--> [25]<!--ref:Zhang2023--><!--anchor:none:--> [26]<!--ref:Javed2025--><!--anchor:none:-->.
Pretrained genomic language models provide transferable sequence representations
[27]<!--ref:DallaTorre2025--><!--anchor:none:--> [28]<!--ref:Ji2024--><!--anchor:none:-->
that can be adapted parameter-efficiently [29]<!--ref:Hu2022--><!--anchor:none:-->.

Two features of the prevailing evaluation practice leave the central question
open. First, cross-cell generalisation is usually tested by holding out cell
types while keeping every candidate locus in training. A model with enough
capacity can then recognise a locus seen during development rather than learn a
sequence rule that transfers to new loci [30]<!--ref:Whalen2022--><!--anchor:none:-->.
Second, sequence models are usually
compared with other sequence models rather than with an explicit
nucleotide-composition baseline, even though composition is itself a sequence
property and a strong candidate explanation for any transferable signal.

Here we ask two questions together. Does a pretrained sequence model outperform
an explicit composition baseline in cell types that were never used in training?
And does any such advantage survive when the candidate loci themselves are held
out, and does it reflect a shared or a cell-specific component of enhancer
activity? We address them with an 11-cell-type held-out design in which 20% of
FANTOM5 candidate loci were excluded from all development-cell training, and with
three model families spanning more than five orders of magnitude in parameter
count: a 21-parameter composition model, a 307,265-parameter from-scratch
convolutional network, and a 50-million-parameter pretrained Nucleotide
Transformer with rank-8 low-rank adaptation.

Composition alone transferred strongly to unseen loci (5.66-fold over
prevalence), and the from-scratch network did not exceed it (5.06-fold). The
pretrained model exceeded a logistic composition model (+0.0175 AUPRC over all
unseen loci; +0.0810 within development-active loci), but the same 20
composition features modelled nonlinearly with XGBoost reached a higher AUPRC
than the pretrained model on unseen development-active loci (0.243 versus
0.241). The difference between the transformer and composition-XGBoost was
−0.0011 (95% confidence interval −0.0358 to +0.0365), and the transformer was
better in only 6 of 11 cell types. Motif counts did not improve on
composition-XGBoost, and no model improved on the strongest composition
baseline for loci that were never active in any development cell. The
transferable signal is therefore concentrated in a shared activity component
captured by nucleotide composition; the apparent advantage of the pretrained
model over the logistic baseline does not survive a nonlinear readout of the
same composition features.

## 2. Results

### 2.1 A locus-held-out, cell-held-out design

The candidate universe was the FANTOM5 hg38 enhancer atlas (63,285 loci). We
held out 20% of loci (12,360 of 61,752 quality-passing loci), stratified by
chromosome and GC content, and excluded them from every development-cell
training, validation and test partition. The remaining 80% contributed 108,252
training, 12,621 validation and 24,247 test records across K562, HepG2 and
GM12878 (Table 1). All three model families were retrained on this reduced set
(Fig. 1). The holdout removes loci, not their sequence neighborhoods: 29.2% of
candidate loci share at least part of their 2 kb window with another candidate
locus (Methods 4.12).

Evaluation used 11 primary cell types that were never used in training; their
labels were read only after all models had been frozen. Each held-out cell type
contributed 58,045 to 61,891 label-eligible loci and 191 to 1,742 positives
(prevalence 0.31% to 2.90%). For each cell type we report performance separately
on the 20% unseen loci and on the 80% retained loci, and, within each, on loci
that were active in at least one development cell ("development-active") or in
none ("development-inactive").

![Figure 1](figures_v3/fig1_design.svg)

### 2.2 Nucleotide composition is a strong transferable baseline

The 21-parameter composition model (four mononucleotide and sixteen dinucleotide
frequencies, logistic regression, no cell-identity input) reached a 5.66-fold
lift over prevalence on unseen loci and 6.21-fold on retained loci. Within
development-active loci it reached 3.96-fold (unseen) and 5.70-fold (retained);
within development-inactive loci it reached 2.74-fold and 2.70-fold.

Composition is therefore a strong transferable signal in its own right, on
loci that the model has never seen. It is a necessary reference; the nonlinear
composition readout below is the stronger one.

### 2.3 A from-scratch convolutional network does not exceed composition

The 307,265-parameter DNA convolutional network with FiLM cell conditioning
(five-seed ensemble) reached 5.06-fold on unseen loci and 5.61-fold on retained
loci. Within development-active loci it reached 3.69-fold and 4.96-fold; within
development-inactive loci, 2.55-fold and 2.46-fold. At no stratum did it exceed
the composition baseline.

### 2.4 A pretrained model exceeds the logistic composition baseline on unseen loci

The pretrained Nucleotide Transformer v2 50M with rank-8 low-rank adaptation
(five-seed ensemble) reached 8.35-fold on unseen loci and 9.49-fold on retained
loci. The paired gain over the logistic composition model was +0.0175 AUPRC on unseen loci (95%
confidence interval 0.0101 to 0.0291) and +0.0487 on retained
development-active loci (0.0203 to 0.0747; Table 2).

Within unseen development-active loci the gain was +0.0810 (0.0501 to 0.1219),
and the per-cell point estimate was positive in all 11 held-out cell types
(range +0.0280 to +0.1719); per-cell 95% intervals excluded zero in 6 of 11.
The per-seed NT−composition difference ranged from +0.078 to +0.095, so the
five-seed ensemble did not create the gain. The NT model therefore adds
information beyond a linear readout of nucleotide composition, and this
contrast remains positive on loci that were never seen during training. The
next sections test whether the gain survives a nonlinear readout of the same
composition features.

### 2.5 The logistic-baseline advantage is present on unseen loci

Relative to the logistic composition baseline, the gain was present on unseen
development-active loci (+0.0810; 95% CI 0.0501 to 0.1219), so the contrast is
not limited to the loci used in training. The unseen-minus-retained difference
was +0.032 (95% CI −0.004 to +0.077; bootstrap p = 0.074), which does not
resolve whether the gain is larger on unseen loci. This contrast alone
therefore does not establish that the gain is free of locus-level memory; the
nonlinear baseline below provides the stronger test.

![Figure 2](figures_v3/fig2_four_strata.svg)

### 2.6 A nonlinear composition baseline removes the NT advantage over logistic composition

The comparison above used a logistic composition model. To test whether the NT
gain reflected sequence information beyond composition, we trained XGBoost on
the same 20 composition features (composition-XGBoost), a motif-only model on
879 JASPAR 2024 vertebrate motifs [36]<!--ref:JASPAR2024--><!--anchor:none:-->
(motif-XGBoost), and a combined model on
composition plus motif counts (composition+motif-XGBoost).

On unseen development-active loci, composition-XGBoost reached 0.243 AUPRC
(Table 3), compared with 0.161 for logistic composition and 0.241 for NT. The
NT was indistinguishable from composition-XGBoost (−0.0011; 95% CI −0.0358 to
+0.0365), and was better in 6 of 11 cell types. On all unseen loci, the NT
retained a small advantage over composition-XGBoost (+0.0106; 0.0030 to
0.0216), but this did not extend to development-active loci. Motif-only was
lower than composition in every stratum (0.103 versus 0.161 on unseen
development-active), and adding motif counts to composition-XGBoost did not
improve it (0.212 versus 0.243). On retained development-active loci, the NT
exceeded composition-XGBoost (+0.0389; 0.0155 to 0.0614), consistent with a
locus-specific component that does not transfer to unseen loci. On unseen
development-inactive loci, neither the NT nor the motif model improved on
composition-XGBoost (NT − composition-XGBoost −0.0020, 95% CI −0.0138 to
+0.0007), although composition-XGBoost itself was slightly above logistic
composition (0.0117 versus 0.0103). This is the central comparison: the
apparent NT advantage over logistic composition does not survive a nonlinear
readout of the same 20 composition features, and the advantage that remains on
retained loci is consistent with locus-specific signal rather than transferable
sequence grammar.

### 2.7 No advantage is detected in the cell-specific component

The development-active and development-inactive subsets behaved differently. On
unseen development-inactive loci, where activity in the held-out cell has no
counterpart in any development cell, we did not detect an improvement of the NT
model over logistic composition (−0.0007; 95% confidence interval −0.0073 to 0.0013) and reached
2.22-fold versus 2.74-fold for logistic composition. Both models remained above chance.
A nonlinear composition model showed the same pattern: composition-XGBoost
reached 0.0117 AUPRC versus 0.0103 for logistic composition, and the NT did not
exceed composition-XGBoost (−0.0020; 95% CI −0.0138 to +0.0007).

The transferable sequence signal is therefore consistent with a shared activity
component: it predicts activity at loci that are active in at least one
development cell, and a nonlinear composition model ranks that component as well
as the NT. For the
cell-specific component, represented by loci that are inactive in all
development cells but active in the held-out cell, neither model showed a
detectable improvement over the strongest composition baseline.

![Figure 3](figures_v3/fig3_per_cell.svg)

### 2.8 Chromosome-level splits are low-powered and unstable

As a sensitivity analysis we split development-active loci by chromosome
provenance. On retained loci the gain was +0.0562 on training chromosomes and
+0.0151 on held-out chromosomes; on unseen loci it was +0.0835 and +0.1307,
respectively (the held-out-chromosome subset contained about 14 positives per
cell in the unseen stratum and about 46 per cell in the retained stratum). The
magnitude varied substantially across the locus subsets, and the
held-out-chromosome subsets are too small to support a primary claim.
Chromosome-level holdouts are reported as a sensitivity analysis; the
locus-level holdout is the primary generalisation test. The chromosome-level
results are shown in Supplementary Fig. S1.

### 2.9 In the logistic comparison, the unseen-locus contrast requires pretrained initialization under this budget

To test whether the NT gain reflected pretrained weights rather than
architecture or capacity, we trained the same NT configuration from random
initialisation with the same tokenizer, LoRA modules, data, optimizer and
four-epoch schedule (five seeds). On retained development-active loci, the
random-init ensemble reached 6.75-fold (AUPRC 0.258), close to the 6.77-fold
(AUPRC 0.273) of the pretrained ensemble (Δ +0.0153; 95% CI −0.0147 to
+0.0465). On unseen development-active loci, it fell to 2.82-fold (AUPRC
0.114), compared with 5.94-fold (AUPRC 0.241) for the pretrained ensemble and
3.96-fold (AUPRC 0.161) for logistic composition. The pretrained-minus-random-init
difference was +0.1275 (95% CI +0.0917 to +0.1812) and was positive in 11 of
11 cell types; the random-init model did not improve on logistic composition
(−0.0465; 95% CI −0.0723 to −0.0303). The retained-to-unseen dissociation is
consistent with locus-level memorization in the random-init model and a
pretraining-dependent contrast against logistic composition under this training
budget (Fig. 4). It does not establish an advantage over composition-XGBoost.

A cell-free control, in which the cell index was set to UNK for every training
record (five seeds), gave a similar gain over logistic composition on unseen
development-active loci (+0.0842; 95% CI 0.0551 to 0.1243), so the gain over the
logistic baseline is not explained by cell conditioning. This control does not
change the comparison with composition-XGBoost, against which no gain was
detected (Supplementary Table S8).

![Figure 4](figures_v3/fig4_control.svg)

### 2.10 The logistic-baseline contrast survives a sequence-disjoint sensitivity analysis

Because the locus holdout removes annotated loci rather than their sequence
neighborhoods, we repeated the comparison after removing held-out loci whose
2 kb windows overlapped a retained window (Methods 4.12). Of the 12,360
held-out loci, 3,028 (24.5%) overlapped a retained window, and 445 of 1,142
unseen positives (39.0%) did so. Removing every overlapping held-out locus
attenuated the development-active gain over logistic composition from +0.0810 to +0.0667 (95% CI 0.0208
to 0.1131); removing loci with at least 500 bp of overlap gave +0.0543 (0.0065
to 0.0997). The gain over logistic composition on all unseen loci was essentially unchanged (+0.0184,
0.0092 to 0.0324; with the 500 bp rule, +0.0175, 0.0084 to 0.0317), and no gain
was detected on development-inactive loci under either rule (+0.0002, −0.0033
to 0.0021 for the 500 bp rule; 7 cell types). The NT advantage over logistic
composition therefore persists when shared sequence neighborhoods are removed,
but up to about one third of the development-active effect is attributable to
loci whose input windows overlap training sequence (Supplementary Fig. S3 and
Table S6). The NT advantage over logistic composition therefore persists when
shared sequence neighborhoods are removed, but this sensitivity analysis
addresses the overlap caveat only for the logistic comparison. It does not
change the comparison with composition-XGBoost, which remains the primary test
of transferable sequence information.

## 3. Discussion

The central result is a baseline effect. Against logistic composition, the
pretrained Nucleotide Transformer looked clearly better: +0.0175 AUPRC over all
unseen loci and +0.0810 within development-active loci. That comparison,
however, used a linear readout of 20 composition features. When the same
features were modelled with XGBoost, composition-XGBoost reached 0.243 AUPRC on
unseen development-active loci, compared with 0.161 for logistic composition and
0.241 for the NT. The NT was indistinguishable from composition-XGBoost
(−0.0011; 95% CI −0.0358 to +0.0365; better in 6 of 11 cell types). The
apparent advantage of the pretrained model over logistic composition therefore
reflects the linearity of that baseline, not a demonstrated advantage in
transferable sequence grammar.

The shared-versus-cell-specific decomposition locates the signal. On loci that
are active in at least one development cell, a flexible function of base and
dinucleotide frequencies ranks held-out-cell activity as well as the pretrained
transformer. On loci that are inactive in every development cell, neither the
NT nor composition-XGBoost improved on the other, and both stayed close to the
low prevalence of that stratum (NT − composition-XGBoost −0.0020; 95% CI
−0.0138 to +0.0007). The candidate universe contains a large shared component:
among sequence-eligible loci, 49% of held-out-cell positives are active in at
least one of the three development cells when pooled across cell types
(per-cell mean 56%). The remaining cell-specific component is precisely the set
of loci that are inactive in all development cells but active in the held-out
cell, and no cell-agnostic sequence model in this study captured it.

Motif counts did not change this picture. A motif-only model was lower than
composition in every stratum, and adding 879 JASPAR motifs to composition-XGBoost
did not improve it (0.212 versus 0.243 on unseen development-active loci). The
shared component is therefore composition-like rather than motif-vocabulary-like
at the resolution tested here: a flexible function of short nucleotide
frequencies captures it, and adding motif counts does not add measurable
transferable information. This is a statement about predictive sufficiency in
this benchmark, not a claim that transcription-factor syntax is biologically
irrelevant.

The retained-locus result points in the same direction. The NT retained an
advantage on retained development-active loci (+0.0389; 95% CI 0.0155 to
0.0614) but not on unseen loci. That pattern is consistent with a locus-specific
component learned during training, rather than a generalisable sequence rule.
The unseen-minus-retained difference was +0.032 (95% CI −0.004 to +0.077), so a
partial contribution of locus-level memory cannot be excluded. Removing held-out
loci whose 2 kb windows overlapped training sequence attenuated the gain over
logistic composition from +0.0810 to +0.0667 (any overlap) and +0.0543 (at least
500 bp), but this sensitivity analysis does not create an advantage over
composition-XGBoost.

The random-init control sharpens the interpretation of the logistic comparison.
On retained development-active loci, the random-init model was close to the
pretrained model (AUPRC 0.258 versus 0.273; Δ +0.015, 95% CI −0.015 to +0.046),
but on unseen loci it fell to 0.114 versus 0.241 (Δ +0.127, 95% CI +0.092 to
+0.181) and was below logistic composition. Under this training budget,
pretrained initialization matters relative to a random-init NT, but it does not
establish an advantage over a nonlinear composition model. The cell-free
control gives the same lesson: a gain over logistic composition persists when
cell identity is removed (+0.0842; 0.0551 to 0.1243), but no gain is detected
against composition-XGBoost (Supplementary Table S8).

The practical implication is methodological. A 21-parameter logistic model
required no accelerator and no cell identity and reached 5.66-fold over
prevalence on unseen loci; a nonlinear model on the same 20 features reached
0.243 AUPRC on unseen development-active loci. Nucleotide composition should
therefore be reported as a standard baseline, and a nonlinear readout of
composition should be included before a gain is attributed to sequence grammar.
Sequence models evaluated only against other sequence models, or only against a
linear composition model, may attribute information captured by base and
dinucleotide frequencies to regulatory syntax
[20]<!--ref:Yue2023--><!--anchor:none:--> [30]<!--ref:Whalen2022--><!--anchor:none:-->.

The result also helps explain why cross-cell benchmarks that hold out cell types
but share all candidate loci may overstate generalisation. In such designs,
models that learn the shared activity component are evaluated largely on loci
whose activity they have already seen in another cell type. Holding out loci
separates the shared component from the cell-specific component and shows that
the advantage over logistic composition is confined to the former, where a
nonlinear composition model matches the NT.

For the cell-specific component, no improvement over the strongest composition
baseline was detected. Models that receive cell state (for example, matched
chromatin accessibility) are the natural next step, but that comparison lies
outside the sequence-only design used here. The result should not be read as
evidence that sequence models cannot learn cell-specific regulatory syntax; it
shows that the present cross-cell, cell-agnostic evaluation does not demonstrate
such learning beyond a nonlinear composition baseline.

Several limitations constrain these conclusions. The study uses one pretrained
backbone, one label resource and one assay; the candidate universe is limited to
FANTOM5 loci; windows are fixed at 2 kb; and activity is binary. The
development-inactive subset contains about 54 positives per cell, so its
interval, although it excludes positive effects above 0.0013, is estimated from
a limited number of positives. No experimental validation of predicted loci was
performed, and no attribution analysis was used to identify the sequence
features driving NT scores. These are scope constraints; the locus holdout is
not sequence-disjoint (Methods 4.12), so the comparison is paired but “unseen”
means unseen annotated loci rather than unseen sequence neighbourhoods.

## 4. Methods

### 4.1 Data resources

FANTOM5 hg38 enhancer atlas: `F5.hg38.enhancers.bed.gz` (63,285 loci),
`F5.hg38.enhancers.expression.usage.matrix.gz` (1,829 libraries) and
`F5.hg38.enhancers.expression.tpm.matrix.gz`. Reference: hg38 primary assembly
FASTA and the ENCODE hg38 blacklist v2 (`ENCFF356LFX`). Inputs were frozen with
recorded MD5 and SHA256 (`data/raw/metadata.tsv`, `data/raw/checksums.*`).
Accessibility data are not an input to any model reported here.

### 4.2 Candidate loci and windows

All FANTOM5 hg38 loci were candidates; no accessibility filter was applied.
Exclusions were: contig outside chr1–chr22; window extending beyond the contig;
overlap with an ENCODE blacklist interval; N fraction above 1%. Each retained
locus contributed the 2,000 bp window `[mid − 1000, mid + 1000)` with
`mid = floor((start + end) / 2)`, and the SHA256 of the window was recorded.

### 4.3 Activity labels

For a cell state with N usable CAGE libraries, a locus was labelled positive if
at least ⌈2N/3⌉ libraries had usage = 1, negative if all N had usage = 0, and
excluded otherwise. Development cell lines have N = 3, giving the 2-of-3 rule;
held-out primary cell types have N = 4 to 13, giving thresholds of 3 to 9
libraries. Labelling used only the usage matrix.

### 4.4 Development cells and held-out cells

Development cells (labels exposed): K562 (CNhs12334–CNhs12336), HepG2
(CNhs12328–CNhs12330) and GM12878 (CNhs12331–CNhs12333). Held-out cells (labels
never read during training or selection): 11 primary cell types selected by a
metadata-only rule — at least four libraries per cell type; exclusion of
disease, tumour, treatment, infection and whole-blood categories and of lineages
matching the development cells; within each lineage stratum, the cell type with
the most libraries. The set covers adipocyte, endothelial, epithelial
(two types), fibroblast (two types), immune (B and T cells), muscle (skeletal and
smooth) and bone-marrow mesenchymal stem cells.

One held-out cell type, CD19+ B cells, shares the B lineage with the GM12878
cell line. GM12878 is an immortalised lymphoblastoid line rather than a primary
B-cell sample, so this pair was retained; a sensitivity analysis excluding
CD19+ B cells is reported in Supplementary Table S5.

### 4.5 Locus holdout

Quality-passing loci were grouped by chromosome and sorted by GC content; every
fifth locus was assigned to the held-out set (12,360 of 61,752, 20.0%). Held-out
loci were excluded from all development-cell training, validation and test
partitions. The model-selection chromosomes were unchanged (validation chr18,
chr19, chr14; test chr3, chr5, chr11; remaining autosomes for training).

### 4.6 Model families

**Composition (21 parameters).** Four mononucleotide and sixteen dinucleotide
frequencies, `LogisticRegression(max_iter=3000, C=1.0,
class_weight="balanced")`, trained on development-cell records and applied
without cell identity.

**DNA-CNN (307,265 parameters).** One-hot 2,000 bp input; a 32-dimensional
cell-identity embedding (four states: K562, HepG2, GM12878 and UNK) is mapped
by FiLM (`Linear(32 → 512)`, gain and bias); Conv1d(4→64, k=11) → maxpool 4 →
Conv1d(64→128, k=7) → maxpool 4 → Conv1d(128→256, k=5) → maxpool 5; FiLM is
applied after the third convolution; concatenation of max- and mean-pooled
activations → `Linear(512→128)` → GELU → Dropout(0.2) → `Linear(128→1)`.

**Nucleotide Transformer v2 50M with LoRA.** Backbone
`nucleotide-transformer-v2-50m-multi-species` (revision `81b29e5`), loaded via
`AutoModelForMaskedLM.from_pretrained(..., trust_remote_code=True).esm`. LoRA
rank 8, alpha 16, dropout 0.05, target modules `query` and `value`. Head:
attention-mask-weighted mean of the last hidden state concatenated with a
32-dimensional cell embedding → `Linear(hidden + 32 → 256)` → GELU →
Dropout(0.2) → `Linear(256 → 1)`.

**Composition-XGBoost.** The same 20 composition features used by the logistic
model were modelled with XGBoost (`max_depth=4`, learning rate 0.05, 2,000
trees with early stopping on validation AUPRC, subsample 0.8,
`colsample_bytree=0.8`, `scale_pos_weight = n_neg / n_pos`).

**Motif-XGBoost.** JASPAR 2024 CORE vertebrates non-redundant PFMs (879 motifs)
[36]<!--ref:JASPAR2024--><!--anchor:none:-->
were scanned against each 2 kb window with MOODS at p = 1e-3; per-motif hit
counts were used as features in an XGBoost model (`max_depth=6`, learning rate
0.05, 2,000 trees, early stopping on validation AUPRC).

**Composition+motif-XGBoost.** The 20 composition features and 879 motif counts
were concatenated and modelled with the same XGBoost configuration as
motif-XGBoost.

### 4.7 Training

CNN: 15 epochs, batch 64, AdamW, learning rate 1e-3, weight decay 0.01,
BCEWithLogits with `pos_weight = n_neg / n_pos`, checkpoint selected by
validation AUPRC, cell index replaced by UNK with probability 0.3; five seeds
(13, 37, 71, 101, 137).

NT-LoRA: four epochs, batch 32 with gradient accumulation 2, AdamW with weight
decay 0.01, head learning rate 2e-4, LoRA learning rate 2e-5, BCEWithLogits with
`pos_weight`, cell index replaced by UNK with probability 0.3, checkpoint
selected by validation AUPRC; the same five seeds. All models were trained on
the 80% retained loci only.

To separate pretraining from architecture and capacity, the same NT
configuration was also initialised from the model config without pretrained
weights (`AutoModelForMaskedLM.from_config`), with the same tokenizer, LoRA
modules, data, optimizer and four-epoch schedule, for the same five seeds. This
control is referred to as NT-random-init.

A cell-free control used the same NT configuration and training schedule but set
the cell index to UNK for every training record; five seeds were trained.

No hyperparameter search was performed. The CNN and NT budgets were fixed before
the locus-holdout analysis, so the comparison reflects these prespecified
budgets rather than independently tuned optima. The random-init control was not
tuned separately; it isolates the effect of pretrained initialization under the
same optimization budget, not the best achievable performance of a randomly
initialized model.

### 4.8 Evaluation

Each model scored every locus once. The two cell-conditioned models used UNK as
the cell index; the composition model has no cell input. This produces a
cell-agnostic ranking, which was then evaluated against each held-out cell's
label set. Because the deep models were trained with a cell embedding and
evaluated with UNK, their scores reflect a learned average-cell representation
rather than a strictly cell-free model. For each cell type and stratum we report
AUPRC and lift (AUPRC / prevalence); AUROC and the fraction of the ranking
needed for 50% recall are reported in Supplementary Table S7. AUPRC was
prespecified as the primary metric because precision-recall behaviour is
sensitive to class imbalance [31]<!--ref:Williams2021--><!--anchor:none:-->.
Strata were defined by locus provenance (unseen versus retained) and by
development-cell activity (active in at least one development cell versus none).
Cell-type values were macro-averaged.

### 4.9 Statistics

Paired comparisons used a 1 Mb genomic-block bootstrap (blocks defined by
chromosome and `floor(midpoint / 1e6)`), with 1,000 resamples per cell type and
3,000 hierarchical resamples across cell types for macro intervals. Resamples
without a positive were skipped. 95% percentile intervals are reported.
Per-cell intervals are descriptive and were not adjusted for multiplicity
across the 11 cell types and four strata; the macro interval is the primary
inferential quantity. The
unseen-minus-retained difference was tested by resampling 1 Mb blocks
independently within each subset, pairing the resulting per-cell differences,
and applying the same hierarchical cell-level bootstrap for the macro interval.
The
protocol's earlier 10,000-replicate specification was executed as 2,000 in the
full-locus analysis of the previous stage; the locus-holdout analysis reported
here uses 1,000 per-cell and 3,000 macro resamples.

### 4.10 Software and reproducibility

Python 3.11; PyTorch 2.13.0+cu130; Transformers 4.57.1; PEFT 0.20.0;
scikit-learn 1.9.0; NumPy 2.4.6; pandas 2.3.3. Training used one NVIDIA A100
40 GB with TF32 enabled for the locus-holdout runs. Code and trained weights are
not publicly released. They are available from the corresponding author on
reasonable request for non-commercial research, subject to the licence of the
base Nucleotide Transformer model. Aggregate outputs needed to reproduce the
reported tables and figures are provided as Supplementary Tables S1–S8.

### 4.11 Protocol notes

The locus holdout was introduced after an audit that examined chromosome-level
splits; the audit's chromosome-level conclusion was not reproduced under the
locus holdout and is retained as a sensitivity analysis. The full-locus model
reported in the previous stage reached 10.333-fold over all held-out loci; it is
consistent with the present result but shares candidate loci between training
and evaluation, so it is not used as a generalisation test.

### 4.12 Sequence-neighborhood overlap

The holdout is defined at the level of annotated loci, not at the level of
disjoint sequence windows. Across the 61,752 quality-passing loci, 18,047
(29.2%) have a 2 kb window that overlaps at least one other candidate window;
the median maximal overlap among overlapping loci is 1,083 bp. In the frozen
holdout, 3,028 of 12,360 held-out loci (24.5%) had a window overlapping a
retained window, 2,410 overlapped by at least 500 bp and 1,617 by at least
1,000 bp. Among 1,142 unseen positives, 445 (39.0%; per-cell mean 39.7%, range
31.2–50.0%) overlapped a retained window, and 32.8% overlapped by at least
500 bp. The sensitivity analysis removed either every overlapping held-out
locus or only those with at least 500 bp of overlap, from both the positive and
negative pools, and recomputed AUPRC and the paired 1 Mb block-bootstrap
interval.

## 5. Limitations

- One pretrained backbone; the baseline-dependent conclusion is demonstrated for
  this backbone and these two composition readouts.
- The NT advantage over logistic composition is not reproduced against XGBoost
  on the same 20 composition features on unseen development-active loci. The
  primary comparison should therefore be read with both logistic and nonlinear
  composition baselines.
- The random-init control uses the same four-epoch schedule and was not
  separately tuned; it isolates the benefit of pretrained initialisation under
  this budget, not pretraining in general.
- One label resource and one assay (FANTOM5 CAGE usage); cross-resource
  generalisation is untested.
- Candidate universe restricted to FANTOM5 loci; this is not genome-wide de novo
  enhancer discovery.
- Fixed 2 kb windows and binary labels; the 1-of-3-active class is excluded.
- Development-inactive loci contribute about 54 positives per cell; their
  interval excludes positive effects above 0.0013 but is estimated from a
  limited number of positives.
- The locus holdout is not a sequence-disjoint holdout: 39.0% of unseen
  positives overlap a retained window. Removing overlapping loci attenuated the
  development-active gain from +0.0810 to +0.0667 (any overlap) and +0.0543
  (at least 500 bp), while both intervals remained above zero; the
  development-inactive subset retained only 7 cell types after removal.
- One held-out cell type, CD19+ B cells, shares the B lineage with GM12878;
  excluding it changed the macro ΔAUPRC by less than 0.002 (Supplementary
  Table S5).
- The CNN and NT models were trained with cell identity available for 70% of
  records and evaluated with the UNK index; the composition model never received
  cell identity. This training-signal asymmetry may affect the model comparison.
- No experimental validation of predicted loci; independent assays such as
STARR-seq [32]<!--ref:Hansen2022--><!--anchor:none:--> carry their own
systematic errors [33]<!--ref:Das2023--><!--anchor:none:-->.
- No attribution or motif analysis was performed for the locus-holdout models
[34]<!--ref:Koo2021--><!--anchor:none:--> [35]<!--ref:Majdandzic2022--><!--anchor:none:-->.
- The model has no cell-state input at evaluation. The cell index is set to UNK,
  so the score is cell-agnostic at evaluation; the study measures cross-cell
  transfer of a shared activity component, not cell-type-specific prediction.

## Data and code availability

Code and trained weights are not publicly released. They are available from the
corresponding author on reasonable request for non-commercial research, subject
to the licence of the base Nucleotide Transformer model. Aggregate outputs
needed to reproduce the reported tables and figures, including cohort counts,
per-cell metrics, per-seed metrics, the CD19-exclusion sensitivity analysis and
the sequence-neighborhood audit, are provided as Supplementary Tables S1–S8.
Upstream data are the FANTOM5 hg38 enhancer atlas, the hg38 primary assembly and
the ENCODE hg38 blacklist v2, with per-file URLs and checksums in the
supplementary data package.

## Figure legends

**Fig. 1 | A genomic locus anchors the locus-held-out, cell-held-out design.**
(a) A representative FANTOM5 enhancer (chr12:54,304,149–54,304,755; 606 bp)
within its 2 kb window. ATAC tracks are shown for context only and are not model
inputs; the locus is ATAC-open in all three development cells (filled blocks)
but CAGE-defined active only in K562 (orange peak and filled marker). The blocks
mark called intervals, not signal amplitudes. The locus
belongs to the 20% held-out set, is shown for label definition only and is
excluded from training; the bottom shows a 16-bp sequence excerpt.
(b) CAGE-defined activity of 24 FANTOM5 loci (rows) across three development
cells and 11 held-out primary cell types (columns), grouped by locus split and
by the cell set in which each locus is active. (c) The three model families
compared on identical inputs and splits; stronger composition and motif
baselines are compared in Fig. 2. (d) Evaluation crosses locus
provenance (retained versus unseen) with development-cell activity, giving four
strata.

**Fig. 2 | Stronger composition baselines remove the apparent NT advantage.**
Macro AUPRC for five models across four locus strata: logistic composition,
composition-XGBoost, motif-XGBoost, composition+motif-XGBoost and NT-LoRA.
Bars are macro means over held-out cell types; the number of cell types is
shown above each panel. (a) Unseen loci, all. (b) Unseen loci,
development-active. (c) Retained loci, development-active. (d) Unseen loci,
development-inactive (9 cell types with at least 10 positives). Δ is the macro
NT − composition-XGBoost difference with its 95% paired 1 Mb block-bootstrap
confidence interval. The dashed horizontal line marks the prevalence baseline.

**Fig. 3 | The NT is indistinguishable from composition-XGBoost on unseen
development-active loci.** (a) Per-cell paired AUPRC on unseen
development-active loci; x is composition-XGBoost and y is NT-LoRA. The dashed
diagonal marks equal performance; point size is the number of positives.
Filled points are above the diagonal (NT better), open points below. (b) Macro
NT − composition-XGBoost ΔAUPRC with 95% paired 1 Mb block-bootstrap confidence
intervals for the four strata. The diamond and vertical line show the macro
estimate and interval; the dashed vertical line marks zero. The macro
difference on unseen development-active loci is −0.0011 (95% CI −0.0358 to
+0.0365). Per-cell point estimates are shown for the comparison of primary
interest; per-cell intervals were not used for inference.

**Fig. 4 | Pretrained initialization is required for the logistic-baseline contrast on unseen loci under this budget.** (a) Macro
AUPRC on retained and unseen development-active loci for logistic composition (grey
dashed), random-init NT (teal) and pretrained NT (navy); thin lines are
individual cell types and thick lines are macro means. (b) Per-cell
difference in AUPRC between pretrained and random-init NT on retained (x) and
unseen (y) loci; the dashed diagonal marks equal gains and all 11 cell types
lie above it. n = 11 primary cell types; exact values and intervals are in
Table S4. The random-init model matches the pretrained model on retained loci
but falls below composition on unseen loci, consistent with locus-level
memorization and a pretraining-dependent contrast against the logistic
composition baseline.

**Supplementary Fig. S1 | Chromosome-level sensitivity.** (a) Macro ΔAUPRC by
locus provenance and chromosome provenance; labels give the number of cell
types with a positive Δ over the number evaluated. (b) Per-cell ΔAUPRC in the
same four groups; black bars are group means. Held-out-chromosome subsets
contain about 14 to 46 positives per cell and are underpowered; the locus-level
holdout is the primary generalisation test.

**Supplementary Fig. S2 | Overlap of active-locus sets across development and
held-out cells.** UpSet plot of the top 25 multi-set intersections of active
loci across the three development cell lines and 11 held-out primary cell
types. Navy intersections contain both development and held-out cells, grey
intersections are development-only and gold intersections are held-out-only;
singleton sets are omitted. Bar charts show intersection and set sizes. The
plot visualizes the shared and cell-specific components used in the main text.

**Supplementary Fig. S3 | Sequence-disjoint sensitivity.** (a) Macro
NT−composition ΔAUPRC on unseen development-active loci for all held-out loci
and after removing held-out loci whose 2 kb windows overlapped a retained
window (any overlap or at least 500 bp). (b) Same for all unseen loci. Points
are observed macro means and bars are 95% paired 1 Mb block-bootstrap
intervals; n = 11 cell types. The development-inactive subset retained 7 cell
types after overlap removal (Table S6).

**Table 1 | Cohort and label counts.** Development-cell counts after the 20%
locus holdout, and counts for the 11 held-out primary cell types (library
number, activity threshold, positives, negatives and prevalence).

**Table 2 | Four-stratum performance.** Mean AUPRC and lift for the three
models in each stratum, with per-cell means for the number of loci and positives,
and the macro ΔAUPRC (NT − logistic composition; point estimate of the per-cell mean)
with its 95% confidence interval.

**Table 3 | Stronger composition and motif baselines.** AUPRC for logistic
composition, composition-XGBoost, motif-XGBoost, composition+motif-XGBoost and
NT-LoRA in the four strata; Δ is the macro NT − composition-XGBoost difference
with its 95% confidence interval.

**Supplementary Table S8 | Cell-free NT control.** Macro NT − logistic
composition differences for an NT trained with UNK for every record.

## References

References are numbered in order of first appearance and formatted in Vancouver style.

1. Yang JH, Hansen AS. Enhancer selectivity in space and time: from enhancer-promoter interactions to promoter activation. *Nat Rev Mol Cell Biol*. 2024;25:574-591. doi:10.1038/s41580-024-00710-6
2. Blayney JW, Francis H, Rampasekova A, Camellato B, Mitchell L, Stolper R, et al. Super-enhancers include classical enhancers and facilitators to fully activate gene expression. *Cell*. 2023;186:5826-5839.e18. doi:10.1016/j.cell.2023.11.030
3. Andersson R, et al. An atlas of active enhancers across human cell types and tissues. *Nature*. 2014;507:455-461. doi:10.1038/nature12787
4. Kawaji H, Kasukawa T, Forrest A, Carninci P, Hayashizaki Y. The FANTOM5 collection, a data series underpinning mammalian transcriptome atlases in diverse cell types. *Sci Data*. 2017;4:170113. doi:10.1038/sdata.2017.113
5. Forrest ARR, et al. A promoter-level mammalian expression atlas. *Nature*. 2014;507:462-470. doi:10.1038/nature13182
6. Narita T, Ito S, Higashijima Y, Chu WK, Neumann K, Walter J, et al. Enhancers are activated by p300/CBP activity-dependent PIC assembly, RNAPII recruitment, and pause release. *Mol Cell*. 2021;81:2166-2182.e6. doi:10.1016/j.molcel.2021.03.008
7. Fong N, Sheridan RM, Ramachandran S, Bentley DL. The pausing zone and control of RNA polymerase II elongation by Spt5: implications for the pause-release model. *Mol Cell*. 2022;82:3632-3645.e4. doi:10.1016/j.molcel.2022.09.001
8. Mimoso CA, Goldman SR. PRO-seq: precise mapping of engaged RNA Pol II at single-nucleotide resolution. *Curr Protoc*. 2023;3:e961. doi:10.1002/cpz1.961
9. Han Z, Li W. Enhancer RNA: what we know and what we can achieve. *Cell Prolif*. 2022;55:e13202. doi:10.1111/cpr.13202
10. Grandi FC, Modi H, Kampman L, Corces MR. Chromatin accessibility profiling by ATAC-seq. *Nat Protoc*. 2022;17:1518-1552. doi:10.1038/s41596-022-00692-9
11. ENCODE Project Consortium. Expanded encyclopaedias of DNA elements in the human and mouse genomes. *Nature*. 2020;583:699-710. doi:10.1038/s41586-020-2493-4
12. Kang Y, Kim YW, Kang J, Kim A. Histone H3K4me1 and H3K27ac play roles in nucleosome eviction and eRNA transcription, respectively, at enhancers. *FASEB J*. 2021;35:e21781. doi:10.1096/fj.202100488r
13. Penzar D, Nogina D, Noskova E, Zinkevich A, Meshcheryakov G, Lando A, et al. LegNet: a best-in-class deep learning model for short DNA regulatory regions. *Bioinformatics*. 2023;39:btad457. doi:10.1093/bioinformatics/btad457
14. Liu Q, Hua K, Zhang X, Wong WH, Jiang R. DeepCAGE: incorporating transcription factors in genome-wide prediction of chromatin accessibility. *Genomics Proteomics Bioinformatics*. 2022;20:496-507. doi:10.1016/j.gpb.2021.08.015
15. Smith GD, Ching WH, Cornejo-Páramo P, Wong ES. Decoding enhancer complexity with machine learning and high-throughput discovery. *Genome Biol*. 2023;24:116. doi:10.1186/s13059-023-02955-4
16. Avsec Ž, et al. Effective gene expression prediction from sequence by integrating long-range interactions. *Nat Methods*. 2021;18:1196-1203. doi:10.1038/s41592-021-01252-x
17. Avsec Ž, Weilert M, Shrikumar A, Krueger S, Alexandari A, Dalal K, et al. Base-resolution models of transcription-factor binding reveal soft motif syntax. *Nat Genet*. 2021;53:354-366. doi:10.1038/s41588-021-00782-6
18. Chen KM, Wong AK, Troyanskaya OG, Zhou J. A sequence-based global map of regulatory activity for deciphering human genetics. *Nat Genet*. 2022;54:940-949. doi:10.1038/s41588-022-01102-2
19. Yuan H, Kelley DR. scBasset: Sequence-based modeling of single-cell ATAC-seq using convolutional neural networks. *Nat Methods*. 2022;19:1088-1096. doi:10.1038/s41592-022-01562-8
20. Yue T, Wang Y, Zhang L, Gu C, Xue H, Wang W, et al. Deep learning for genomics: from early neural nets to modern large language models. *Int J Mol Sci*. 2023;24:15858. doi:10.3390/ijms242115858
21. Routhier E, Mozziconacci J. Genomics enters the deep learning era. *PeerJ*. 2022;10:e13613. doi:10.7717/peerj.13613
22. Benegas G, Ye C, Albors C, Li JC, Song YS. Genomic language models: Opportunities and challenges. *Trends Genet*. 2025;41:286-302. doi:10.1016/j.tig.2024.11.013
23. Chen S, Gan M, Lv H, Jiang R. DeepCAPE: A deep convolutional neural network for the accurate prediction of enhancers. *Genomics Proteomics Bioinformatics*. 2021;19:565-577. doi:10.1016/j.gpb.2019.04.006
24. Thibodeau A, Khetan S, Eroglu A, Tewhey R, Stitzel ML, Ucar D. CoRE-ATAC: A deep learning model for the functional classification of regulatory elements from single cell and bulk ATAC-seq data. *PLoS Comput Biol*. 2021;17:e1009670. doi:10.1371/journal.pcbi.1009670
25. Zhang Z, Feng F, Qiu Y, Liu J. A generalizable framework to comprehensively predict epigenome, chromatin organization, and transcriptome. *Nucleic Acids Res*. 2023;51:5931-5947. doi:10.1093/nar/gkad436
26. Javed N, et al. A multi-modal transformer for cell type-agnostic regulatory predictions. *Cell Genom*. 2025;5:100762. doi:10.1016/j.xgen.2025.100762
27. Dalla-Torre H, et al. Nucleotide Transformer: Building and evaluating robust foundation models for human genomics. *Nat Methods*. 2025;22:287-297. doi:10.1038/s41592-024-02523-z
28. Ji Y, Zhou Z, Liu H, Davuluri RV. DNABERT-2: Efficient foundation model and benchmark for multi-species genome. In: *International Conference on Learning Representations*; 2024. Available from: https://arxiv.org/abs/2306.15006.
29. Hu EJ, et al. LoRA: Low-rank adaptation of large language models. In: *International Conference on Learning Representations*; 2022. Available from: https://openreview.net/forum?id=nZeVKeeFYf9.
30. Whalen S, Schreiber J, Noble WS, Pollard KS. Navigating the pitfalls of applying machine learning in genomics. *Nat Rev Genet*. 2022;23:169-181. doi:10.1038/s41576-021-00434-9
31. Williams CKI. The effect of class imbalance on precision-recall curves. *Neural Comput*. 2021;33:853-857. doi:10.1162/neco_a_01362
32. Hansen TJ, Hodges E. ATAC-STARR-seq reveals transcription factor-bound activators and silencers within chromatin-accessible regions of the human genome. *Genome Res*. 2022;32:1529-1541. doi:10.1101/gr.276766.122
33. Das M, Hossain A, Banerjee D, Praul CA, Girirajan S. Challenges and considerations for reproducibility of STARR-seq assays. *Genome Res*. 2023;33:479-495. doi:10.1101/gr.277204.122
34. Koo PK, Majdandzic A, Ploenzke M, Anand P, Paul SB. Global importance analysis: an interpretability method to quantify importance of genomic features in deep neural networks. *PLoS Comput Biol*. 2021;17:e1008925. doi:10.1371/journal.pcbi.1008925
35. Majdandzic A, Rajesh C, Tang A, Toneyan S, Labelson EL, Tripathy RK, et al. Selecting deep neural networks that yield consistent attribution-based interpretations for genomics. In: *Proceedings of Machine Learning Research*; 2022. p. 131-149. Available from: https://proceedings.mlr.press/v200/majdandzic22a.html.
36. Rauluseviciute I, Riudavets-Puig R, Blanc-Mathieu R, Castro-Mondragon JA, Ferenc K, Kumar V, et al. JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles. *Nucleic Acids Res*. 2024;52:D174-D182. doi:10.1093/nar/gkad1059

---

## 起草状态（定稿前删除）

- 主结果：五种子位点留出实验（`v2/locus_holdout_20261018/`）。
- 所有数字来自 `LH_posthoc_metrics.tsv`、`LH_posthoc_macro.tsv`、
  `LH_posthoc_deltas.tsv`、`LH_results.tsv`。
- 图表已生成：`paper/figures_v3/`（PNG/PDF/SVG/TIFF）与
  `paper/tables_v3/`（Table 1、Table 2）。
- 35 条引用均收录于 `paper/references.bib`，并于 2026-10-18 通过 Crossref 与
  PubMed/Europe PMC 逐条核验；其中 31 条发表于 2021 年及以后，保留的 4 条旧文献
  为 FANTOM5 原始数据（3 条）与 ENCODE 资源（1 条）。
