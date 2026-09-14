# Stage 9 region-selection rule

## Frozen rule

Before any Stage 9 pathway or cell-type result was generated, the regional
interpretation set was fixed as the **top 10 AAL3 parcels by the
Stage 7 `robustness_rank` column**. Ties would be resolved by ascending
`region_id`. No anatomical name, Parkinson relevance, Stage 9 pathway result,
or Stage 8 validation value enters selection.

The Stage 7 result itself is not re-estimated. A high rank does not imply that
a parcel passed joint FDR; the original `spatial_robustness_label` is retained
and reported. This is therefore a ranked exploratory interpretation of frozen
regions, not a new discovery threshold.

Within each selected parcel, the top 20 genes are defined
by descending additive contribution `expression × L2G / sum(L2G)`. The heatmap
uses the 15 genes with the highest mean contribution across these
parcels. These constants are set in `src/stage_09_biological_interpretation.py`.

## Selected regions

| robustness_rank | region_id | region_name | stage6_z | stage6_fdr | spatial_null_fdr | spatial_robustness_label |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 28 | OFCant R | 3.0566367753650683 | 0.1793820617938206 | 0.7623737626237376 | not_robust |
| 2 | 20 | Frontal Sup Medial R | 2.633686530481111 | 0.2690730926907309 | 0.7623737626237376 | not_robust |
| 3 | 141 | Thal MGN L | 2.383839881049032 | 0.3242675732426757 | 0.7623737626237376 | not_robust |
| 4 | 146 | Thal PuM R | 2.3598755927506923 | 0.3242675732426757 | 0.7623737626237376 | not_robust |
| 5 | 4 | Frontal Sup 2 R | 2.093181760563216 | 0.5077892210778922 | 0.8484618204846182 | not_robust |
| 6 | 37 | Cingulate Mid L | 1.7301977851209205 | 0.6684664866846648 | 0.8484618204846182 | not_robust |
| 7 | 45 | Amygdala L | 1.5387176061019034 | 0.7155570157269987 | 0.8484618204846182 | not_robust |
| 8 | 128 | Thal VL R | 1.8210236317945232 | 0.6684664866846648 | 0.8484618204846182 | not_robust |
| 9 | 41 | Hippocampus L | 1.741246693186499 | 0.6684664866846648 | 0.8484618204846182 | not_robust |
| 10 | 11 | Frontal Inf Orb 2 L | 1.3514219379593178 | 0.7155570157269987 | 0.8484618204846182 | not_robust |
