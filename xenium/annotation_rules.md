# Xenium annotation rules and input lineage

All cell labels are study-specific curated segment identities. “Endothelial-containing segment” does not imply every segmentation unit contains only one endothelial cell. `annotation_cluster_maps.tsv` lists the explicit cluster-to-label vectors. The selected code cells are preserved under `source_evidence/` with exact provenance; they depend on the corresponding preclustered objects and are not an automatic annotation pipeline for new data.

## Preprocessing shared across organs

The source loads Young and Aged Xenium outputs separately, retains cell segmentations and centroid coordinates, omits molecule coordinates, prefixes cell IDs by organ/age, and merges the two ages. It retains units with `nCount_Xenium > 5 & nFeature_Xenium > 3`, then runs `NormalizeData()`, `ScaleData(features = rownames(object))` and `RunPCA(features = rownames(object), npcs = 30)` on the Xenium assay. The notebooks compare dimensionality/clustering settings before saving a working object. The selected full-organ settings are Heart PCs 1–20/resolution 1.0, Muscle PCs 1–20/resolution 1.2, Adipose and Liver PCs 1–15/resolution 1.2. The source does not contain a documented quantitative optimization criterion for these choices; none is inferred here.

The normalization notebooks patch the body of Seurat's `LoadXenium()` to remove an `Unassigned Codeword` mapping. This is an original environment-specific workaround. The proximity runner starts from final annotated RDS files and does not patch any package. Raw-data reconstruction needs the original compatible loader and a validated import before these downstream steps.

## Heart

`Heart_normalize_naming.ipynb`, cell 17, applies the 30-cluster map and flags either clusters 6/21 or normalized Ptprc ≥ 4 as hematopoietic. This subset is normalized, scaled and reduced again; PCs 1–20 and resolution 0.4 define its selected subclusters (cells 20–22).

Cell 27 starts hematopoietic labels as `Monocyte/macrophage lineage` and then applies these ordered rules:

| Subcluster restriction | Marker rule | Label |
|---|---|---|
| 3 | Ms4a1 ≥ 5.0 | B cell lineage |
| 3 | Cd3d ≥ 4.5 and Ms4a1 < 5.0 | T cell lineage |
| 10 or 13 | Ppbp ≥ 6.0 | Platelet-rich |
| 10 or 13 | Ppbp < 6.0 and Csf3r ≥ 5.0 | Neutrophil |
| 5 or 12 | Itgax ≥ 4.5 | DC-associated |
| 5 or 12 | Cd209a ≥ 5.0 | moDC-like |
| 5 or 12 | Siglech ≥ 4.5 and Cd209a < 5.0 | pDC-associated |

Later matching rules overwrite earlier ones. Cell 31 returns the subset labels to the full object by cell ID. These values are empirical annotation cutoffs in normalized expression units; their original selection rationale was not recorded in the source cells.

`Heart_normalize_naming_Endocardium_refined.ipynb`, cells 5–9, uses manual CSV exports under `/work/Xenium/Heart_endocardium`. Each CSV must have a `Cell ID` column; filenames beginning `young_`/`aged_` determine the added `Young_Heart_`/`Aged_Heart_` prefix. Only cells originally labeled `Venous-like large-vessel EC-rich` are eligible. Selected arterial cells become `Artery-like large-vessel EC-rich`; selected venous cells retain `Venous-like large-vessel EC-rich`; remaining eligible cells become `Endocardial cell-rich`. The final object is saved in cell 20 as `heart_qc_annotated_refined_updated.rds`. The selection CSVs are necessary to reproduce this manual refinement and are not supplied in the code folder.

## Muscle

`Muscle_normalize_naming.ipynb`, cell 14, applies its 28-cluster map. The hematopoietic subset is renormalized and reclustered using PCs 1–15; resolution 0.8 is used in the final subcluster label vector (cell 22). Cells 23–25 then set Cd3d ≥ 5.0 to `T cell lineage`, Ms4a1 ≥ 5.5 to `B cell lineage` (the latter assignment takes precedence), and rename `Hematopoietic-associated mixed` to `B cell-associated mixed`. Cell 32 returns final labels to the full object; cell 37 saves `muscle_qc_annotated_refined.rds`.

The proximity analysis includes only the exact `Monocyte/macrophage lineage` label, not every macrophage-associated mixed segment or Trem2-enriched subtype. Expanding that population would change the analysis and is not done in the consolidated runner.

## Adipose

`Adipose_normalize_naming.ipynb`, cell 8, applies its 24-cluster map. The hematopoietic subset is renormalized and its selected graph uses PCs 1–20, resolution 1.0 (cells 17, 20–21). Cell 24 applies the 16-subcluster map, and cell 26 returns the labels to the full object. Cell 32 saves the refined object. `Adipose_IFN_Pira2.ipynb` subsequently adds expression/module-score information and saves `adipose_qc_annotated_refined_updated.rds`; that is the input used by the proximity notebooks.

Only `Monocyte/macrophage lineage` is included as the query population; mixed labels are excluded as in the original final endothelial proximity analysis.

## Liver recorded outside the default manuscript runner

`Liver_normalize_naming.ipynb` uses a 25-cluster map, then a hematopoietic subset with PCs 1–15 and resolution 0.8. Its later rules include final Cd3d ≥ 3.0 and Ms4a1 ≥ 4.0 thresholds (cell 41), followed by a relabeling of hematopoietic subclusters 0/16 as `Fibroblast-rich` (cells 55–56). The full-organ label set uses `LSEC-rich`, `Kupffer cell` and `Macrophage lineage`. However, `Liver_IFN_Pira2_proximity.ipynb`, cell 5, queries `Monocyte/macrophage lineage`. The intended query population and final figure scope need confirmation; this inconsistency is not silently resolved.

## Input identity and verification

The eight Xenium output directory labels are listed in `sample_manifest.tsv`. They are sample identifiers from filenames, not a verified cross-organ animal mapping. The one-mouse-per-organ-age statement follows the response letter. The manuscript's animal correspondence table remains the authority for matching organs from the same mouse.

Saved notebook outputs establish which cutoffs were calculated in those sessions but do not prove that the current final RDS files have unchanged cells, expression layers or annotations. The consolidated runner exports annotation counts and the recalculated cutoff summary so that these can be compared after the original inputs are located.
