# Input objects and analysis boundaries

| Input | Origin and purpose | Required contents |
|---|---|---|
| Capillary/sinusoidal age-comparison Seurat RDS | Final capillary/sinusoidal subset, normal diet, young and aged; organ-specific aging SCENIC | SCT corrected counts; unique cell IDs; Tissue and Age metadata |
| Capillary/sinusoidal diet-comparison Seurat RDS | Final capillary/sinusoidal subset, aged normal diet and HFD; diet SCENIC | SCT corrected counts; Tissue and Diet metadata |
| Artery age/diet comparison RDS | Final artery subset for corresponding contrasts | SCT corrected counts; Tissue and Age or Diet metadata |
| Four-organ young capillary/sinusoidal RDS | Organ-specific regulon comparison in normal young cells | SCT corrected counts; Tissue metadata |
| All-sample capillary/sinusoidal RDS | Combined analysis of all organ/condition groups | SCT corrected counts; Type metadata |
| Organ endothelial RDS | Four age/diet groups in an organ; subset to capillary or sinusoidal before SCENIC | SCT corrected counts; Vessel_Type and Condition metadata |
| Endothelial RDS for IFN scores | All vascular endothelial subtypes; select Diet = Normal and score separately within each organ | RNA normalized data; Tissue, Diet, Condition, Vessel_Type |
| Aging-signature RDS | Aging-comparison capillary/sinusoidal object | SCT normalized data; Type metadata; all signature/control genes |
| Reactome-annotated endothelial RDS | Output of per-organ scGSVA/UCell pathway scoring and SCENIC AUC assembly | `regulonAUC_extended` assay; pathway metadata columns named in `Reactome_annotations.csv` |
| Organ `Merged_CellChat_prep.rds` | Organ-wide cell-type-prepared object, including nonendothelial populations | SCT corrected counts; Celltype_CellChat and Condition metadata |

The SCENIC manifest maps each run to its source-relative input path. The CellChat manifest records the input basename and condition; the runner accepts the full input path as an argument.

Cells and score/label arrays must be joined by cell identifier, not by an assumed shared row order. The generalized functions check or enforce those joins. Original biological replicates and sample identifiers should remain in metadata even where a descriptive cell-level figure uses Condition or Tissue for grouping.
