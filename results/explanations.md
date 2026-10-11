# Explanations for top-ranked bispecific target pairs

Provenance: deterministic

```
CD3E x TIGIT - checkpoint blockade with effector co-stimulation
---------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.914 (percentile among checkpoint|effector pairs, scored as immune combination;
  evidence coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 0.7% of malignant cells co-express both antigens (strongest in
  diffuse large B-cell lymphoma; measured across 4 dataset(s)). Co-incidence lift 52.36x versus
  independent expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 0.2% carry both at nTPM>=10. Individually the arms peak in myeloma
  (CD3E) and lymphoma (TIGIT).
SAFETY: the AND-gate peak in normal tissue is 4.57 (log1p nCPM); 2 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 1.48 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): t-cells (CD3E 425, TIGIT 96 nCPM); nk-cells
  (CD3E 137, TIGIT 59 nCPM); plasma cells (CD3E 8, TIGIT 10 nCPM); thymocytes (CD3E 84, TIGIT 3
  nCPM).
  Within the tumour itself, double-positive rate is -6.0% higher in malignant than in non-malignant
  cells - i.e. NOT selective for malignant cells in these samples.
MODEL: top additive contributions - cls_membrane__mean +239.319, cls_membrane__min -162.127,
  normal_n_expressed__mean +88.195, tumour_best_frac__min -58.629, taa_tumour_max_log +58.027,
  taa_window_vital_log +37.972.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
ERBB2 x MET - dual tumour-antigen targeting (avidity gain and escape prevention)
--------------------------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.902 (percentile among TAA|TAA pairs, scored as pathway redundancy; evidence
  coverage 84.0% of rubric weight).
EFFICACY: in single tumour cells, 36.2% of malignant cells co-express both antigens (strongest in
  colorectal cancer; measured across 5 dataset(s)). Co-incidence lift 1.73x versus independent
  expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 65.3% carry both at nTPM>=10. Individually the arms peak in
  bile_duct_cancer (ERBB2) and bile_duct_cancer (MET).
SAFETY: the AND-gate peak in normal tissue is 4.79 (log1p nCPM); 79 of 154 normal cell types express
  both above threshold, of which 14 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 2.14 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): renal connecting tubule cells (ERBB2 128, MET
  119 nCPM); renal collecting duct intercalated cells (ERBB2 109, MET 161 nCPM); papillary tip
  epithelial cells (ERBB2 108, MET 152 nCPM); respiratory secretory cells (ERBB2 99, MET 210 nCPM).
  Within the tumour itself, double-positive rate is 35.8% higher in malignant than in non-malignant
  cells - a favourable in-situ margin.
MODEL: top additive contributions - cls_membrane__min -78.346, tumour_mean_frac__min +46.701,
  tumour_n_types_hit__min -42.733, line_double_pos_frac +35.854, normal_max_log__mean +20.769,
  window_vital_log__mean -19.348.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
EGFR x MET - dual tumour-antigen targeting (avidity gain and escape prevention)
-------------------------------------------------------------------------------
STATUS: already documented - amivantamab (approved, dual signalling blockade). Recovered by the
  ranking rather than proposed by it.
Within-class score 0.890 (percentile among TAA|TAA pairs, scored as pathway redundancy; evidence
  coverage 84.0% of rubric weight).
EFFICACY: in single tumour cells, 23.5% of malignant cells co-express both antigens (strongest in
  colorectal cancer; measured across 5 dataset(s)). Co-incidence lift 3.25x versus independent
  expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 59.3% carry both at nTPM>=10. Individually the arms peak in
  bile_duct_cancer (EGFR) and bile_duct_cancer (MET).
SAFETY: the AND-gate peak in normal tissue is 6.26 (log1p nCPM); 84 of 154 normal cell types express
  both above threshold, of which 14 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 1.04 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): salivary basal cells (EGFR 657, MET 521
  nCPM); ocular epithelial cells (EGFR 739, MET 366 nCPM); epididymal basal cells (EGFR 266, MET 347
  nCPM); endometrial glandular cells (EGFR 241, MET 1018 nCPM).
  Within the tumour itself, double-positive rate is 23.2% higher in malignant than in non-malignant
  cells - a favourable in-situ margin.
MODEL: top additive contributions - cls_membrane__min -93.331, cls_membrane__mean +44.050,
  tumour_profile_corr +38.747, line_double_pos_frac +34.939, normal_max_log__mean +24.887,
  cls_transporter__mean +22.162.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
CD47 x TNFRSF17 - tumour-targeted checkpoint blockade (localises immune release to tumour)
------------------------------------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.882 (percentile among TAA|checkpoint pairs, scored as tumour localised
  checkpoint; evidence coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 8.0% of malignant cells co-express both antigens (strongest in
  diffuse large B-cell lymphoma; measured across 4 dataset(s)). Co-incidence lift 1.49x versus
  independent expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 7.4% carry both at nTPM>=10. Individually the arms peak in
  bile_duct_cancer (CD47) and myeloma (TNFRSF17).
SAFETY: the AND-gate peak in normal tissue is 4.90 (log1p nCPM); 7 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 1.53 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): plasma cells (CD47 134, TNFRSF17 566 nCPM);
  foveolar cells (CD47 71, TNFRSF17 28 nCPM); gastric chief cells (CD47 166, TNFRSF17 20 nCPM); pdcs
  (CD47 160, TNFRSF17 15 nCPM).
  Within the tumour itself, double-positive rate is 3.5% higher in malignant than in non-malignant
  cells - a favourable in-situ margin.
MODEL: top additive contributions - cls_membrane__min -62.298, cls_membrane__mean +39.091,
  window_vital_log__mean +25.355, normal_max_log__mean +16.899, taa_normal_vital_n_expressed
  +13.122, normal_tau__min +12.575.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
NOTCH2NLA x SIRPA - tumour-targeted checkpoint blockade (localises immune release to tumour)
--------------------------------------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.879 (percentile among TAA|checkpoint pairs, scored as tumour localised
  checkpoint; evidence coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 22.3% of malignant cells co-express both antigens (strongest in
  glioblastoma; measured across 4 dataset(s)). Co-incidence lift 1.56x versus independent expression
  - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 43.7% carry both at nTPM>=10. Individually the arms peak in
  bone_cancer (NOTCH2NLA) and kidney_cancer (SIRPA).
SAFETY: the AND-gate peak in normal tissue is 2.53 (log1p nCPM); 1 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 4.16 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): tuft cells (NOTCH2NLA 23, SIRPA 12 nCPM);
  medullary thymic epithelial cells (NOTCH2NLA 2, SIRPA 10 nCPM); pdcs (NOTCH2NLA 1, SIRPA 2 nCPM);
  conjunctival goblet cells (NOTCH2NLA 1, SIRPA 15 nCPM).
  Within the tumour itself, double-positive rate is 1.6% higher in malignant than in non-malignant
  cells - a favourable in-situ margin.
MODEL: top additive contributions - cls_membrane__min +48.775, normal_max_log__mean -39.029,
  cls_membrane__absdiff +23.369, window_vital_log__mean +9.404, taa_tumour_max_log -8.855,
  cls_membrane__mean -8.409.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
NOTCH2NLA x VSIR - tumour-targeted checkpoint blockade (localises immune release to tumour)
-------------------------------------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.879 (percentile among TAA|checkpoint pairs, scored as tumour localised
  checkpoint; evidence coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 4.2% of malignant cells co-express both antigens (strongest in
  glioblastoma; measured across 4 dataset(s)). Co-incidence lift 5.14x versus independent expression
  - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 20.6% carry both at nTPM>=10. Individually the arms peak in
  bone_cancer (NOTCH2NLA) and head_and_neck_cancer (VSIR).
SAFETY: the AND-gate peak in normal tissue is 3.18 (log1p nCPM); 1 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 3.95 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): tuft cells (NOTCH2NLA 23, VSIR 69 nCPM); pdcs
  (NOTCH2NLA 1, VSIR 50 nCPM); conjunctival goblet cells (NOTCH2NLA 1, VSIR 26 nCPM); podocytes
  (NOTCH2NLA 0, VSIR 38 nCPM).
  Within the tumour itself, double-positive rate is -16.3% higher in malignant than in non-malignant
  cells - i.e. NOT selective for malignant cells in these samples.
MODEL: top additive contributions - cls_membrane__min +29.975, window_vital_log__mean -13.998,
  sc_pboth_mean -13.268, cls_membrane__absdiff +13.075, window_vital_log__absdiff +10.884,
  sc_pboth__bladder_SCBO8 -8.803.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
HAVCR2 x NOTCH2NLA - tumour-targeted checkpoint blockade (localises immune release to tumour)
---------------------------------------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.878 (percentile among TAA|checkpoint pairs, scored as tumour localised
  checkpoint; evidence coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 1.2% of malignant cells co-express both antigens (strongest in
  glioblastoma; measured across 4 dataset(s)). Co-incidence lift 1.47x versus independent expression
  - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 1.6% carry both at nTPM>=10. Individually the arms peak in
  lymphoma (HAVCR2) and bone_cancer (NOTCH2NLA).
SAFETY: the AND-gate peak in normal tissue is 1.22 (log1p nCPM); 0 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 4.18 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): medullary thymic epithelial cells (HAVCR2 2,
  NOTCH2NLA 2 nCPM); pdcs (HAVCR2 6, NOTCH2NLA 1 nCPM); conjunctival goblet cells (HAVCR2 1,
  NOTCH2NLA 1 nCPM); podocytes (HAVCR2 19, NOTCH2NLA 0 nCPM).
  Within the tumour itself, double-positive rate is -17.8% higher in malignant than in non-malignant
  cells - i.e. NOT selective for malignant cells in these samples.
MODEL: top additive contributions - cls_membrane__min +63.660, cls_membrane__absdiff +27.759,
  normal_max_log__max -19.367, window_vital_log__mean -17.523, normal_max_log__mean -16.732,
  taa_tumour_max_log -14.854.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
CD47 x CD79B - tumour-targeted checkpoint blockade (localises immune release to tumour)
---------------------------------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.877 (percentile among TAA|checkpoint pairs, scored as tumour localised
  checkpoint; evidence coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 39.7% of malignant cells co-express both antigens (strongest in
  diffuse large B-cell lymphoma; measured across 4 dataset(s)). Co-incidence lift 2.17x versus
  independent expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 11.8% carry both at nTPM>=10. Individually the arms peak in
  bile_duct_cancer (CD47) and myeloma (CD79B).
SAFETY: the AND-gate peak in normal tissue is 2.98 (log1p nCPM); 1 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 3.45 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): b-cells (CD47 222, CD79B 19 nCPM); plasma
  cells (CD47 134, CD79B 2 nCPM); hematopoietic stem cells (CD47 119, CD79B 2 nCPM); monocytes (CD47
  125, CD79B 1 nCPM).
  Within the tumour itself, double-positive rate is 18.4% higher in malignant than in non-malignant
  cells - a favourable in-situ margin.
MODEL: top additive contributions - cls_membrane__min -84.325, normal_max_log__mean -63.828,
  sc_pboth__dlbcl +35.745, window_vital_log__mean +35.231, cls_membrane__mean +31.694,
  window_log__max +23.990.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
CD2 x TIGIT - checkpoint blockade with effector co-stimulation
--------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.868 (percentile among checkpoint|effector pairs, scored as immune combination;
  evidence coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 0.7% of malignant cells co-express both antigens (strongest in
  diffuse large B-cell lymphoma; measured across 4 dataset(s)). Co-incidence lift 32.22x versus
  independent expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 0.2% carry both at nTPM>=10. Individually the arms peak in
  leukemia (CD2) and lymphoma (TIGIT).
SAFETY: the AND-gate peak in normal tissue is 4.57 (log1p nCPM); 2 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 1.69 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): t-cells (CD2 525, TIGIT 96 nCPM); nk-cells
  (CD2 142, TIGIT 59 nCPM); plasma cells (CD2 7, TIGIT 10 nCPM); thymocytes (CD2 42, TIGIT 3 nCPM).
  Within the tumour itself, double-positive rate is -5.7% higher in malignant than in non-malignant
  cells - i.e. NOT selective for malignant cells in these samples.
MODEL: top additive contributions - cls_membrane__mean +203.482, cls_membrane__min -145.652,
  normal_n_expressed__mean +86.730, taa_tumour_max_log +59.610, tumour_best_frac__min -57.937,
  taa_window_vital_log +35.983.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
EGFR x TNFRSF1A - dual tumour-antigen targeting (avidity gain and escape prevention)
------------------------------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.867 (percentile among TAA|TAA pairs, scored as pathway redundancy; evidence
  coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 39.3% of malignant cells co-express both antigens (strongest in
  glioblastoma; measured across 5 dataset(s)). Co-incidence lift 3.67x versus independent expression
  - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 62.0% carry both at nTPM>=10. Individually the arms peak in
  bile_duct_cancer (EGFR) and bladder_cancer (TNFRSF1A).
SAFETY: the AND-gate peak in normal tissue is 5.19 (log1p nCPM); 109 of 154 normal cell types
  express both above threshold, of which 15 are dose-limiting lineages. Pairing lowers peak normal
  exposure by 2.11 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): syncytiotrophoblasts (EGFR 1482, TNFRSF1A 179
  nCPM); decidual stromal cells (EGFR 178, TNFRSF1A 197 nCPM); hepatocytes (EGFR 298, TNFRSF1A 170
  nCPM); esophageal basal cells (EGFR 276, TNFRSF1A 151 nCPM).
  Within the tumour itself, double-positive rate is 37.2% higher in malignant than in non-malignant
  cells - a favourable in-situ margin.
MODEL: top additive contributions - cls_membrane__min -79.329, cls_membrane__mean +42.556,
  normal_max_log__mean +20.041, window_vital_log__mean +16.618, taa_normal_frac_expressed -14.732,
  tumour_profile_corr +13.586.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
SDHA x SDHB - dual tumour-antigen targeting (avidity gain and escape prevention)
--------------------------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.863 (percentile among TAA|TAA pairs, scored as pathway redundancy; evidence
  coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 35.1% of malignant cells co-express both antigens (strongest in
  colorectal cancer; measured across 5 dataset(s)). Co-incidence lift 1.59x versus independent
  expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 100.0% carry both at nTPM>=10. Individually the arms peak in
  bile_duct_cancer (SDHA) and bile_duct_cancer (SDHB).
SAFETY: the AND-gate peak in normal tissue is 4.28 (log1p nCPM); 20 of 154 normal cell types express
  both above threshold, of which 10 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 1.58 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): distal convoluted tubule cells (SDHA 76, SDHB
  71 nCPM); renal collecting duct intercalated cells (SDHA 91, SDHB 63 nCPM); choroid plexus
  epithelial cells (SDHA 62, SDHB 67 nCPM); renal connecting tubule cells (SDHA 76, SDHB 60 nCPM).
  Within the tumour itself, double-positive rate is 21.0% higher in malignant than in non-malignant
  cells - a favourable in-situ margin.
MODEL: top additive contributions - cls_membrane__min +14.548, window_vital_log__mean -11.815,
  tumour_best_frac__min +9.329, line_and_max_log +8.674, cls_membrane__absdiff -7.816,
  window_log__mean +7.259.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
MET x TNFRSF1A - dual tumour-antigen targeting (avidity gain and escape prevention)
-----------------------------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.862 (percentile among TAA|TAA pairs, scored as pathway redundancy; evidence
  coverage 84.0% of rubric weight).
EFFICACY: in single tumour cells, 37.0% of malignant cells co-express both antigens (strongest in
  colorectal cancer; measured across 5 dataset(s)). Co-incidence lift 1.75x versus independent
  expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 69.4% carry both at nTPM>=10. Individually the arms peak in
  bile_duct_cancer (MET) and bladder_cancer (TNFRSF1A).
SAFETY: the AND-gate peak in normal tissue is 4.80 (log1p nCPM); 86 of 154 normal cell types express
  both above threshold, of which 15 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 2.12 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): epididymal basal cells (MET 347, TNFRSF1A 121
  nCPM); fallopian secretory cells (MET 116, TNFRSF1A 142 nCPM); vascular endothelial cells (MET
  114, TNFRSF1A 152 nCPM); foveolar cells (MET 130, TNFRSF1A 113 nCPM).
  Within the tumour itself, double-positive rate is 36.2% higher in malignant than in non-malignant
  cells - a favourable in-situ margin.
MODEL: top additive contributions - cls_membrane__min -66.503, tumour_n_types_hit__min -42.834,
  tumour_mean_frac__min +38.299, line_double_pos_frac +35.756, normal_max_log__mean +22.180,
  tumour_best_frac__min +15.882.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
FCGR3A x PDCD1 - checkpoint blockade with effector co-stimulation
-----------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.862 (percentile among checkpoint|effector pairs, scored as immune combination;
  evidence coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 0.2% of malignant cells co-express both antigens (strongest in
  glioblastoma; measured across 4 dataset(s)). Co-incidence lift 24.80x versus independent
  expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 0.0% carry both at nTPM>=10. Individually the arms peak in
  bile_duct_cancer (FCGR3A) and lymphoma (PDCD1).
SAFETY: the AND-gate peak in normal tissue is 1.10 (log1p nCPM); 0 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 6.36 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): t-cells (FCGR3A 22, PDCD1 2 nCPM); b-cells
  (FCGR3A 2, PDCD1 0 nCPM); loop of henle epithelial cells (FCGR3A 0, PDCD1 0 nCPM); microglia
  (FCGR3A 6, PDCD1 0 nCPM).
  Within the tumour itself, double-positive rate is -6.8% higher in malignant than in non-malignant
  cells - i.e. NOT selective for malignant cells in these samples.
MODEL: top additive contributions - cls_membrane__mean +235.481, cls_membrane__min -179.333,
  normal_n_expressed__mean +87.786, normal_max_log__mean -87.510, window_vital_log__mean +54.739,
  taa_tumour_max_log +52.679.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
CTLA4 x ICOS - checkpoint blockade with effector co-stimulation
---------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.861 (percentile among checkpoint|effector pairs, scored as immune combination;
  evidence coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 0.3% of malignant cells co-express both antigens (strongest in
  diffuse large B-cell lymphoma; measured across 4 dataset(s)). Co-incidence lift 134.04x versus
  independent expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 0.2% carry both at nTPM>=10. Individually the arms peak in
  skin_cancer (CTLA4) and lymphoma (ICOS).
SAFETY: the AND-gate peak in normal tissue is 4.48 (log1p nCPM); 1 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 1.22 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): t-cells (CTLA4 87, ICOS 300 nCPM); cdc (CTLA4
  7, ICOS 5 nCPM); neutrophil progenitors (CTLA4 3, ICOS 4 nCPM); nk-cells (CTLA4 3, ICOS 23 nCPM).
  Within the tumour itself, double-positive rate is -3.1% higher in malignant than in non-malignant
  cells - i.e. NOT selective for malignant cells in these samples.
MODEL: top additive contributions - normal_tau__absdiff -206.613, cls_membrane__min -179.991,
  cls_membrane__mean +87.159, tumour_best_frac__min -53.716, window_vital_log__mean -36.920,
  tumour_mean_log__min -29.638.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
PPP1R18 x TNFRSF4 - tumour-antigen x immune-effector engager (redirects effector cells to tumour)
-------------------------------------------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.851 (percentile among TAA|effector pairs, scored as effector engager; evidence
  coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 1.1% of malignant cells co-express both antigens (strongest in
  diffuse large B-cell lymphoma; measured across 4 dataset(s)). Co-incidence lift 1.87x versus
  independent expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 1.2% carry both at nTPM>=10. Individually the arms peak in
  head_and_neck_cancer (PPP1R18) and lymphoma (TNFRSF4).
SAFETY: the AND-gate peak in normal tissue is 2.10 (log1p nCPM); 0 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 2.75 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): nk-cells (PPP1R18 7, TNFRSF4 20 nCPM);
  t-cells (PPP1R18 3, TNFRSF4 34 nCPM); monocytes (PPP1R18 4, TNFRSF4 3 nCPM); macrophages (PPP1R18
  2, TNFRSF4 15 nCPM).
  Within the tumour itself, double-positive rate is -2.3% higher in malignant than in non-malignant
  cells - i.e. NOT selective for malignant cells in these samples.
MODEL: top additive contributions - cls_membrane__min +63.509, cls_membrane__absdiff +28.701,
  window_vital_log__mean -16.920, taa_tumour_max_log -14.666, normal_max_log__mean -10.979,
  normal_max_log__max -10.943.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
CD28 x PPP1R18 - tumour-antigen x immune-effector engager (redirects effector cells to tumour)
----------------------------------------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.851 (percentile among TAA|effector pairs, scored as effector engager; evidence
  coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 0.9% of malignant cells co-express both antigens (strongest in
  diffuse large B-cell lymphoma; measured across 4 dataset(s)). Co-incidence lift 6.21x versus
  independent expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 3.8% carry both at nTPM>=10. Individually the arms peak in myeloma
  (CD28) and head_and_neck_cancer (PPP1R18).
SAFETY: the AND-gate peak in normal tissue is 1.63 (log1p nCPM); 0 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 3.37 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): nk-cells (CD28 4, PPP1R18 7 nCPM); t-cells
  (CD28 147, PPP1R18 3 nCPM); platelets (CD28 8, PPP1R18 3 nCPM); macrophages (CD28 16, PPP1R18 2
  nCPM).
  Within the tumour itself, double-positive rate is -4.8% higher in malignant than in non-malignant
  cells - i.e. NOT selective for malignant cells in these samples.
MODEL: top additive contributions - cls_membrane__min +52.259, normal_max_log__mean -26.742,
  cls_membrane__absdiff +24.206, normal_max_log__max -21.833, window_vital_log__mean +14.247,
  taa_window_vital_log +11.034.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
KLRK1 x PPP1R18 - tumour-antigen x immune-effector engager (redirects effector cells to tumour)
-----------------------------------------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.851 (percentile among TAA|effector pairs, scored as effector engager; evidence
  coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 0.9% of malignant cells co-express both antigens (strongest in
  diffuse large B-cell lymphoma; measured across 4 dataset(s)). Co-incidence lift 1.09x versus
  independent expression - little evidence of same-cell preference.
  Across 1,206 cancer cell lines, 3.1% carry both at nTPM>=10. Individually the arms peak in
  lymphoma (KLRK1) and head_and_neck_cancer (PPP1R18).
SAFETY: the AND-gate peak in normal tissue is 2.93 (log1p nCPM); 1 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 0.86 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): oligodendrocyte progenitor cells (KLRK1 25,
  PPP1R18 18 nCPM); nk-cells (KLRK1 22, PPP1R18 7 nCPM); t-cells (KLRK1 12, PPP1R18 3 nCPM);
  oligodendrocytes (KLRK1 8, PPP1R18 2 nCPM).
  Within the tumour itself, double-positive rate is -3.9% higher in malignant than in non-malignant
  cells - i.e. NOT selective for malignant cells in these samples.
MODEL: top additive contributions - cls_membrane__min +50.492, cls_membrane__absdiff +23.906,
  normal_max_log__mean -23.815, normal_max_log__max -12.805, window_vital_log__mean +10.776,
  taa_tumour_max_log -10.389.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
CD3E x PPP1R18 - tumour-antigen x immune-effector engager (redirects effector cells to tumour)
----------------------------------------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.851 (percentile among TAA|effector pairs, scored as effector engager; evidence
  coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 2.0% of malignant cells co-express both antigens (strongest in
  diffuse large B-cell lymphoma; measured across 4 dataset(s)). Co-incidence lift 2.33x versus
  independent expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 3.2% carry both at nTPM>=10. Individually the arms peak in myeloma
  (CD3E) and head_and_neck_cancer (PPP1R18).
SAFETY: the AND-gate peak in normal tissue is 2.10 (log1p nCPM); 0 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 3.95 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): nk-cells (CD3E 137, PPP1R18 7 nCPM); t-cells
  (CD3E 425, PPP1R18 3 nCPM); platelets (CD3E 3, PPP1R18 3 nCPM); neutrophils (CD3E 2, PPP1R18 4
  nCPM).
  Within the tumour itself, double-positive rate is -8.4% higher in malignant than in non-malignant
  cells - i.e. NOT selective for malignant cells in these samples.
MODEL: top additive contributions - normal_max_log__mean -47.083, cls_membrane__min +44.715,
  cls_membrane__absdiff +21.799, taa_window_vital_log +15.906, window_vital_log__mean +15.091,
  window_vital_log__min +10.704.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
CD40 x PPP1R18 - tumour-antigen x immune-effector engager (redirects effector cells to tumour)
----------------------------------------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.851 (percentile among TAA|effector pairs, scored as effector engager; evidence
  coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 24.5% of malignant cells co-express both antigens (strongest in
  diffuse large B-cell lymphoma; measured across 5 dataset(s)). Co-incidence lift 2.10x versus
  independent expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 27.2% carry both at nTPM>=10. Individually the arms peak in
  kidney_cancer (CD40) and head_and_neck_cancer (PPP1R18).
SAFETY: the AND-gate peak in normal tissue is 3.24 (log1p nCPM); 1 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 1.79 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): microglia (CD40 24, PPP1R18 43 nCPM);
  bergmann glia (CD40 8, PPP1R18 4 nCPM); late spermatids (CD40 17, PPP1R18 4 nCPM); papillary tip
  epithelial cells (CD40 16, PPP1R18 4 nCPM).
  Within the tumour itself, double-positive rate is 6.4% higher in malignant than in non-malignant
  cells - a favourable in-situ margin.
MODEL: top additive contributions - cls_membrane__min +47.830, sc_pboth__dlbcl +33.638,
  normal_max_log__mean -31.047, cls_membrane__absdiff +24.703, normal_max_log__max -13.988,
  taa_window_vital_log +11.481.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
CD3E x CTLA4 - checkpoint blockade with effector co-stimulation
---------------------------------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.850 (percentile among checkpoint|effector pairs, scored as immune combination;
  evidence coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 0.4% of malignant cells co-express both antigens (strongest in
  diffuse large B-cell lymphoma; measured across 4 dataset(s)). Co-incidence lift 83.77x versus
  independent expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 0.2% carry both at nTPM>=10. Individually the arms peak in myeloma
  (CD3E) and skin_cancer (CTLA4).
SAFETY: the AND-gate peak in normal tissue is 4.48 (log1p nCPM); 1 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 1.57 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): t-cells (CD3E 425, CTLA4 87 nCPM); neutrophil
  progenitors (CD3E 5, CTLA4 3 nCPM); cdc (CD3E 3, CTLA4 7 nCPM); nk-cells (CD3E 137, CTLA4 3 nCPM).
  Within the tumour itself, double-positive rate is -2.9% higher in malignant than in non-malignant
  cells - i.e. NOT selective for malignant cells in these samples.
MODEL: top additive contributions - cls_membrane__min -127.606, cls_membrane__mean +51.033,
  tumour_best_frac__min -49.931, normal_tau__absdiff -36.698, tumour_mean_log__min -25.514,
  cls_membrane__absdiff -22.527.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
HAVCR2 x PDCD1 - dual checkpoint blockade
-----------------------------------------
STATUS: already documented - lomvastomig (clinical, dual checkpoint). Recovered by the ranking
  rather than proposed by it.
Within-class score 0.826 (percentile among checkpoint|checkpoint pairs, scored as immune
  combination; evidence coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 0.3% of malignant cells co-express both antigens (strongest in
  glioblastoma; measured across 4 dataset(s)). Co-incidence lift 18.60x versus independent
  expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 0.0% carry both at nTPM>=10. Individually the arms peak in
  lymphoma (HAVCR2) and lymphoma (PDCD1).
SAFETY: the AND-gate peak in normal tissue is 1.10 (log1p nCPM); 0 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 4.30 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): t-cells (HAVCR2 32, PDCD1 2 nCPM); brain
  inhibitory neurons (HAVCR2 15, PDCD1 0 nCPM); distal convoluted tubule cells (HAVCR2 13, PDCD1 0
  nCPM); b-cells (HAVCR2 12, PDCD1 0 nCPM).
  Within the tumour itself, double-positive rate is -7.1% higher in malignant than in non-malignant
  cells - i.e. NOT selective for malignant cells in these samples.
MODEL: top additive contributions - cls_membrane__mean +235.855, cls_membrane__min -173.250,
  normal_n_expressed__mean +87.114, normal_max_log__mean -74.209, window_vital_log__mean +59.174,
  taa_tumour_max_log +54.141.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
PDCD1 x TIGIT - dual checkpoint blockade
----------------------------------------
STATUS: already documented - rilvegostomig (clinical, dual checkpoint). Recovered by the ranking
  rather than proposed by it.
Within-class score 0.812 (percentile among checkpoint|checkpoint pairs, scored as immune
  combination; evidence coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 0.4% of malignant cells co-express both antigens (strongest in
  diffuse large B-cell lymphoma; measured across 4 dataset(s)). Co-incidence lift 6.49x versus
  independent expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 0.0% carry both at nTPM>=10. Individually the arms peak in
  lymphoma (PDCD1) and lymphoma (TIGIT).
SAFETY: the AND-gate peak in normal tissue is 1.10 (log1p nCPM); 0 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 3.48 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): t-cells (PDCD1 2, TIGIT 96 nCPM); brain
  inhibitory neurons (PDCD1 0, TIGIT 0 nCPM); distal convoluted tubule cells (PDCD1 0, TIGIT 0
  nCPM); b-cells (PDCD1 0, TIGIT 4 nCPM).
  Within the tumour itself, double-positive rate is -2.9% higher in malignant than in non-malignant
  cells - i.e. NOT selective for malignant cells in these samples.
MODEL: top additive contributions - cls_membrane__mean +238.314, cls_membrane__min -192.412,
  normal_n_expressed__mean +89.707, normal_max_log__mean -80.028, taa_tumour_max_log +51.027,
  window_vital_log__mean +49.788.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
CD96 x TIGIT - dual checkpoint blockade
---------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.794 (percentile among checkpoint|checkpoint pairs, scored as immune
  combination; evidence coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 0.7% of malignant cells co-express both antigens (strongest in
  diffuse large B-cell lymphoma; measured across 4 dataset(s)). Co-incidence lift 5.79x versus
  independent expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 0.1% carry both at nTPM>=10. Individually the arms peak in
  leukemia (CD96) and lymphoma (TIGIT).
SAFETY: the AND-gate peak in normal tissue is 4.57 (log1p nCPM); 3 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 2.15 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): t-cells (CD96 830, TIGIT 96 nCPM); nk-cells
  (CD96 656, TIGIT 59 nCPM); epicardial cells (CD96 43, TIGIT 14 nCPM); plasma cells (CD96 13, TIGIT
  10 nCPM).
  Within the tumour itself, double-positive rate is -5.9% higher in malignant than in non-malignant
  cells - i.e. NOT selective for malignant cells in these samples.
MODEL: top additive contributions - window_vital_log__max -173.220, cls_membrane__min -143.255,
  tumour_best_frac__min -52.034, cls_membrane__mean +49.550, window_vital_log__mean -47.270,
  tumour_mean_log__min -28.814.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
LAG3 x TIGIT - dual checkpoint blockade
---------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.790 (percentile among checkpoint|checkpoint pairs, scored as immune
  combination; evidence coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 0.3% of malignant cells co-express both antigens (strongest in
  diffuse large B-cell lymphoma; measured across 4 dataset(s)). Co-incidence lift 13.30x versus
  independent expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 0.1% carry both at nTPM>=10. Individually the arms peak in
  lymphoma (LAG3) and lymphoma (TIGIT).
SAFETY: the AND-gate peak in normal tissue is 4.03 (log1p nCPM); 2 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 0.55 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): t-cells (LAG3 55, TIGIT 96 nCPM); nk-cells
  (LAG3 11, TIGIT 59 nCPM); endometrial stromal cells (LAG3 4, TIGIT 5 nCPM); epicardial cells (LAG3
  2, TIGIT 14 nCPM).
  Within the tumour itself, double-positive rate is -2.9% higher in malignant than in non-malignant
  cells - i.e. NOT selective for malignant cells in these samples.
MODEL: top additive contributions - cls_membrane__min -181.617, window_vital_log__max -171.290,
  cls_membrane__mean +81.012, normal_max_log__mean -73.082, cls_membrane__absdiff -39.772,
  tumour_best_frac__min -21.423.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

```
CTLA4 x TIGIT - dual checkpoint blockade
----------------------------------------
STATUS: not among the documented bispecific pairs in this project's reference list - i.e. a genuine
  suggestion, with all the caveats below.
Within-class score 0.749 (percentile among checkpoint|checkpoint pairs, scored as immune
  combination; evidence coverage 100.0% of rubric weight).
EFFICACY: in single tumour cells, 0.5% of malignant cells co-express both antigens (strongest in
  diffuse large B-cell lymphoma; measured across 4 dataset(s)). Co-incidence lift 5.13x versus
  independent expression - the two antigens tend to occupy the same cells.
  Across 1,206 cancer cell lines, 0.1% carry both at nTPM>=10. Individually the arms peak in
  skin_cancer (CTLA4) and lymphoma (TIGIT).
SAFETY: the AND-gate peak in normal tissue is 4.48 (log1p nCPM); 1 of 154 normal cell types express
  both above threshold, of which 0 are dose-limiting lineages. Pairing lowers peak normal exposure
  by 0.09 log units versus the worse single arm.
  Normal lineages at risk (highest shared expression): t-cells (CTLA4 87, TIGIT 96 nCPM); nk-cells
  (CTLA4 3, TIGIT 59 nCPM); neutrophil progenitors (CTLA4 3, TIGIT 3 nCPM); innate lymphoid cells
  (CTLA4 2, TIGIT 2 nCPM).
  Within the tumour itself, double-positive rate is -2.3% higher in malignant than in non-malignant
  cells - i.e. NOT selective for malignant cells in these samples.
MODEL: top additive contributions - cls_membrane__mean +237.606, cls_membrane__min -194.969,
  normal_n_expressed__mean +89.447, normal_max_log__mean -75.736, taa_tumour_max_log +51.378,
  window_vital_log__mean +44.781.
CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection an under-
  estimate; cell lines lack microenvironment and often lose lineage markers. This is target triage
  to prioritise experiments, not evidence of efficacy.
```

