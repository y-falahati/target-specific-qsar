# Data Dictionary — EGFR ChEMBL Activity Dataset

**File:** `data/raw/chembl_egfr_activities.json`
**Target:** EGFR (`CHEMBL203`)
**Records:** 58,847
**Columns:** 47

This data dictionary describes the structure and role of each field in the raw ChEMBL EGFR activity dataset. The raw dataset is kept unchanged; decisions about filtering, transformation, and exclusion are deferred to the Data Curation stage.

| #  | Column                      | Meaning / Description                                                                    | Current Data Type | Missing |
| -- | --------------------------- | ---------------------------------------------------------------------------------------- | ----------------- | ------: |
| 1  | `action_type`               | Information about the type/nature of action associated with the activity                 | object / dict     |  89.97% |
| 2  | `activity_comment`          | Textual comments describing the activity result                                          | string            |  55.00% |
| 3  | `activity_id`               | Unique identifier of the ChEMBL activity record                                          | integer           |      0% |
| 4  | `activity_properties`       | Additional properties associated with the activity record                                | list              |      0% |
| 5  | `assay_chembl_id`           | Unique ChEMBL identifier of the assay                                                    | string            |      0% |
| 6  | `assay_description`         | Description of the assay and experimental context                                        | string            |      0% |
| 7  | `assay_type`                | High-level assay classification, e.g. B = Binding, F = Functional                        | string            |      0% |
| 8  | `assay_variant_accession`   | Accession identifier associated with the protein/variant used in the assay               | string            |  80.73% |
| 9  | `assay_variant_mutation`    | Reported EGFR mutation used in the assay, e.g. `L858R` or `T790M`                        | string            |  76.28% |
| 10 | `bao_endpoint`              | BioAssay Ontology endpoint identifier                                                    | string            |      0% |
| 11 | `bao_format`                | BioAssay Ontology format identifier                                                      | string            |      0% |
| 12 | `bao_label`                 | Descriptive BioAssay Ontology label, e.g. `single protein format` or `cell-based format` | string            |      0% |
| 13 | `canonical_smiles`          | Canonical SMILES representation of the compound structure                                | string            |   0.18% |
| 14 | `data_validity_comment`     | Comment indicating a potential data-validity issue                                       | string            |  98.84% |
| 15 | `data_validity_description` | Description of a potential data-validity issue                                           | string            |  98.84% |
| 16 | `document_chembl_id`        | ChEMBL identifier of the source document/publication                                     | string            |      0% |
| 17 | `document_journal`          | Journal in which the source publication appeared                                         | string            |  50.61% |
| 18 | `document_year`             | Publication year of the source document                                                  | float             |  10.31% |
| 19 | `ligand_efficiency`         | Ligand-efficiency-related properties associated with the activity                        | object / dict     |  66.26% |
| 20 | `modality`                  | Modality/type classification of the molecule or agent                                    | string            |  97.54% |
| 21 | `molecule_chembl_id`        | Unique ChEMBL identifier of the molecule                                                 | string            |      0% |
| 22 | `molecule_pref_name`        | Preferred/known name of the molecule                                                     | string            |  81.09% |
| 23 | `parent_molecule_chembl_id` | ChEMBL identifier of the parent molecule                                                 | string            |      0% |
| 24 | `pchembl_value`             | Normalized activity value expressed on the pChEMBL scale                                 | object / string   |  63.73% |
| 25 | `potential_duplicate`       | Flag indicating a potentially duplicated activity record                                 | integer           |      0% |
| 26 | `qudt_units`                | Unit identifier represented according to QUDT                                            | string            |  28.30% |
| 27 | `record_id`                 | Identifier of the underlying data record                                                 | integer           |      0% |
| 28 | `relation`                  | Relationship between the reported activity and its value, e.g. `=`, `>`, `<`             | string            |   7.89% |
| 29 | `src_id`                    | Identifier of the data source within ChEMBL                                              | integer           |      0% |
| 30 | `standard_flag`             | Flag associated with ChEMBL activity standardization                                     | integer           |      0% |
| 31 | `standard_relation`         | Relationship associated with the standardized activity value, e.g. `=`, `>`, `<`         | string            |   7.89% |
| 32 | `standard_text_value`       | Standardized textual activity value when a numerical value is not available              | string            |  94.04% |
| 33 | `standard_type`             | Standardized activity endpoint, e.g. `IC50`, `Ki`, `Kd`, `Inhibition`                    | string            |      0% |
| 34 | `standard_units`            | Standardized unit of the activity measurement, e.g. `nM`, `uM`, `%`                      | string            |  16.39% |
| 35 | `standard_upper_value`      | Upper bound of the standardized activity value for interval/censored measurements        | object            |    100% |
| 36 | `standard_value`            | Standardized numerical activity value; a key candidate for the QSAR target variable      | object / string   |  31.02% |
| 37 | `target_chembl_id`          | ChEMBL identifier of the biological target; `CHEMBL203` in this project                  | string            |      0% |
| 38 | `target_organism`           | Organism associated with the biological target; `Homo sapiens` in this project           | string            |      0% |
| 39 | `target_pref_name`          | Preferred name of the biological target; EGFR in this project                            | string            |      0% |
| 40 | `target_tax_id`             | Taxonomy identifier of the target organism                                               | string            |      0% |
| 41 | `text_value`                | Original activity value represented as text                                              | string            |  94.04% |
| 42 | `toid`                      | Identifier associated with ontology/target mapping                                       | object            |    100% |
| 43 | `type`                      | Original activity type in the source data                                                | string            |      0% |
| 44 | `units`                     | Original activity measurement unit                                                       | string            |  10.96% |
| 45 | `uo_units`                  | Unit identifier according to the Unified Ontology                                        | string            |  28.29% |
| 46 | `upper_value`               | Upper bound of the original activity value for interval/censored measurements            | string            |  99.80% |
| 47 | `value`                     | Original activity value before/alongside standardization                                 | string            |  31.02% |

## Key Fields for the Project

### Compound Identity

* `molecule_chembl_id`
* `parent_molecule_chembl_id`
* `canonical_smiles`

These fields identify the chemical entity and provide the molecular structure required for structure-based analysis.

### Activity

* `standard_type`
* `standard_value`
* `standard_units`
* `standard_relation`
* `pchembl_value`

These fields describe the biological activity measurement and will be central to defining the modeling endpoint during Data Curation.

### Assay Context

* `assay_chembl_id`
* `assay_description`
* `assay_type`
* `bao_format`
* `bao_label`

These fields provide experimental context and are important because activity values can originate from substantially different assay formats.

### Mutation Context

* `assay_variant_accession`
* `assay_variant_mutation`

These fields describe whether a specific target variant or mutation was reported for the assay.

**Important:** A missing `assay_variant_mutation` value means that no mutation was reported in this field. It should not automatically be interpreted as proof that the assay used wild-type EGFR.

### Target Information

* `target_chembl_id`
* `target_pref_name`
* `target_organism`
* `target_tax_id`

These fields identify and describe the biological target.

### Quality Control

* `data_validity_comment`
* `data_validity_description`
* `potential_duplicate`
* `standard_flag`
* `relation`

These fields can be used during quality assessment and data curation.

## Data-Curation Notes

The raw dataset should remain unchanged at this stage.

Missing values should **not** be handled with a blanket `dropna()` operation. Different fields have different meanings and therefore require different decisions during curation.

In particular:

* Missing `standard_value` means there is no numerical activity value available in that field.
* Missing `standard_units` prevents straightforward numerical standardization unless another valid quantitative field can be used.
* Missing `canonical_smiles` prevents structure-based QSAR for that record.
* Missing `pchembl_value` does not necessarily invalidate a record if a valid standardized activity value and unit are available.
* Missing `assay_description` should be treated as missing assay-context information rather than automatically removing the record.
* Missing `assay_variant_mutation` should not automatically be classified as wild-type EGFR.

Fields such as `standard_value`, `value`, and `pchembl_value` are represented as strings/objects in the raw JSON-derived dataframe even when their contents are numerical. Type conversion should therefore be handled explicitly during Data Curation.

**Current project stage:**

`Raw ChEMBL Data → Data Dictionary → Initial EDA → Data Curation → QSAR Dataset → Modeling`

The current Data Dictionary documents the raw dataset only. No filtering or modeling decision is encoded in this document.
