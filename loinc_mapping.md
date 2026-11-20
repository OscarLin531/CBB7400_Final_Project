# LOINC Mapping

LOINC codes for the two lab tests in this project.

| Test | NHANES variable | LOINC | Long Common Name | System | Example UCUM Units |
|---|---|---|---|---|---|
| HbA1c | `LBXGH` | [4548-4](https://loinc.org/4548-4/) | Hemoglobin A1c/Hemoglobin.total in Blood | Bld (blood) | % |
| Fasting plasma glucose | `LBXGLU` → `fpg_harmonized` | [1558-6](https://loinc.org/1558-6/) | Fasting glucose [Mass/volume] in Serum or Plasma | Ser/Plas (serum or plasma) | mg/dL |

**Note:** `fpg_harmonized` applies the CDC bridging equation to 2013–2014 values; later cycles use `LBXGLU` as published.

**Source:** LOINC codes and names from loinc.org (https://loinc.org/4548-4/, https://loinc.org/1558-6/; accessed November 19, 2026).

This material contains content from LOINC® (http://loinc.org). LOINC is copyright © Regenstrief Institute, Inc. and the Logical Observation Identifiers Names and Codes (LOINC) Committee and is available at no cost under the license at http://loinc.org/license.
