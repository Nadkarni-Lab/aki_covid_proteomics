## <b>"Proteomic Characterization of Acute Kidney Injury in Patients Hospitalized with SARS-CoV2 Infection"</b>
###### <b>Code repo</b>

<b>Files:</b>

- <i>post_discharge_eGFR.R</i> : Find proteins associated with post-discharge eGFR
- <i>extract_eGFR_measurements.R</i> :  Get all eGFR measurements taken after discharge for patients in the cohort.
- <i>AKI_associated_proteins.R</i>: Find proteins associated with in-hospital AKI. Uses Limma R package
- <i>CreateProteinInteractionNetwork.ipynb</i>: Create protein-protein interaction network for AKI associated proteins. 
  - This file uses the annotations in Suhre_somascan.txt and somalogic_interactions.txt files
