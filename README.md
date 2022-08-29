## <b>"Proteomic Characterization of Acute Kidney Injury in Patients Hospitalized with SARS-CoV2 Infection"</b>

##  <b>Citation:</b>

<i>Proteomic Characterization of Acute Kidney Injury in Patients Hospitalized with SARS-CoV2 Infection.</i> 
Ishan Paranjpe, Pushkala Jayaraman, Chen-Yang Su, Sirui Zhou, Steven Chen, Ryan Thompson, Diane Marie Del Valle, Ephraim Kenigsberg, Shan Zhao, Suraj   Jaladanki, Kumardeep Chaudhary, Steven Ascolillo, Akhil Vaid, Arvind Kumar, Manish Paranjpe, Ross O’Hagan, Samir Kamat, Faris F. Gulamali, Hui Xie, Joceyln Harris, Manishkumar Patel, Kimberly Argueta, Craig Batchelor, Kai Nie, Sergio Dellepiane, Leisha Scott, Matthew A Levin, John Cijiang He, Steven G Coca, Lili Chan, Evren U Azeloglu, Eric Schadt, Noam Beckmann, Sacha Gnjatic, Miram Merad, Seunghee Kim-Schulze, Brent Richards, Benjamin S Glicksberg, Alexander W Charney, Girish N Nadkarni 
<i>medRxiv 2021.12.09.21267548; doi: https://doi.org/10.1101/2021.12.09.21267548 </i>

## <b>Code repo Files:</b>

- <i>post_discharge_eGFR.R</i> : Find proteins associated with post-discharge eGFR
- <i>extract_eGFR_measurements.R</i> :  Get all eGFR measurements taken after discharge for patients in the cohort.
- <i>AKI_associated_proteins.R</i>: Find proteins associated with in-hospital AKI. Uses Limma R package
- Data Files: <i>Suhre_somascan.txt</i> and <i>somalogic_interactions.txt</i>
- <i>CreateProteinInteractionNetwork.ipynb</i>: Create protein-protein interaction network for AKI associated proteins. 
  - This file uses the annotations in Suhre_somascan.txt and somalogic_interactions.txt files
