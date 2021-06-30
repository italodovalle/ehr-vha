# Data


`toy_paths.txt`

* randomly generated disease paths for test runs

`outpat_males_edgelist_unfiltered` -> `outpatient_males_cooccurrence.csv`

* disease pairs and number of total co-occurrence in patients

`bigram_final` -> `outpatient_males_disease_steps.csv`

* number of patients in which the disease in `source` is preceeded by the disease in `target`

`icd9_chapter_mapping.csv`

* ICD-9: codes, descriptions, names, chapters, and sub-chapters

`outpat_males_disease_prevalence.csv`

* prevalence of each disease in the VHA-EHR database


# Methods

## Process Paths

* Input: file containing disease paths
    * each patient is represented by a sequence of ICD codes, separatered by '-'
    * example (toy example, randomly generated):
    
```
650-256-526-667-650-283-101-954-589-817-820-231
508-817-630-626-526-164-965-316-191-757-643-771
458-803-502-630
292-508-458-954
820-141-601-537-879
```

* Output: edge lists
    * each line contains a disease pair and the counts of patients in which the disease i preceded disease j

```
index,count,source,target
650-256,599,650,256
256-526,620,256,526
526-667,607,526,667
667-650,619,667,650
650-283,616,650,283
```

```
usage: process_paths.py [-h] -infile INFILE -b BIGRAM [-t TRIGRAM]

Get edge lists

optional arguments:
  -h, --help      show this help message and exit
  -infile INFILE  infile
  -b BIGRAM       output bigram
  -t TRIGRAM      output trigram (optional)
```

`./process_paths.py -infile ../data/toy_paths.txt -b bigram.csv`


## Phi-coefficient and Fisher Exact Test

Calculating $\phi$ coefficient and its respective p-value


```
usage: compute_fisher_phi.py [-h] [-infile INFILE] [-tmp_dir TMP_DIR] [-outfile OUTFILE] [-ncpus N CPUS]

Compute phi and Fisher

optional arguments:
  -h, --help        show this help message and exit
  -infile INFILE    infile
  -tmp_dir TMP_DIR  tmp dir
  -outfile OUTFILE  outfile
  -ncpus N CPUS     number of cpus
```


* infile: `data/outpatient_males_cooccurrence.csv`
* outfile: `output/fisher/fisher_directional_bh_fdr.csv`


```
./compute_fisher_phi.py -infile ../data/outpatient_males_cooccurrence.csv -tmp_dir tmp/ -ncpus 10

```

## Filter Network

Criteria to obtain the final network:

* Remove diagnoses for which prevalence < 10th percentile
* Remove links with < 100 patients
* Remove links for which a disease is succeeded by a diagnosis that records death.
    * 348: Other conditions of brain (category that contains the diagnosis Brain death)
    * 656: Other known or suspected fetal and placental problems affecting management of mother (catery that contains the diagnosis Intrauterine death affecting management of mother)
    * 768: Intrauterine hypoxia and birth asphyxia
    * 798: Sudden death cause unknown
    
    
Jupyter Notebook: [../notebooks/FilterNetwork.ipynb]
* inputs:
    * `data/outpat_males_disease_prevalence.csv`
    * `data/outpatient_males_disease_steps.csv`
    * `output/fisher/fisher_directional_bh_fdr.csv`
    * `data/icd9_chapter_mapping.csv`
* output:
    * `output/network_filtered/bigram_filtered.csv`
    
    
