# Pipeline Borrelia MLST / NCBI / ospC – projet GENAspe

Ce dépôt contient un pipeline bioinformatique pour analyser des données de séquençage de *Borrelia*, Le pipeline combine :

- un typage MLST à partir des allèles PubMLST ;
- une seconde identification MLST à partir d'une base NCBI locale ;
- une analyse phylogénétique MLST ;
- une analyse du gène *ospC* ;
- des analyses statistiques et des figures finales en R.

Le pipeline est piloté par `main.sh` depuis la racine du dépôt.

---

## 1. Organisation générale du dépôt

Structure attendue :

```text
.
├── main.sh
├── scripts/
│   ├── qc.sh
│   ├── loci_specificity.py
│   ├── typing-mapping_with_kma.py
│   ├── typing_kma_ncbi.py
│   ├── build_final_species_table.py
│   ├── make_out_table.py
│   ├── concatenate_loci.py
│   ├── ospc_build_candidates.sh
│   ├── ospc_pipeline_helpers.py
│   ├── run_ospc_selection_and_tree.sh
│   ├── select_ospc_for_phylogeny.py
│   ├── analyse_stats_genASPE.R
│   ├── analyse_stats_genASPE_posthoc.R
│   ├── repartition_borrelia_year.R
│   ├── repartition_borrelia_altitude.R
│   ├── tree_MLST_heatmap_alt_year.R
│   ├── tree_MLST_loci.R
│   ├── ospc_iNEXT.R
│   ├── tree_ospC.R
│   └── face_to_face_MLST_ospC.R
├── data/
│   ├── fastq_raw/
│   ├── fastq/
│   ├── metadata_genaspe.csv
│   ├── ospC_ena_myannot.fasta
│   ├── pubMLST_alleles/
│   │   ├── raw/
│   │   │   ├── clpA.fas
│   │   │   ├── clpX.fas
│   │   │   ├── nifS.fas
│   │   │   ├── pepX.fas
│   │   │   ├── pyrG.fas
│   │   │   ├── recG.fas
│   │   │   ├── rplB.fas
│   │   │   └── uvrA.fas
│   │   ├── alleles.fasta
│   │   └── borrelia_burgdorferi_sensu_lato.*
│   ├── pubMLST_profile/
│   │   ├── BIGSdb_3429203_6349790279_07966.csv
│   │   └── borrelia_spp
│   └── ncbi/
│       ├── Borrelia/extract_borrelia.sh
│       ├── Borreliella/extract_borreliella.sh
│       └── db_prep.sh
└── results/
```

Les dossiers `results/` et `data/fastq/` peuvent être créés automatiquement par le pipeline. Le dossier `data/fastq_raw/` doit contenir les FASTQ bruts.

---

## 2. Environnements conda

Le projet utilise trois environnements conda :

```text
qc_env                    FastQC + MultiQC
trimming_env              Cutadapt
mlst_ospc_pipeline_env    Python, KMA, seqkit, MAFFT, IQ-TREE, BBTools, SPAdes, HybPiper, BLAST, BWA, R et packages R
```

Les fichiers YAML fournis peuvent être utilisés pour les créer :

```bash
mamba env create -f env/qc_env.yml
mamba env create -f env/trimming_env.yml
mamba env create -f env/mlst_ospc_pipeline_env.yml
```

Si `mamba` n'est pas disponible :

```bash
conda env create -f env/qc_env.yml
conda env create -f env/trimming_env.yml
conda env create -f env/mlst_ospc_pipeline_env.yml
```

---

## 3. Données d'entrée attendues

### 3.1 FASTQ bruts

Les FASTQ bruts appariés doivent être placés dans :

```text
data/fastq_raw/
```

Nommage attendu :

```text
sample_R1.fastq.gz
sample_R2.fastq.gz
```

Après QC/trimming, les fichiers filtrés sont écrits dans :

```text
data/fastq/
```

Ces fichiers trimés sont ensuite utilisés par le pipeline.

### 3.2 Métadonnées

Le fichier de métadonnées doit être :

```text
data/metadata_genaspe.csv
```

Il doit contenir les colonnes permettant de relier les échantillons aux sites et aux dates de collecte :

```text
n°
n°ADN EPIA
Type échantillon
date collecte
```

### 3.3 Références PubMLST

Les allèles MLST bruts doivent être placés dans :

```text
data/pubMLST_alleles/raw/
```

avec les huit loci :

```text
clpA.fas
clpX.fas
nifS.fas
pepX.fas
pyrG.fas
recG.fas
rplB.fas
uvrA.fas
```

Si l'index KMA PubMLST est absent, `main.sh` il est reconstruit automatiquement.
Les fichiers attendus après indexation sont :

```text
borrelia_burgdorferi_sensu_lato.name
borrelia_burgdorferi_sensu_lato.comp.b
borrelia_burgdorferi_sensu_lato.length.b
borrelia_burgdorferi_sensu_lato.seq.b
```

### 3.4 Profils PubMLST / BIGSdb

Le fichier BIGSdb utilisé pour relier allèles, ST et espèces doit être :

```text
data/pubMLST_profile/BIGSdb_3429203_6349790279_07966.csv
```

Le fichier des profils ST doit être :

```text
data/pubMLST_profile/borrelia_spp
```

### 3.5 Base NCBI locale

La base NCBI est préparée à partir de :

```text
data/ncbi/Borrelia/extract_borrelia.sh
data/ncbi/Borreliella/extract_borreliella.sh
data/ncbi/db_prep.sh
```

Si l'index KMA NCBI est absent, `main.sh` le lance automatiquement :
L'index attendu est :

```text
data/ncbi/borrelia_ncbi_kma_index.*
```

### 3.6 Références ospC

L'analyse ospC utilise par défaut :

```text
data/ospC_ena_myannot.fasta
```

Il est possible de changer ce fichier au lancement :

```bash
OSPC_REF_FASTA=data/ospC_refs.fasta bash main.sh
```

---

## 4. Lancement du pipeline

### 4.1 Pipeline complet

Depuis la racine du dépôt :

```bash
bash main.sh
```

Le pipeline complet lance :

```text
1. QC, trimming et QC post-trim
2. indexation PubMLST KMA si nécessaire
3. spécificité des allèles PubMLST
4. typage PubMLST par KMA
5. préparation de la base NCBI si nécessaire
6. typage complémentaire NCBI par KMA
7. construction du tableau final des espèces
8. table présence/absence pour les analyses de prévalence
9. alignements MLST avec MAFFT
10. concaténation des loci MLST
11. arbre MLST avec IQ-TREE
12. analyse ospC
13. scripts R de statistiques et figures
```

### 4.2 Lancement sans ospC

Pour désactiver toute la partie ospC :

```bash
bash main.sh -no_ospc
```

Dans ce mode, le pipeline ignore :

```text
scripts/ospc_build_candidates.sh
scripts/run_ospc_selection_and_tree.sh
scripts/ospc_iNEXT.R
scripts/tree_ospC.R
scripts/face_to_face_MLST_ospC.R
```

Les analyses MLST, NCBI, statistiques, figures de répartition et arbres MLST restent exécutées.

### 4.3 Variables utiles

variables modifiables au lancement :

nombre de threads :
```bash
THREADS=8 bash main.sh
```

références *osp*C
```bash
OSPC_REF_FASTA=data/ospC_refs.fasta bash main.sh
```

environnement conda principale :
```bash
ENV_NAME=mlst_ospc_pipeline_env bash main.sh
```

---

## 5. Sorties

### 5.1 Typage PubMLST

```text
results/pubMLST_typing/<sample>/alleles.tsv # locus correspondant par échantillon
results/pubMLST_typing/<sample>/st_sp.tsv # espèce et ST
results/pubMLST_typing/<sample>/perfect_alleles_specificity.tsv  # loci spécifique présent pour l'échantillon
results/pubMLST_typing/mafft_input/<locus>.fasta  # séquence des loci de l'échantillon
```

### 5.2 Typage NCBI

```text
results/ncbi_typing/<sample>/<sample>.res # fichier résultats de l'outil KMA
results/ncbi_typing/<sample>/<sample>.fsa # séquence des loci de l'échantillon
```

### 5.3 Tableau final

```text
results/metadata_genaspe_species_final.tsv
results/table_infection.tsv
```
`metadata_genaspe_species_final.tsv` contient l'identification finale retenue pour chaque échantillon.
`table_infection.tsv` contient une table présence/absence par espèce, utilisée par les scripts statistiques et graphiques.

### 5.4 Phylogénie MLST

```text 
results/phylogeny/MLST_concat.fasta # séquences utilisées en phylogénie MLST
results/phylogeny/MLST_iqtree.treefile # fichier d'arbre phylogénétique MLST
results/ospc/ospc_selection/selected_ospC_iqtree.treefile  # fichier d'arbre phylogénétique *osp*C
```

### 5.5 Analyses statistiques et figures R

Exemples de sorties :

```text
stats_genASPE_Borrelia/ # analyses statistiques
stats_genASPE_posthoc_bacteries_only/ # analyses statistiques
results/repartition_2023_2024.tsv # graphique à barre de répartition selon l'année
results/repartition_2023_2024.png # graphique à barre de répartition selon l'année
results/repartition_altitude.tsv # graphique à barre de répartition selon l'altitude
results/repartition_altitude.png # graphique à barre de répartition selon l'altitude
tree_altitude_year_heatmaps.png # Arbre MLST annoté avec l'année et l'altitude
tree_MLST_loci.png # Arbre MLST avec loci présent par *Borrelia* 
ospC_tree__A4__heatmap.pdf
face_to_face_MLST_ospC.png # Arbre *MLST*osp*C annoté avec l'année et l'altitude
results/ospc_iNEXT_final/ospc_iNEXT_probabilites_couverture.tsv # analyse de raréfaction
results/ospc_iNEXT_final/figures/ # Courbes de raréfactions
```

---

## 6. Rôle de chaque script

### 6.1 Script principal

#### `main.sh`

Script d'exécution du pipeline.

Il active les environnement conda, lance les étapes dans le bon ordre, gère les index KMA PubMLST et NCBI si absents, puis déclenche les analyses Python, shell et R.

---

### 6.2 QC et trimming

#### `scripts/qc.sh`

Réalise :

```text
1. FastQC sur les FASTQ bruts ;
2. trimming avec Cutadapt ;
3. FastQC sur les FASTQ trimés ;
4. synthèse MultiQC.
```

Les outils sont appelés dans des environnements séparés :

```text
FastQC / MultiQC → qc_env
Cutadapt         → trimming_env
```

Entrées :

```text
data/fastq_raw/*_R1.fastq.gz
data/fastq_raw/*_R2.fastq.gz
```

Sorties :

```text
data/fastq/*_trimmed_R1.fastq.gz
data/fastq/*_trimmed_R2.fastq.gz
results/qc/
```

---

### 6.3 Typage PubMLST

#### `scripts/loci_specificity.py`

Construit le fichier de spécificité allèle/espèce à partir de l'export BIGSdb.
Ce fichier indique, pour chaque couple locus/allèle, l'espèce majoritaire, le nombre d'isolats, la pureté et le détail des espèces observées dans BIGSdb.

#### `scripts/typing-mapping_with_kma.py`

Réalise le typage MLST PubMLST par KMA.

Fonctions :

```text
1. mapping des reads sur la base PubMLST ;
2. extraction du meilleur allèle par locus ;
3. détection des allèles parfaits 100 % identité / 100 % couverture ;
4. recherche des profils ST compatibles ;
5. écriture des FASTA par locus pour la phylogénie MLST ;
```

Entrées :

```text
data/fastq/
data/pubMLST_alleles/borrelia_burgdorferi_sensu_lato.*
data/pubMLST_profile/borrelia_spp
data/pubMLST_profile/mlst_species_specificity/all_alleles_specificity.tsv
```

Sorties :

```text
results/pubMLST_typing/<sample>/alleles.tsv
results/pubMLST_typing/<sample>/st_sp.tsv
results/pubMLST_typing/<sample>/perfect_alleles_specificity.tsv
results/pubMLST_typing/mafft_input/<locus>.fasta
```

---

### 6.4 Base et typage NCBI

#### `data/ncbi/Borrelia/extract_borrelia.sh`

Extrait les séquences de référence du groupe *Borrelia* depuis les données NCBI locales.

#### `data/ncbi/Borreliella/extract_borreliella.sh`

Extrait les séquences de référence du groupe *Borreliella* depuis les données NCBI locales.

#### `data/ncbi/db_prep.sh`

Fusionne les références NCBI extraites et construit l'index KMA :

```text
data/ncbi/borrelia_ncbi_kma_index.*
```

#### `scripts/typing_kma_ncbi.py`

Réalise le mapping contre la base NCBI locale avec KMA.

Seuil de lancement KMA :

```text
-ID 95
```

---

### 6.5 Tableau final et table d'infection

#### `scripts/build_final_species_table.py`

Construit le tableau final d'identification des espèces.

Hiérarchie de décision :

```text
1. PubMLST via les allèles parfaits et leur spécificité d'espèce ;
2. NCBI strict ;
3. NCBI relâché ;
4. inconnu si aucune information exploitable.
```

Entrées :

```text
data/metadata_genaspe.csv
results/pubMLST_typing/
results/ncbi_typing/
```

Sortie :

```text
results/metadata_genaspe_species_final.tsv
```

#### `scripts/make_out_table.py`

Construit une table présence/absence par espèce à partir du tableau final.

Entrée :

```text
results/metadata_genaspe_species_final.tsv
```

Sortie :

```text
results/table_infection.tsv
```

---

### 6.6 Phylogénie MLST

#### `scripts/concatenate_loci.py`

Concatène les alignements MLST locus par locus.

Sorties :

```text
results/phylogeny/MLST_concat.fasta
results/phylogeny/MLST_partitions.txt
```

L'arbre MLST est ensuite produit par `main.sh` avec IQ-TREE :

```bash
iqtree3 \
  -st DNA \
  -s results/phylogeny/MLST_concat.fasta \
  -p results/phylogeny/MLST_partitions.txt \
  -m TESTMERGE \
  -bb 1000 \
  -nt AUTO \
  -pre results/phylogeny/MLST_iqtree \
  -redo
```

Sortie :

```text
results/phylogeny/MLST_iqtree.treefile
```

---

### 6.7 Analyse ospC
#### `scripts/ospc_build_candidates.sh`

Construit des séquebces ospC à partir depuis plusieurs approches d'assemblage et mapping.

Sorties :

```text
results/ospc/ospc_kma/
results/ospc/ospc_seed_spades/
results/ospc/ospc_kma_guided_spades/
results/ospc/hybpiper_ospc/
results/ospc/ospc_personalized/
```

#### `scripts/ospc_pipeline_helpers.py`

Script Python appelé par `ospc_build_candidates.sh`.

Il gère :

```text
- préparation des références ospC ;
- nettoyage des headers FASTA ;
- extraction de la meilleure référence KMA ;
- choix du bait guidé par KMA ;
- construction des graines conservées ;
- extraction / fusion des contigs ;
- construction des bases KMA par échantillon.
```

#### `scripts/run_ospc_selection_and_tree.sh`

Lance la sélection finale des séquences ospC, l'alignement MAFFT et l'arbre IQ-TREE.

Sorties :

```text
results/ospc/ospc_selection/ospc_sample_decisions.tsv
results/ospc/ospc_selection/ospc_candidate_summary.tsv
results/ospc/ospc_selection/selected_ospC_consensus_oriented.fasta
results/ospc/ospc_selection/selected_ospC_consensus_oriented.mafft.fasta
results/ospc/ospc_selection/selected_ospC_iqtree.treefile
```

#### `scripts/select_ospc_for_phylogeny.py`

Sélectionne une séquence ospC par échantillon.
Le script produit une table de décision par échantillon et le FASTA final utilisé pour l'arbre ospC.

---

### 6.8 Scripts R statistiques et figures générales

#### `scripts/analyse_stats_genASPE.R`

Réalise les analyses statistiques globales à partir de `results/table_infection.tsv`.

Sortie :

```text
stats_genASPE_Borrelia/
```

#### `scripts/analyse_stats_genASPE_posthoc.R`

Réalise les analyses post-hoc et comparaisons complémentaires à partir des résultats statistiques.

Sortie :

```text
stats_genASPE_posthoc_bacteries_only/
```

#### `scripts/repartition_borrelia_year.R`

Produit un graphique de répartition proportionnelle des espèces entre 2023 et 2024.

Entrée :

```text
results/table_infection.tsv
```

Sorties :

```text
results/repartition_2023_2024.tsv
results/repartition_2023_2024.png
results/repartition_2023_2024.pdf
```

#### `scripts/repartition_borrelia_altitude.R`

Produit un graphique de répartition proportionnelle des espèces selon l'altitude/site.

Entrée :

```text
results/table_infection.tsv
```

Sorties :

```text
results/repartition_altitude.tsv
results/repartition_altitude.png
results/repartition_altitude.pdf
```

#### `scripts/tree_MLST_heatmap_alt_year.R`

Produit un arbre MLST annoté avec :

```text
- espèce ;
- altitude ;
- année.
```

Sorties :

```text
tree_altitude_year_heatmaps.png
tree_altitude_year_heatmaps.pdf
```

#### `scripts/tree_MLST_loci.R`

Produit un arbre MLST avec une heatmap de présence/absence des huit loci MLST.

Sorties :

```text
tree_MLST_loci.png
tree_MLST_loci.pdf
```

---

### 6.9 Scripts R dépendants de ospC

Ces scripts ne sont lancés que si l'analyse ospC est activée.

#### `scripts/ospc_iNEXT.R`

Utilise les sorties ospC, regroupe les séquences ospC en groupes à 95 % d'identité et produit les figures de raréfaction/extrapolation.

Sorties principales :

```text
results/ospc_iNEXT_final/ospc_iNEXT_probabilites_couverture.tsv
results/ospc_iNEXT_final/figures/*.png
results/ospc_iNEXT_final/figures/*.pdf
```

#### `scripts/tree_ospC.R`

Produit un arbre ospC annoté avec les espèces MLST, l'altitude et l'année.

Sorties :

```text
ospC_tree_A4_heatmap.svg
ospC_tree_A4_heatmap.pdf
```

#### `scripts/face_to_face_MLST_ospC.R`

Produit une figure comparant l'arbre MLST et l'arbre ospC face à face, avec des traits reliant les échantillons communs.

Sorties :

```text
face_to_face_MLST_ospC.png
face_to_face_MLST_ospC.jpg
face_to_face_MLST_ospC.svg
face_to_face_MLST_ospC.tsv
```

---

## 7. Nettoyage / relance

Le pipeline réutilise certaines bases déjà indexées si elles existent :

```text
data/pubMLST_alleles/borrelia_burgdorferi_sensu_lato.*
data/ncbi/borrelia_ncbi_kma_index.*
```

Pour forcer leur reconstruction, supprimer les fichiers d'index correspondants avant de relancer.

Exemple pour PubMLST :

```bash
rm data/pubMLST_alleles/borrelia_burgdorferi_sensu_lato.*
bash main.sh
```

Exemple pour NCBI :

```bash
rm data/ncbi/borrelia_ncbi_kma_index.*
bash main.sh
```

---

