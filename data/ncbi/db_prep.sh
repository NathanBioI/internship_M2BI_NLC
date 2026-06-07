#!/bin/bash

# Conda
CONDA_PATH=$(conda info --base)
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate mlst_ospc_pipeline_env

# Paramètres
DIR1="Borreliella"
DIR2="Borrelia"
OUT_FILE="borrelia_db_ncbi.fa"
DB_NAME="borrelia_ncbi_kma_index"

# Nettoyage si un ancien fichier existe
rm -f "$OUT_FILE"

# Fonction pour ajouter les fichiers d'un dossier
process_dir() {
    for filepath in "$1"/*.fa; do
        filename=$(basename "$filepath")
        
        # On vérifie si c'est ospC.fa pour l'ignorer
        if [ "$filename" == "ospC.fa" ]; then
            continue
        fi
        
        # On ajoute le contenu au fichier final
        if [ -f "$filepath" ]; then
            cat "$filepath" >> "$OUT_FILE"
            # saut de ligne entre chaque fichier
            echo "" >> "$OUT_FILE"
        fi
    done
}

# Exécution pour les deux dossiers
process_dir "$DIR1"
process_dir "$DIR2"

# Indexation KMA
if [ -s "$OUT_FILE" ]; then
    kma index -i "$OUT_FILE" -o "$DB_NAME"
    echo "Index KMA généré"
else
    echo "Erreur"
fi
