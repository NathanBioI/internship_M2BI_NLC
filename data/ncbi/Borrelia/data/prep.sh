#!/bin/bash

# Nom du fichier
INFO_FILE="data_summary.tsv"

# extraire la col 1 (nom) et la col 6 (identifiant GCF)
# pipe vers une boucle while
tail -n +2 "$INFO_FILE" | awk -F'\t' '{print $1 "|" $6}' | while IFS="|" read -r full_name accession; do
    
    # Nettoyage des espaces de fin de ligne (souvent invisibles)
    accession=$(echo "$accession" | tr -d '\r\n ' )
    full_name=$(echo "$full_name" | xargs) # Enlève les espaces inutiles

    # Si l'accession est vide, on passe
    [[ -z "$accession" ]] && continue

    # Logique pour extraire l'espèce (2ème ou 3ème mot)
    if [[ "$full_name" == "Candidatus"* ]]; then
        species=$(echo "$full_name" | awk '{print $3}')
    else
        species=$(echo "$full_name" | awk '{print $2}')
    fi

    # --- DIAGNOSTIC ---
    # Décommente la ligne suivante si tu veux voir ce que le script lit :
    # echo "Analyse de : $accession | Espèce : $species"

    # suppression (espèce .sp ou sp.)
    if [[ "$species" == "sp." ]] || [[ "$species" == ".sp" ]]; then
        if [[ -d "$accession" ]]; then
            echo "🗑️  Suppression : $accession (genre $species)"
            rm -rf "$accession"
        fi

    # renommage
    elif [[ -d "$accession" ]]; then
        new_name="${species}_${accession}"
        echo "Renommage : $accession -> $new_name"
        mv "$accession" "$new_name"
    fi
done

echo "Terminé"
