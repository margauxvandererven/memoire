#!/bin/bash

# serveur="student@10.149.144.46"
# dossier_dist="/home/student/SPECTRUM/Turbospectrum_NLTE-20.0/COM/syntspec/"
# # dossier_dist="/home/student/SPECTRUM/TSFitPy/output_files"
# dossier_local="/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/syntspec/BD-221742b/"

LOCAL_FOLDER="/Users/margauxvandererven/OneDrive - Université Libre de Bruxelles/memoire/BD-221742"
REMOTE_USER="student"
REMOTE_HOST="10.149.144.46"
REMOTE_FOLDER="/home/student/SPECTRUM/Turbospectrum_NLTE-20.0/COM/syntspec/synth_margaux"

rsync -avz "$REMOTE_USER@$REMOTE_HOST:$REMOTE_FOLDER" "$LOCAL_FOLDER/" 

if [ $? -eq 0 ]; then
    echo "Synchronisation réussie."
else
    echo "Erreur lors de la synchronisation."
fi


# echo "Téléchargement des fichiers générés récemment..."
 
# read -p "Combien de fichiers récents souhaitez-vous télécharger ? " nombre_fichiers

# # Vérifier que l'entrée est un nombre valide
# if ! [[ "$nombre_fichiers" =~ ^[0-9]+$ ]]; then
#     echo "Veuillez entrer un nombre valide."
#     exit 1
# fi

# fichiers_recents=$(ssh "$serveur" "ls -td $dossier_dist/*Oabu_8.54.conv | head -n $nombre_fichiers")
# for fichier in $fichiers_recents; do
#     scp -v "$serveur:$fichier" "$dossier_local"
# done

# echo "Les spectres récents ont été téléchargés avec succès."

# #Définir l'adresse du serveur et les répertoires
# serveur="student@10.149.144.46"
# dossier_dist="/home/student/SPECTRUM/TSFitPy/output_files"
# dossier_local="/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/syntspec/BD-221742b/syntspec_TSFit/"

# # Trouver le dernier dossier créé
# dernier_dossier=$(ssh "$serveur" "ls -td $dossier_dist/*/ | head -n 1")

# # Vérifier si un dossier a été trouvé
# if [ -n "$dernier_dossier" ]; then
#     echo "Dernier dossier trouvé : $dernier_dossier"

#     # Télécharger le dossier en entier avec SCP
#     scp -r -v "$serveur:$dernier_dossier" "$dossier_local"

#     echo "Le dernier dossier a été téléchargé avec succès."
# else
#     echo "Aucun dossier récent trouvé dans $dossier_dist."
# fi