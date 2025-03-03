#!/bin/bash

serveur="student@10.149.144.46"

#faire synthèse avec et sans mol, puis utiliser find peaks et trouver des raies OH, 12CO, 13CO, 12CN, 13CN
#faire un fichier config avec raie et abu de O
#lancer le fichier script-OH, faire la minimisation avec script python et sortir abu de O
#dans fichier config faire varier abu de C, garder abu de 


script_com="/home/student/SPECTRUM/Turbospectrum_NLTE-20.0/COM/script-IR_BD-221742_OH.com"

# Connexion SSH et exécution du script .com
echo "Exécution du script sur le serveur distant..."
ssh "$serveur" "cd $(dirname "$script_com"); ./$(basename "$script_com")"