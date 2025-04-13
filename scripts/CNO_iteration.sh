#!/bin/bash

# Configuration
serveur="student@10.149.144.46"
script_OH="/home/student/SPECTRUM/Turbospectrum_NLTE-20.0/COM/script-IR_BD-221742_OH_it.com"
script_CO="/home/student/SPECTRUM/Turbospectrum_NLTE-20.0/COM/script-IR_BD-221742_CO_it.com"
transfert_script="./transfert_spectres.sh"
tolerance=0.05  # Critère de convergence

# Fonction pour créer un range d'abondances
create_abundance_range() {
    local initial_value=$1
    local step=0.1
    local num_points=10
    local values=()
    
    for ((i=-5; i<=5; i++)); do
        values+=("$(bc <<< "$initial_value + $i * $step")")
    done
    echo "${values[@]}"
}

# Fonction pour lancer les scripts distants
run_remote_script() {
    local script=$1
    local wavelength=$2
    local abundance_var=$3
    local abundance_fix=$4
    local iteration=$5
    
    echo "Lancement du script $script avec abondance $abundance_var pour la raie $wavelength"
    ssh "$serveur" "cd $(dirname "$script"); ./$(basename "$script") $wavelength $abundance_var $abundance_fix $iteration" 
}

# Fonction principale d'itération
iterate_CNO() {
    local initial_O=$1
    local initial_C=$2
    local wavelengths_OH=$3
    local wavelengths_CO=$4
    local iteration=1
    
    # Convertir les chaînes de longueurs d'onde en tableaux
    IFS=',' read -ra OH_lines <<< "$wavelengths_OH"
    IFS=',' read -ra CO_lines <<< "$wavelengths_CO"
    
    # while true; do
        echo "Itération $iteration"
        
        # # Création du range pour O
        O_range=($(create_abundance_range $initial_O))
        
        # Boucle sur les abondances de O
        for O_abu in "${O_range[@]}"; do
            # Boucle sur chaque raie OH
            for raie in "${OH_lines[@]}"; do
                run_remote_script "$script_OH" "$raie" "$O_abu" "$initial_C" "$iteration"
            done
        done

        bash "$transfert_script"

        O_range_str=$(IFS=,; echo "${O_range[*]}")
        echo "Abondances O: $O_range_str"
        new_O=$(python3 calculate_O_abundance.py "$iteration" "$O_range_str")
        echo $new_O
        
        # Création du range pour C
        C_range=($(create_abundance_range $initial_C))
        
        # # Boucle sur les abondances de C
        # for C_abu in "${C_range[@]}"; do
        #     # Boucle sur chaque raie CO
        #     for raie in "${CO_lines[@]}"; do
        #         run_remote_script "$script_CO" "$raie" "$C_abu" "$new_O" "$iteration"
        #     done
        # done

        # bash "$transfert_script"

        
        # C_range_str=$(IFS=,; echo "${C_range[*]}")
        # new_C=$(python3 calculate_C_abundance.py "$iteration" "$C_range_str")
        
        # # Vérification de la convergence
        # if (( $(bc <<< "sqrt(($new_O - $initial_O)^2 + ($new_C - $initial_C)^2) < $tolerance") )); then
        #     echo "Convergence atteinte!"
        #     echo "Abondance finale O: $new_O"
        #     echo "Abondance finale C: $new_C"
        #     break
        # fi
        
        # initial_O=$new_O
        # initial_C=$new_C
        # iteration=$((iteration + 1))
        
        # if [ $iteration -gt 10 ]; then
        #     echo "Maximum d'itérations atteint"
        #     break
        # fi
    # done
}

main() {
    local initial_O=8.66 #solaire aspl 2007
    local initial_C=8.44 #shetye 20..
    local wavelengths_OH=$(python3 -c "
import json
with open('../data_lines/raies_moleculaires_test.txt', 'r') as f:
    data = json.load(f)
print(','.join(data['OH_lines'].keys()))
")
    
    local wavelengths_CO=$(python3 -c "
import json
with open('../data_lines/raies_moleculaires_test.txt', 'r') as f:
    data = json.load(f)
print(','.join(data['12C16O_lines'].keys()))
")
    
    iterate_CNO "$initial_O" "$initial_C" "$wavelengths_OH" "$wavelengths_CO"
}

main
