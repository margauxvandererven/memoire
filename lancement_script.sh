#!/bin/bash

serveur="student@10.149.144.46"
script_com="/home/student/SPECTRUM/Turbospectrum_NLTE-20.0/COM/script-NLTE-IR_BD-221742_tot.com"

# Connexion SSH et exécution du script .com
echo "Exécution du script sur le serveur distant..."
ssh "$serveur" "cd $(dirname "$script_com"); ./$(basename "$script_com")"

# Commande pour envoyer une notification
osascript -e 'display notification "Votre script a terminé son exécution." with title "Script terminé" sound name "Funk"'

echo "Notification envoyée."


# "Basso"
# "Blow"
# "Bottle"*
# "Frog"
# "Funk"**
# "Glass" (celui que vous utilisez actuellement)
# "Hero"
# "Morse"
# "Ping"
# "Pop"
# "Purr"
# "Sosumi"
# "Submarine"
# "Tink"
