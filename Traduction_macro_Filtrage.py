#%%
"""FONCTION DU SCRIPT
-----------------------
Ce fichier .py est la traduction en Python de la macro VBA 
utilisée dans le fichier Excel fourni en exemple pour effectuer l'opération de filtrage
du spectre (en fréquence) sismique.

Le script inclut une routine d’affichage interactif qui permet de tester
rapidement plusieurs combinaisons de paramètres (sensibilité, zone de
recherche, etc.) et de visualiser immédiatement :
    - l’enveloppe finale filtrée (retirant les plateaux, zones
    linéaires, etc.),
    - le nombre de points conservés par l’opération de filtrage.


PRÉ-REQUIS TECHNIQUES
------------------------
- Avoir Python installé (version 3.x recommandée).
- Avoir les bibliothèques Python suivantes installées :
  - numpy
  - pandas
  - matplotlib
  - openpyxl
  - Un éditeur de code avec shell intégré. Je conseille VS Code, qui propose une fenêtre
  interactive plus ergonomique (clic droit sur le script->"Run current file in interactive window")

COMMENT PRÉPARER LE FICHIER D’ENTRÉE EXCEL (nom par défaut : "IN_Test.xlsx") :
-------------------------
Le fichier doit se trouver dans le même dossier que ce script et contenir au minimum une feuille nommée "A" :
   - La colonne 1 (A) contient les fréquences (en Hz), à partir de la ligne 3.
   - La colonne 2 (B) contient les valeurs d’accélération correspondantes (en m/s² ou en g) ,
     également à partir de la ligne 3.
   - Les lignes 1 et 2 peuvent contenir des titres ou être laissées vides.

Auteur : Laura Brémont - Stagiaire service Calculs NFM Systems (jan-fev 2026)
Dernière version : 10/02/26
"""

"""-------------------------------Partie à compléter------------------------------------"""

nom_fichier_xlsx = "IN_Test.xlsx" #doit se trouver dans le même répertoire que ce script (xls(x) acceptés)

#/!\ Assurez-vous que le fichier n'est pas ouvert avant de lancer le script

"""-------------------------------------------------------------------------------------"""


import pandas as pd
from openpyxl import load_workbook
# Pour conserver les macros, utiliser openpyxl ou xlwings 
# uniquement, sans pandas pour écrire les données directement
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 14


###     Définitions utiles
"""le logarithme naturel (base e) est np.log(x, le décimal np.log10(x)"""
L_Start = 1
Nb_passe = "Nombre d'itérations"   
CVX_taille = "Taille de la zone de recherche des segments convexes"  
LINEAIRE_taille = 'Taille de la zone de recherche des segments linéaires'
PLATEAU_taille = 'Taille de la zone de recherche des plateaux'  
Sigma_X = "Sensibilité de sélection des plateaux"   
Sigma_D_X = "Sensibilité de sélection des zones linéaires"

nom_X = 'X'
unite_X = 'm/s^2'

couleur = ['maroon','crimson','gold','magenta','blue','white']

parametres = {LINEAIRE_taille : "de la taille de la zone de recherche des segments linéaires",
              CVX_taille : "de la taille de la zone de recherche des segments convexes",
              Sigma_X : "de la sensibilité de sélection des plateaux",
              Sigma_D_X : "de la sensibilité de sélection des zones linéaires",
              PLATEAU_taille : "de la taille de la zone de recherche des plateaux",
              Nb_passe : "du nombre d'itérations"}

defaut = {LINEAIRE_taille : 6,
              PLATEAU_taille : 4,
              CVX_taille : 5,
              Sigma_X : 0.01,
              Sigma_D_X : 4,
              Nb_passe : 3}

valeurs = {Nb_passe : [1,2,3,4,5,6],                                # on fera 20 itérations du processus de sélection de point 
                CVX_taille : [3,4,5,6,7,8],                         # plage de recherche des parties cvx
                LINEAIRE_taille : [3,4,5,6,7,8],                    # plage de recherche des parties linéaires
                PLATEAU_taille : [3,4,5,6,7,8],                     # plage de recherche des plateaux
                Sigma_X : [0.001, 0.005, 0.01, 0.02, 0.03, 0.04],   # sensibilité de détection des plateaux
                Sigma_D_X : [0.5, 1, 2,2.5,3,3.5]}                  # sensibilité de détection des parties cvx


def Run_Enveloppe(nom_fichier_xlsx , L_Start, choix):
    ###     Import des donées via pandas 
    """ajouter qqch pour choisir la table"""
    WS_Spectres = pd.read_excel(nom_fichier_xlsx, usecols="A,B,D,E")
    
    FREQ1= WS_Spectres[WS_Spectres.columns[0]]
    X1= WS_Spectres[WS_Spectres.columns[1]]
    FREQ = FREQ1.values
    X = X1.values

    RES_FREQ = np.log(WS_Spectres[WS_Spectres.columns[2]])
    RES_X = WS_Spectres[WS_Spectres.columns[3]]
    
    NB_val = len(FREQ)

    """
    V0 est l'historique des sélections : en abscisses, la valeur X de l'échantillon
    fourni pas le client, en ordonnée l'état de sélection (0 ou 1) à l'étape N.
    
    Les colonnes de sorties seront donc FREQ et V0[Nb_passe + 1]*X ;
    pour le tracé, deux solutions : 
        - copier les données en éliminant les points d'ordonnée zéros
        - tracer un diagremmes en barres (histogramme ou spectre de raies)
    """
    V0 = np.ones([1 + int(choix[Nb_passe]), NB_val])  

    for N in range(1, int(choix[Nb_passe]) + 1):    # Début de l'étape de sélection n°N
        
        #print("----------ETAPE %d----------"%N)
        V0[N] = V0[N-1]
        ##  Protocole de sélection des points

        """-----------------------------DETECTION PLATEAUX-----------------------------------"""
        ##  Vérification de la présence d'un plateau sur Plateau_taille éléments encore sélectionnés
        
        for j in range(NB_val - int(choix[PLATEAU_taille])):
            moy = 0
            sigma = 0
            
            marker = np.array([0]*NB_val) #pour stoker les abscisses des points du plateau
            curs = j    #curseur de recherche des PLATEAU_taille prochains éléments "encore en jeu"
            
            while curs >= 0 and curs <= NB_val and (curs-j) < choix[PLATEAU_taille]:
                       
                if V0[N-1, curs] == 1:   #l'élément est encore 'dans la course'
                    #il compte dans le calcul de la moyenne et est sauvegardé pour l'écart type
                        moy += X[curs]/choix[PLATEAU_taille]
                        marker[curs] = 1

                curs +=1
            marker[curs] = 1
                                    
            sigma = np.sum(((X - moy)*marker)**2)
            
            if (sigma/choix[PLATEAU_taille])**0.5 < choix[Sigma_X] and sum(marker) >= 3:   #les données sont suffisement proches pour être redondantes
                #print("Plateau détecté à f = {}".format(j))
                
                #mettre les données intermédiaires à 0 dans le tableau V0
                f = j + 1
                compt = choix[PLATEAU_taille]
                while compt > 1:
                    if marker[f] == 1:
                        V0[N, f] = 0
                        compt -= 1
                    f += 1
                    
            
        """-----------------------------DETECTION PARTIE LINEAIRE-----------------------------------"""
        ## Vérification de la présence d'une partie linéaire
        
        #Calcul dérivée prmière en i
        D_FREQ = np.gradient(X, FREQ)
        
        for j in range(1, NB_val - int(choix[LINEAIRE_taille])-1):
            
            moy_D = 0
            sigma_D = 0
            
            marker_D = np.array([0]*NB_val) #pour stoker les abscisses des points du plateau
            curs_D = j    #curseur de recherche des PLATEAU_taille prochains éléments "encore en jeu"
            
            while curs_D >= 0 and curs_D <= NB_val and (curs_D - j) < choix[LINEAIRE_taille]:
                       
                if V0[N-1, curs_D] == 1:   #l'élément est encore 'dans la course'
                    #il compte dans le calcul de la moyenne et est sauvegardé pour l'écart type
                        moy_D += D_FREQ[curs_D]/choix[LINEAIRE_taille]
                        marker_D[curs_D] = 1
                curs_D +=1
            marker_D[curs_D] = 1
                                    
            sigma_D = np.sum(((X - moy_D)*marker_D)**2)

            if (sigma_D/choix[LINEAIRE_taille])**0.5 < choix[Sigma_D_X] and sum(marker_D) >= 3:   #les données sont suffisement linéaires pour être redondantes
                #print("Linéarité détectée à f = {}".format(j))
                
                #mettre les données intermédiaires à 0 dans le tableau V0
                g = j + 1
                compt_D = choix[LINEAIRE_taille]
                while compt_D > 1 and g < NB_val:
                    if marker_D[g] == 1:
                        V0[N, g] = 0
                        compt_D -= 1
                    g += 1        
        
        """-----------------------------DETECTION PARTIE CONVEXE-----------------------------------"""
        ## Vérification de la présence d'une partie convexe
            
        #Calcul dérivée seconde en i
        D2_FREQ = np.gradient(D_FREQ)
        
        for i in range(1, NB_val - int(choix[CVX_taille])):
            est_cvx = (D2_FREQ[i] >= 0)
            f = i + 1
            compt = 0
            while est_cvx and (compt < choix[CVX_taille]):
                #on sort si la partie n'est pas cvx ou qu'elle l'est, mais excède la sensibilité donnée
                est_cvx = (D2_FREQ[f] >= 0)
                compt += 1
                f += 1
                
            if est_cvx: #si le segment de taille CVX_taille est convexe, on supprime les données intermédiaires
                for r in range(i + 1, i + int(choix[CVX_taille]) + 1):
                    V0[N, r] = 0                
                 

    ###     Stockage
    
    RES_X = [V0[int(choix[Nb_passe])][i]*X[i] for i in range(NB_val) if V0[int(choix[Nb_passe])][i]*X[i] != 0]
    RES_FREQ = [FREQ[i] for i in range(NB_val) if V0[int(choix[Nb_passe])][i]*X[i] != 0]
    
    return FREQ, X, RES_FREQ, RES_X, sum(V0[int(choix[Nb_passe])])



def Display_Choices(PARAMETRE, choix):
    
    print("Choix {}".format(parametres[PARAMETRE]))
    
    plt.figure()
    for n in range(6):
        choix[PARAMETRE] = valeurs[PARAMETRE][n]
        FREQ, X, RES_FREQ, RES_X, nb_pts_restants = Run_Enveloppe(nom_fichier_xlsx, L_Start, choix)
        
        plt.scatter(RES_FREQ, RES_X, color=couleur[n], marker='o', s=3.2**(6-n), label='{} -> {} pts'.format(valeurs[PARAMETRE][n], nb_pts_restants))
        
    plt.title("Test pour le paramètre : {}".format(PARAMETRE))
    #Affichage comparatif
    plt.plot(FREQ, X, color='orange', label='Données client')
    plt.xscale('log')
    plt.xlabel('Fréquence (Hz)')
    plt.ylabel('{} ({})'.format(nom_X, unite_X))
    plt.legend(loc='best')
    plt.show()
    
    res = float(input("Valeur la plus adaptée : "))
    return res
    
    
    
def main():
    choix = defaut  #vont être mis à jour au fur et à mesure de l'affinement du filtre
    
    for PARAMETRE in parametres.keys():
        
        choix[PARAMETRE] = Display_Choices(PARAMETRE, choix)
        
    FREQ, X, RES_FREQ, RES_X, nb_points = Run_Enveloppe(nom_fichier_xlsx, L_Start, choix)
    
    print("Pensez à fermer toutes les fenêtres pour déclencher l'écriture des résultats.")
    
    
    #Affichage comparatif final
    
    plt.title("Tracé de l'enveloppe avec les paramètres retenus")
    plt.plot(FREQ, X, color='orange', label='Données client')
    print("Données client affichées")
    plt.plot(RES_FREQ, RES_X, color='indigo', linestyle = '--', label='Données filtrées : {} pts'.format(nb_points))
    print("Données filtrées affichées")
    plt.xscale('log')
    plt.xlabel('Fréquence (Hz)')
    plt.ylabel('{} ({})'.format(nom_X, unite_X))
    plt.legend(loc='best')
    plt.show()
    

    return {'FREQ RES': RES_FREQ[1:], 
            'X RES' : RES_X[1:]}
    
  
        
"""------------------Gestion de l'affichage et de l'import/export des données--------------------""" 
 
    
df_acc = pd.read_excel(nom_fichier_xlsx, sheet_name="A", engine="openpyxl")
accel_data = df_acc.iloc[2:, 1].values  # Colonne 2 (les données commencent à la ligne 3)
freq = df_acc.iloc[2:, 0].values  # Colonne 1 (sans la ligne de titre)

# Lancer le calcul (fonction supposée existante)
res = main()  # Résultat attendu sous la forme d'une liste ou d'un tableau


# Charger le fichier existant avec openpyxl pour modification sans écrasement

try:
    # Charger le fichier Excel
    wb = load_workbook(nom_fichier_xlsx)
    
    # Créer/Réinitialiser la feuille "RES"
    if "RES" in wb.sheetnames:
        ws = wb["RES"]
        ws.delete_rows(1, ws.max_row)  # Vider les anciennes données
    else:
        ws = wb.create_sheet(title="RES")
    
    # Ajouter les en-têtes
    headers=["RES FREQ","RES X"]
    for col_idx, header in enumerate(headers, start=1):
        ws.cell(row=1, column=col_idx, value=header)
    
    # Ajouter les données
    for row_idx, row in enumerate(res.values(), start=1):  # Commence à la ligne 2
        for col_idx, value in enumerate(row, start=2):
            ws.cell(row=col_idx, column=row_idx, value=value)   #avoir 2 colonnes (feuille RES) à partir de 2 lignes (dict res)
    
    # Sauvegarder les modifications
    wb.save(nom_fichier_xlsx)
    wb.close() 
    print("Résultats écrits avec succès ! Vous pouvez ouvrir {}".format(nom_fichier_xlsx))
   
     
except Exception as e:
    print(f"Une erreur s'est produite : {e}")
