"""Composition du programme :
-   Fonction sro_nigam
-   Partie lecture d'Excel, appel, et écriture des résultats

Auteur : Laura Brémont
Dernière version : 29/01/26
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 14

###Fréquences de l'O.H à tester
frequencies = np.logspace(-1, 1.5, num=150)




def sro_nigam(frequencies, accelerogram, dt, damping):
    """
    Calcule le Spectre de Réponse d'Oscillateur (SRO) selon la méthode de NIGAM & JENNINGS.
    
    Entrées :
    frequencies : array
        Tableau des fréquences propres (Hz) pour lesquelles calculer le spectre
    accelerogram : array
        Signal d'accélération (accélérogramme)
    dt : float
        Pas de temps du signal (s)
    damping : float
        Amortissement réduit xi
        

    Sortie : dict
        Dictionnaire contenant :
        - 'D': Spectre de déplacement relatif
        - 'V': Spectre de pseudo-vitesse
        - 'A': Spectre de pseudo-accélération
    """
    
    # Conversion en tableaux numpy
    freqs = np.asarray(frequencies)
    acc = np.asarray(accelerogram)
    
    # Initialisation des résultats
    n_freqs = len(freqs)
    spectre_d = np.zeros(n_freqs)
    spectre_v = np.zeros(n_freqs)
    spectre_a = np.zeros(n_freqs)
    
    # Pré-calculs constants
    xi = damping
    xi_sq = np.sqrt(1 - xi**2)  # Amortissement réduit (sqrt(1-xi^2))
    
    # Boucle sur chaque fréquence (comme dans la macro VBA "For fn = 1 To Nf_")
    for i, f in enumerate(freqs):
        omega = 2 * np.pi * f
        omega_d = omega * xi_sq
        w2 = omega**2
        w3 = omega**3
        
        # --- Calcul des matrices A et B (Méthode NIGAM) ---
        # Voir doc R5.05.01 Page 5 
        
        exp_val = np.exp(-xi * omega * dt)
        sin_val = np.sin(omega_d * dt)
        cos_val = np.cos(omega_d * dt)
        
        # Coefficients de la matrice A
        a11 = exp_val * (xi / xi_sq * sin_val + cos_val)
        a12 = exp_val / omega_d * sin_val
        a21 = -omega / xi_sq * exp_val * sin_val
        a22 = exp_val * (cos_val - xi / xi_sq * sin_val)
        
        # Matrice A
        A_mat = np.array([[a11, a12],
                          [a21, a22]])
        
        # Coefficients intermédiaires pour la matrice B
        c1 = (2 * xi**2 - 1) / (w2 * dt)
        c2 = (2 * xi) / (w3 * dt)
        
        # Coefficients de la matrice B (partie liée au chargement)
        # Note: termes identiques à la macro VBA
        b11 = exp_val * ((c1 + xi/omega) * sin_val/omega_d + (c2 + 1/w2) * cos_val) - c2
        b12 = -exp_val * (c1 * sin_val/omega_d + c2 * cos_val) - 1/w2 + c2
        
        b21 = exp_val * ((c1 + xi/omega) * (cos_val - xi/xi_sq * sin_val) - 
                         (c2 + 1/w2) * (omega_d * sin_val + xi*omega * cos_val)) + 1/(w2 * dt)
                         
        b22 = -exp_val * (c1 * (cos_val - xi/xi_sq * sin_val) - 
                          c2 * (omega_d * sin_val + xi*omega * cos_val)) - 1/(w2 * dt)
        
        B_mat = np.array([[b11, b12],
                          [b21, b22]])
        
        # --- Résolution temporelle (Boucle sur le temps) ---
        # Initialisation vecteur état Q = [déplacement, vitesse]
        q = np.zeros(2)
        d_max = 0.0
        
        # On parcourt l'accélérogramme pas à pas
        # Correspond à "For tn = 2 To NB_val_T_"
        for t in range(1, len(acc)):
            alpha_vec = np.array([acc[t-1], acc[t]]) # [alpha(t-dt), alpha(t)]
            
            # Relation de récurrence : Q(t) = A * Q(t-dt) + B * Alpha
            q = A_mat @ q + B_mat @ alpha_vec
            
            # Recherche du maximum absolu du déplacement
            if abs(q[0]) > d_max:
                d_max = abs(q[0])
        
        # Stockage des résultats pour cette fréquence
        spectre_d[i] = d_max                  # Déplacement relatif max
        spectre_v[i] = omega * d_max          # Pseudo-vitesse (omega * SD)
        spectre_a[i] = w2 * d_max             # Pseudo-accélération (omega^2 * SD)

    return {
        "Frequences": freqs,
        "Deplacement": spectre_d,
        "Pseudo_Vitesse": spectre_v,
        "Pseudo_Acceleration": spectre_a
    }


def affichage(temps, accel, resultats):
    """Entrées : 
        -   résulats de la méthode précédente
        -   échelle des temps et entrée utilisées
    """
    plt.subplot(2,2,1)
    plt.grid(True, which="both", ls="-")
    plt.plot(temps, accel)
    plt.xlabel("Temps (s)")
    plt.ylabel('Accélération (m/s² ou g)')
    plt.title('Perturbation en accélération')
    
    plt.subplot(2,2,2)
    plt.grid(True, which="both", ls="-")
    plt.loglog(resultats['Frequences'], resultats['Deplacement'],color='purple', label='Amortissement 5%')
    plt.xlabel("Fréquence de l'O.H (Hz)")
    plt.ylabel('Déplacement (m)')
    plt.title('Spectre de déplacement relatif')
    plt.legend(loc='best')
        
    plt.subplot(2,2,3)
    plt.grid(True, which="both", ls="-")
    plt.loglog(resultats['Frequences'], resultats['Pseudo_Vitesse'],color='green', label='Amortissement 5%')
    plt.xlabel("Fréquence de l'O.H (Hz)")
    plt.ylabel('Pseudo-vitesse (m/s)')
    plt.title('Spectre de pseudo-vitesse')
    plt.legend(loc='best')
    
    plt.subplot(2,2,4)
    plt.grid(True, which="both", ls="-")
    plt.loglog(resultats['Frequences'], resultats['Pseudo_Acceleration'],color='coral', label='Amortissement 5%')
    plt.xlabel("Fréquence de l'O.H (Hz)")
    plt.ylabel('Pseudo-accélération (m/s² ou g)')
    plt.title('Spectre de pseudo-accélération')
    plt.legend(loc='best')
    
    plt.tight_layout()
    plt.show()
    
    
# --- Exemple d'utilisation ---

if __name__ == "__main__":

    dt = 0.01
    temps = np.arange(0, 10, dt)
    # Test avec un sinus dégradé
    accel = 1.0 * np.sin(2 * np.pi * 2.0 * temps) * np.exp(-5 * temps) * temps**5 #(m/s² ou g)
    
    #accel =  1.0 * np.sin(2.0 * temps + 3.0)
    
    # Paramètres de calcul
    amortissement = 0.05  # 5%
    frequences = np.logspace(np.log10(0.1), np.log10(50), 100) # De 0.1 à 50 Hz
    
    # Appel de la fonction
    print("Calcul en cours...")
    resultats = sro_nigam(frequences, accel, dt, amortissement)
    
    # Affichage puis export
    print("Calcul terminé.")
    print(f"SRO Max (Accélération): {np.max(resultats['Pseudo_Acceleration']):.4f}")
    
    affichage(temps, accel, resultats)
    

    
# Lire le fichier Excel (feuille Calculs)
df_calculs = pd.read_excel("Spectre.xlsm", sheet_name="Calculs", header=None)
dt = df_calculs.iloc[10, 3] # Cellule D11 (index 0 donc ligne 10, col 3)
amort = df_calculs.iloc[19, 3] # Cellule D20

# Lire l'accélérogramme (feuille A)
df_acc = pd.read_excel("Spectre.xlsm", sheet_name="A")
accel_data = df_acc.iloc[2:, 1].values # Colonne 2
temps = df_acc.iloc[2:, 0].values # Colonne 1, sans f=0

# Lancer le calcul
res = sro_nigam(frequencies, accel_data, dt, amort)

affichage(temps, accel_data, res)

# Sauvegarder
pd.DataFrame(res).to_excel("Spectre.xlsm", sheet_name="RES")
