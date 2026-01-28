import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# --- Paramètres Globaux ---
dt = 0.1
t_fin = 5.0
TIME = np.arange(0, t_fin, dt)

# Paramètres de l'oscillateur
m = 20          # kg
xi = 0.05       # amortissement
f_Force = 15    # Hz
wf = 1.0        # pulsation forcée (rad/s)
wn = 2*np.pi*f_Force # pulsation propre
wd = wn*np.sqrt(1-xi**2) # pulsation amortie
F0 = 100        # Amplitude

# Définition de la force (Entrée)
# Cas 1 : Sinus
# force = F0 * np.sin(wf * TIME)
# Cas 2 : Dirac (Impulsion unitaire à t=0 pour tester)
F = np.zeros_like(TIME)
F[0] = 1.0 / dt  # Approximativement 1/dt pour que l'intégrale vaille 1

# --- Fonctions Corrigées ---

def Integrande(F_val, t, tau):
    """
    Fonction sous l'intégrale de Duhamel : h(t-tau) * F(tau)
    h(t) = (1/m*wd) * exp(-xi*wn*t) * sin(wd*t)
    """
    if t < tau:
        return 0.0
    return (1/(m*wd)) * np.exp(-xi*wn*(t - tau)) * F_val * np.sin(wd*(t - tau))

def Trapeze_numpy(t, F_array, time_array):
    """Méthode des trapèzes utilisant numpy.trapz"""
    # On intègre de tau=0 à tau=t
    indices = np.where(time_array <= t)[0]
    if len(indices) < 2:
        return 0.0
    
    # On calcule l'intégrande pour tous les tau concernés
    y_vals = [Integrande(F_array[j], t, time_array[j]) for j in indices]
    
    # Intégration
    return np.trapz(y_vals, x=time_array[indices])

def Trapeze_manuel_optimise(TIME, F):
    """
    Méthode des trapèzes optimisée (O(N)).
    On utilise la décomposition trigonométrique pour éviter de tout recalculer.
    """
    U = np.zeros(len(TIME))
    
    # Pré-calcul vectoriel des termes A et B (éviter les boucles lentes)
    # A(tau) = exp(xi*wn*tau) * F(tau) * cos(wd*tau)
    # B(tau) = exp(xi*wn*tau) * F(tau) * sin(wd*tau)
    # Note: On utilise ici l'exponentielle positive pour compenser celle du terme final
    
    exp_term = np.exp(xi * wn * TIME)
    Y_A = exp_term * F * np.cos(wd * TIME)
    Y_B = exp_term * F * np.sin(wd * TIME)
    
    ACum = 0.0
    BCum = 0.0
    
    # Boucle unique cumulative
    for i in range(1, len(TIME)):
        # Aire Trapèze = 0.5 * dt * (y[i] + y[i-1])
        ACum += 0.5 * dt * (Y_A[i] + Y_A[i-1])
        BCum += 0.5 * dt * (Y_B[i] + Y_B[i-1])
        
        # Reconstruction de la solution à l'instant t
        term_common = (1/(m*wd)) * np.exp(-xi * wn * TIME[i])
        U[i] = term_common * (ACum * np.sin(wd * TIME[i]) - BCum * np.cos(wd * TIME[i]))
        
    return U

def main():
    print("Calcul en cours...")
    
    # 1. Méthode Numpy (Lente car recalculée à chaque pas)
    start = datetime.now()
    # On utilise une liste compréhension pour calculer U à chaque instant t
    U_numpy = [Trapeze_numpy(t, F, TIME) for t in TIME]
    end = datetime.now()
    print(f"Temps exécution Numpy (trapz) : {end-start}")
    
    # 2. Méthode Manuelle Optimisée (Rapide)
    start = datetime.now()
    U_manuel = Trapeze_manuel_optimise(TIME, F)
    end = datetime.now()
    print(f"Temps exécution Manuel Optimisé : {end-start}")

    # Affichage
    plt.figure(figsize=(10, 6))
    plt.plot(TIME, U_numpy, 'o', label="Numpy trapz", markersize=2)
    plt.plot(TIME, U_manuel, label="Manuel Optimisé")
    plt.xlabel('Temps (s)')
    plt.ylabel('Déplacement')
    plt.title("Réponse impulsionnelle (Duhamel)")
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
