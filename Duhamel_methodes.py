#%%
"""
Méthodes de calcul de l'intégrae de Duhamel
pour la production d'un spectre de réponse d'un oscillateur harmonique 

Auteur : Laura Brémont (stagiaire)
Dernière version : 26/01/26
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy as sym # type: ignore
import scipy.integrate as integrate
import scipy.linalg as alg
from datetime import datetime


### Entrée : Rampe
"""

fig = plt.figure()
axes = fig.add_axes([0.1,0.1,2,1])
axes.plot(t,p)

axes.set_xlabel('t (sec)')
axes.set_ylabel('Force (N)')
axes.set_title('Entrée en rampe ')
axes.set_xlim([0,t1])
axes.set_ylim([0,P0])
plt.grid()
plt.show()"""



""""---------------------------SOLUTION ANALYTIQUE-----------------------------"""
sym.init_printing() # Ecrire automatiquement des expressions bien formatées

##Import de symboles depuis la bibliothèque Sympy
m = sym.Symbol('m')
wf = sym.Symbol('wf')
wn = sym.Symbol('wn')
F0 = sym.Symbol('F0')
t1 = sym.Symbol('t1')
tau = sym.Symbol('tau')
t = sym.Symbol('t')

##Intégrande de l'intégrale de Duhamel

#Cas s'une rampe :
"""
f = tau * sym.sin(wf*t-wf*tau)
"""
"""
#Cas d'un sinus (harmonique) :
f = sym.sin(wn*tau)*sym.sin(wf*t-wf*tau)


defint = sym.integrate(f, (tau, 0,t)) #tau est la variable d'intégration

print("SymPy a généré cette expression pour l'intégrale définie:")

sym.simplify(defint) #factorise, réduit, bref simplifie l'expression

"""

""""-----------------------AFFICHAGE (application numérique)------------------------"""
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 14

#On redéfinit les constantes pour l'affichage

dt = 0.1
dtho = dt
DECOUPAGE = np.arange(0, 5, 0.1)
TIME = DECOUPAGE*dt

F0 = 1000       #N
wf = 1          #sec
delT = 0.1      #sec
m = 20          #kg

# Entrée en sinus
f_Force = 15                 #Frequence forcée : Hz
wf = 2*np.pi*f_Force         #Pulsation forcée : rads/s
F0 = 100                       #Amplitude : N
force = F0*np.sin(wf*TIME)

#Entrée en dirac
F = np.concatenate(( np.array([1]), np.array([0]*(len(TIME) - 1)) ))

xi = 0.05            #taux d'amortissement de l'oscillateur
wn = 2*np.pi*wf      #pulsation propre associée à T : rads/s
wd = wn*np.sqrt(1-xi**2)    #pulsation amortie : rads/s
k = m*wn**2          #raideur déduite a partir de T et m connues : kg.m/s**2
    


liste_périodes = [3,4,5] #périodes d'entrées harmoniques à tester 


plt.figure()

"""
for pr in liste_périodes:
    wf2 = pr*wf
    xi = 0.05         #taux d'amortissement de l'oscillateur
    wn = 2*np.pi*wf2    #pulsation propre associée à T : rads/s
    wd = wn*np.sqrt(1-xi**2)    #pulsation amortie : rads/s
    k = m*wn**2       #raideur déduite a partir de T et m connues : kg.m/s**2
    
    
    #Si entrée en rampe :
    
    TIME = np.arange(0, t1+delT, delT) #échelle des temps en abscisses
    u = (F0/k)*(TIME -(np.sin(wn*TIME))) #avec la constante prise en compte
    plt.plot(TIME,u/(F0/k), label=f"T = {pr}t1")
    
    
    #Si entrée en sinus (harmonique) :
    
    force = F0*np.sin(wf2*TIME)   #Force imposée : N
    
    #comme wn>0 et wf>0 :
    u = (F0/m/wd)*(wf2*np.sin(wn*TIME) - wn*np.sin(wn*TIME))/(wf2**2 - wn**2) #avec la constante prise en compte
    plt.plot(TIME, u, label=f"wn = {pr}wf")
    
"""
def main():
    #Comparaison des solutions numériques :
    """ 
    start = datetime.now()
    U = [ Trapeze_numpy(TIME[i], F) for i in range(len(TIME)) ]
    end = datetime.now()
    plt.plot(TIME,U, label=f"Trapèzes numpy : tmp exe = {end-start}sec")
    """
    
    start = datetime.now()
    U = Trapeze_manuel(TIME, F)
    end = datetime.now()
    plt.plot(TIME,U, label=f"Trapèzes manuels : tmp exe = {end-start}sec")

    start = datetime.now()
    U = [ Rectangle_numpy(TIME[i], F) for i in range(len(TIME)) ]
    end = datetime.now()
    plt.plot(TIME, U, label=f"Rectangles numpy : tmp exe = {end-start}sec")

    plt.xlabel('t en sec')         
    plt.ylabel('Déplacement')
    plt.title("Réponse d'un oscillateur harmonique à F")
    plt.legend(loc='best')
    plt.xlim([0, TIME[-1]])
    plt.grid()
    plt.show()

""""---------------------------SOLUTIONS NUMERIQUES-----------------------------"""

def Integrande(F, t, j):
    return (1/m/wd)*np.e*(xi*wn*t)*F[j]*np.sin(wd*(t - j*dtho))

def Trapeze_numpy(t, F):
    """Méthode des trapèzes avec la bibliothèque Numpy"""
    
    return np.trapezoid(lambda j : Integrande(F, t, int(j)), x=np.arange(1, int(t/dt), 1) )

    
def Quad_scipy(t, F):
    """Méthode d'intégration de la bibliothèque Scipy"""
    
    return integrate.quad(Integrande, dx=1, args=(F, t))
    
    
def Rectangle_numpy(t, F):    
    """Méthode des rectangles avec la bibliothèque Numpy ; 
    t est un multiple de dt innférieur à len(TIME)*dt = TIME[-1]"""
    
    THO = np.concatenate(( np.array([dtho]*(int(t/dt))), np.array([0]*(len(TIME) - int(t/dt))) ))
    
    INTEGRANDE = np.array([Integrande(F, t, j) for j in range(len(TIME))])
    print(INTEGRANDE)
    return np.dot(INTEGRANDE, THO)
    
    
def Trapeze_manuel(TIME, F):
    """Méthode des trapèzes pour l'intégration discrète"""
    #Initialise les déplacements à 0
    U = np.zeros(len(TIME))
    
    #Dictionnaires pour éviter de calculer deux fois chaque valeur de l'intégrande
    Y_A = {0: np.e**(xi*wn*TIME[0]) * F[0] * np.cos(wd*TIME[0])}
    
    Y_B = {0: np.e**(xi*wn*TIME[0]) * F[0] * np.sin(wd*TIME[0])}
    
    #Initialise pour la somme cumulative à chq étape "+dt"
    ACum_i=0
    BCum_i=0

    for i in range(1, len(TIME)):
            
        if i>0:

            #Calcul de A à t = i*dt
            
            Y_A[i] = np.e**(xi*wn*TIME[i]) * F[i] * np.cos(wd*TIME[i])
        
            ACum_i += 0.5*delT*(Y_A[i]+Y_A[i-1])    #Aire cumulée entre 0 et t
            A_i = (1/(m*wd))*ACum_i                 #Valeur de A à t

            #Calcul de B à t = i*dt
            Y_B[i] = np.e**(xi*wn*TIME[i]) * F[i] * np.sin(wd*TIME[i])
            
            BCum_i += 0.5*delT*(Y_B[i]+Y_B[i-1])    #Aire cumulée entre 0 et t
            B_i = (1/(m*wd))*BCum_i                 #Valeur de B à t

            #Calcule le déplacement à t = i*dt
            U[i] = A_i*np.e**(-xi*wn*TIME[i])*np.sin(wd*TIME[i]) - B_i * np.e**(-xi*wn*TIME[i])*np.cos(wd*TIME[i])

    return U



main()
