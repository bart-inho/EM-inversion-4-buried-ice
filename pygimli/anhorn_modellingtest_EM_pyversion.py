# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 13:23:56 2021

@author: bart
"""

import matplotlib.pyplot as plt
import numpy as np
import pygimli as pg

# Modèle de resistivité :
    
res_map = np.empty((30,200),dtype='float')
res_map[0:10,:] = 200
res_map[10:20,:] = 20
res_map[20:30,:] = 100

position_x = 98 # position initiale du modèle
noiseEM = 1 
nlay = 40
thickness = 1

thk = pg.Vector(nlay - 1, thickness)  # épaisseur de chaque couche (1m)
centroids = np.cumsum(thk)-thk/2 # centroïde des couches

# On sample dans le model "res_map" la résistivité correspondant aux position de couches. 
res_EM = np.ones(nlay) *  res_map[0,position_x]
    
# Inversion des l'indexing
for i in centroids:
    if i < 30:
        ab = 29 - np.round((i-0.01)).astype(int) # index dans le modèle resmap
        ab2 = np.floor(i/thickness).astype(int) # index dans le modèle EM
        res_EM[ab2] = res_map[ab,position_x]
        
# Mise en commun des modèles
        
model = pg.cat(thk, res_EM)
noiseEM = 1

# Espacement des boucles
coilspacing = 10

nf = 10
freq = pg.Vector(nf)
freq[0] = 100
freq[1] = 400
freq[2] = 600
freq[3] = 1000
freq[4] = 5000
freq[5] = 8000
freq[6] = 9000
freq[7] = 10000
freq[8] = 20000
freq[9] = 30000


# Donction du modèle
fEM = pg.core.FDEM1dModelling(nlay, freq, coilspacing) # coilspacing : espacement des boucles
dataEM = fEM(model)

# Bruit possible
for i in range(len(dataEM)):
    dataEM[i] += np.random.random(1)[0] * 0.1
    
    
# Plot du modèle    
plt.scatter(dataEM[:nf], freq, marker = 'x')
plt.scatter(dataEM[nf:nf*2], freq, marker = 'o')
plt.yscale('log')
plt.ylim((min(freq), max(freq)))
plt.grid(which='both')
plt.xlabel('Phase %')
plt.ylabel('$f$ $[Hz]$');

#Une autre méthode pour comparer. Ca a l'air de donner les memes datas
Data_box_FEM = pg.physics.em.FDEM(x=[1],freqs=freq,coilSpacing=coilspacing)
FOP = pg.core.FDEM1dRhoModelling(centroids, Data_box_FEM.freq(), Data_box_FEM.coilSpacing,-Data_box_FEM.height)
dataEM = FOP(res_EM)

plt.scatter(dataEM[:nf], freq, marker = 'x')
plt.scatter(dataEM[nf:nf*2], freq, marker = 'o')
plt.yscale('log')
plt.ylim((min(freq), max(freq)))
plt.grid(which='both')
plt.xlabel('Phase %')
plt.ylabel('$f$ $[Hz]$');

###############
#INVERSION 
###############

transRhoa = pg.trans.TransLog() # paramètre initial - no limit sur les couches
transThk = pg.trans.TransLog()
transRes = pg.trans.TransLogLU(1., 1000.) # limitation des résistivités - entre 1 et 1000 ohmm
transEM = pg.trans.Trans()

# Nombre de couches
nLayer_inv = 20 # nombre de layer du modèle inverse

# Modèle initial, avant inversion, 2m de couches uniformes à 10 Ohmm
thk_inv = pg.Vector(nLayer_inv - 1, 2)  # 2m d'épaisseur par couches
res_inv = np.ones(nLayer_inv) *  10
initial_model = pg.cat(thk_inv, res_inv)  # mise en commun du modèle


# Paramètres d'inversion -> non uniqueness des résultats

fEM_inv = pg.core.FDEM1dModelling(nLayer_inv, freq, coilspacing)
fEM_inv.region(0).setTransModel(transThk)
fEM_inv.region(1).setTransModel(transRes)
lamEM =  100 # paramètre de régulation pour le modèle
invEM = pg.core.Inversion(dataEM, fEM_inv, transEM)
invEM.setModel(initial_model)
invEM.setRelativeError(0.1/100)
invEM.setLambda(lamEM)
invEM.setMarquardtScheme(0.8)
invEM.setDeltaPhiAbortPercent(0.00001)
invEM.setMaxIter(40)
invEM.setBlockyModel(True)

# Finalisation
modelEM = invEM.run()
respEM = invEM.response()

# Plot de l'inversion

plt.semilogy(dataEM[0:nf], freq, 'bx', label='syn In-Phase')
plt.semilogy(dataEM[nf:nf*2], freq, 'bo', label='syn Out Phase')
plt.semilogy(respEM[0:nf], freq, 'rx', label='syn In-Phase - Inv')
plt.semilogy(respEM[nf:nf*2], freq, 'ro', label='syn Out Phase - Inv')
plt.ylim((min(freq), max(freq)))
plt.xlabel("Phase in %")
plt.ylabel("$f$ in Hz")
plt.grid(which='both')
plt.legend(loc="best")

# Retour au modèle initial
model_Final = np.array(modelEM)
cetr = np.cumsum(model_Final[0:nLayer_inv-2])
res_Final = model_Final[nLayer_inv-1:-2]

# plot du modèle initial
fig, ax = plt.subplots(figsize = ((12,8)))
pg.viewer.mpl.drawModel1D(ax, values=res_Final, depths=cetr,plot='semilogx', 
                          color='blue',label='Initial Model')
pg.viewer.mpl.drawModel1D(ax, values=res_EM[:-1], depths=centroids,plot='semilogx', 
                          color='red',label='Inverted Model')
plt.legend();







