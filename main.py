# Task 1 - Acoustic-gravitational wave propagation in the solar interior
# ULL - 19/10/2022
# Sergio Guerra Arencibia

import matplotlib.pyplot as plt
import numpy as np

# Definition of constants
EPS = 0.000025
RADIUS_SUN = 696340  # Solar radius in km 
GAMMA = 5/3          # Adiabatic index
DELTAZ = 23.3        # Model provided has a delta z of 23.3
DELTAX = 23.3        # We set the delta x equal to delta z
FRECUENCIES = [0.002, 0.003, 0.0035, 0.005]   # Frecuencies in Hz
ZLTP = 500           # Z value where we calculate the LTP

# Function for plotting 2 magnitudes
def plotMagnitudes(xValues, yValues, points = []):
    fig, ax = plt.subplots()
    fig.set_figheight(10)
    fig.set_figwidth(15)
    ax.plot(xValues, yValues)
    if (len(points) != 0):
        for i in range(len(points)):
            ax.scatter(points[i], FRECUENCIES[i], s = 20)
    plt.show()

# Function to change from Z values provided in the solar model to Radius values
def modelZtoRadius(zValues):
    return (RADIUS_SUN + zValues) / RADIUS_SUN

# RK4 generic method for a system of 4 equations
def rk4(function, initialValues, t0, tf, step):
    times = np.linspace(t0, tf, step)
    deltaT = times[1] - times[0]
    values = np.zeros(shape=(4, n))
    values[0, 0] = initialValues[0]
    values[1, 0] = initialValues[1]
    values[2, 0] = initialValues[2]
    values[3, 0] = initialValues[3]

    for i in range(1, n):
        k1 = function(t[i-1], *values[:, i-1])
        k2 = function(t[i-1] + 0.5*deltaT, *values[:, i-1] + 0.5*k1*deltaT)
        k3 = function(t[i-1] + 0.5*deltaT, *values[:, i-1] + 0.5*k2*deltaT)
        k4 = function(t[i-1] + deltaT, *values[:, i-1] + k3*deltaT)
        values[:, i] = values[:, i-1] + (1/6)*deltaT*(k1 + 2*k2 + 2*k3 + k4)
    return (values)

# -------------------------------------------------------------------------
# Reading the model data
zArray = np.array([])
pArray = np.array([])
rhoArray = np.array([])
tArray = np.array([])

with open("model_jcd.dat") as openFileObject:
    next(openFileObject)
    for line in openFileObject:
        currentValues = line.split()
        zArray = np.append(zArray, float(currentValues[0]))
        pArray = np.append(pArray, float(currentValues[1]))
        rhoArray = np.append(rhoArray, float(currentValues[2]))
        tArray = np.append(tArray, float(currentValues[3]))

# # Plotting preassure, density and temperature vs Z
# fig, ax = plt.subplots(3, 1)
# fig.set_figheight(10)
# fig.set_figwidth(15)
# ax[0].plot(zArray, tArray)
# ax[1].plot(zArray, pArray)
# ax[2].plot(zArray, rhoArray)
# plt.show()

# -------------------------------------------------------------------------
# Sound velocity
soundVelocity = np.array([])

def soundVelocityFormula(preassure, density):
    cs = (GAMMA * preassure) / density
    return(np.sqrt(cs))
soundVelocity = soundVelocityFormula(pArray, rhoArray)
# plotMagnitudes(modelZtoRadius(zArray), soundVelocity)

# Height scale
sunGravity = 274 # m/s**2
sunGravity = sunGravity * 100 # to c.g.s

def heightScaleFormula(preassure, density):
    h = preassure / (density * sunGravity)
    return h
scaleHeight = heightScaleFormula(pArray, rhoArray)
# plotMagnitudes(modelZtoRadius(zArray), scaleHeight)

# w_c
def wcFormula(cs, h):
    return (cs / (2 * h))
w_c = wcFormula(soundVelocity, scaleHeight)
# plotMagnitudes(modelZtoRadius(zArray), w_c)

# N
def bruntVaisala (h):
    firstMember = sunGravity / h
    return np.sqrt(firstMember * ((GAMMA - 1) / GAMMA))
n = bruntVaisala(scaleHeight)
# plotMagnitudes(modelZtoRadius(zArray), n)

# ---------------------------------------------------------

# Get the external Z at which the wave reflects
externalZ = np.array([])
frecuencyw_c = w_c / (2*np.pi)
for i in range(len(FRECUENCIES)):
    for j in range(len(w_c)):
        if (np.abs(FRECUENCIES[i] - frecuencyw_c[j]) < EPS):
            externalZ = np.append(externalZ, zArray[j])
            break
# plotMagnitudes(modelZtoRadius(zArray), frecuencyw_c, modelZtoRadius(externalZ))

# Get the k_x of the wave (it's constant since in this axis the medium does not vary)
# We get it from the LTP expression, where we know that k_z = 0
k_xArray = np.array([])  # Units will be cm^-1
for i in range(len(FRECUENCIES)):
    k_xArray = np.append(k_xArray, (2*np.pi*FRECUENCIES[i]) / soundVelocity[ZLTP])
print(k_xArray)