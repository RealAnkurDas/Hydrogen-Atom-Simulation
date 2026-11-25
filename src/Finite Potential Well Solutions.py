import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Constants
m = 9.109e-31  # Electron mass (kg)
hbar = 1.055e-34  # Reduced Planck's constant (J·s)
V0 = 0.3 * 1.602e-19  # Convert eV to Joules
L = 10e-9 / 2  # Well width 10 nm, so L = 5 nm

# Define the total energy constant C
C = (2 * m * V0 * L**2) / hbar**2

# Define transcendental equations
def eq1(eta):  # eta tan(eta) = sqrt(C - eta^2)
    return eta * np.tan(eta) - np.sqrt(C - eta**2)

def eq2(eta):  # eta cot(eta) = -sqrt(C - eta^2)
    return -eta / np.tan(eta) - np.sqrt(C - eta**2)

# Find solutions for eta using fsolve
initial_guesses = np.linspace(0.1, np.sqrt(C) - 0.1, 10)  # Increase number of guesses
etas_solutions = []

for guess in initial_guesses:
    try:
        root1 = fsolve(eq1, guess)[0]
        root2 = fsolve(eq2, guess)[0]
        
        if root1 > 0 and root1 not in etas_solutions and np.isfinite(root1):
            etas_solutions.append(root1)
        if root2 > 0 and root2 not in etas_solutions and np.isfinite(root2):
            etas_solutions.append(root2)
    except:
        continue  # Skip problematic values

# Compute unique energy levels
energy_levels = []
for eta in sorted(set(etas_solutions)):  # Remove duplicates
    k = 1.466 / L
    E = (hbar**2 * k**2) / (2 * m) / 1.602e-19  # Convert to eV
    energy_levels.append(E)

# Sort and print energy levels
energy_levels.sort()
print("Energy levels (sorted, in eV):")
for i, E in enumerate(energy_levels):
    print(f"n = {i+1}: E = {E:.4f} eV")

# Plot results
etas = np.linspace(0.01, np.sqrt(C)-0.01, 1000)
xis1 = etas * np.tan(etas)
xis2 = -etas / np.tan(etas)

plt.figure(figsize=(8, 6))
plt.plot(etas, xis1, label=r"$\eta \tan \eta$")
plt.plot(etas, xis2, label=r"$\eta \cot \eta$")
plt.xlabel(r'$\eta$')
plt.ylabel(r'$\xi$')
plt.title('Graphical Solution for Finite Potential Well')
plt.legend()
plt.grid()
plt.show()

t1 = (hbar**3) * (L**2) * m * (V0**3)
t2 = (hbar**2) * (L**4) * (m**2) * (V0**4)
num = (hbar*hbar*V0) - 2*(t1-t2)**0.5
den = 2*hbar*hbar
Jtoev = (num/den)/1.602e-19
print(Jtoev)
