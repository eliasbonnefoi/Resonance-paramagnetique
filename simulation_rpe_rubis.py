# -*- coding: utf-8 -*-
"""
TP numérique de Physique Quantique 2025-2026
Résonance Paramagnétique Électronique (RPE) dans le rubis
Elias Bonnefoi — ESPCI Paris PSL
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.constants import hbar, k, e, m_e, physical_constants

mu_B = physical_constants['Bohr magneton'][0]
h    = physical_constants['Planck constant'][0]
g_e  = 2.0023

D_freq = 5.73e9
D = D_freq

f_cav = 9.19e9

S_val = 3/2
ms_vals = np.array([3/2, 1/2, -1/2, -3/2])

def construire_operateurs_spin(S=3/2):
    dim = int(2*S + 1)
    ms = np.array([S - i for i in range(dim)])

    Sz = np.diag(ms).astype(complex)

    Sp = np.zeros((dim, dim), dtype=complex)
    Sm = np.zeros((dim, dim), dtype=complex)
    for i in range(dim - 1):
        m = ms[i+1]
        Sp[i, i+1] = np.sqrt(S*(S+1) - m*(m+1))
        Sm[i+1, i] = np.sqrt(S*(S+1) - m*(m-1))

    Sx = (Sp + Sm) / 2.
    Sy = (Sp - Sm) / (2.j)
    S2 = Sx @ Sx + Sy @ Sy + Sz @ Sz
    return Sx, Sy, Sz, S2, Sp, Sm

Sx, Sy, Sz, S2, Sp, Sm = construire_operateurs_spin(S=3/2)


def hamiltonien(B0_G, theta_deg):
    B0 = B0_G * 1e-4
    theta = np.radians(theta_deg)

    H_cf = D * (Sz @ Sz - S2 / 3.)

    prefact = g_e * mu_B * B0 / h
    H_Z = prefact * (np.cos(theta)*Sz + np.sin(theta)*Sx)

    return H_cf + H_Z

def diagonaliser(B0_G, theta_deg):
    H = hamiltonien(B0_G, theta_deg)
    vals, vecs = np.linalg.eigh(H)
    idx = np.argsort(vals)
    return vals[idx], vecs[:, idx]


theta_2a = 40.
B0_range = np.linspace(0, 6000, 600)

energies_2a = np.zeros((len(B0_range), 4))
for i, B0 in enumerate(B0_range):
    vals, _ = diagonaliser(B0, theta_2a)
    energies_2a[i, :] = vals / 1e9

fig, ax = plt.subplots(figsize=(8, 5))
couleurs = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
labels   = [r'$E_1$', r'$E_2$', r'$E_3$', r'$E_4$']
for j in range(4):
    ax.plot(B0_range, energies_2a[:, j], color=couleurs[j], label=labels[j])
ax.axhline(f_cav/1e9, color='k', ls='--', lw=1.2,
           label=fr'$f_{{\rm cav}}={f_cav/1e9}$ GHz')
ax.set_xlabel(r'$B_0$ (G)', fontsize=13)
ax.set_ylabel(r'Énergie / $h$ (GHz)', fontsize=13)
ax.set_title(fr"Niveaux d'énergie du Cr$^{{3+}}$ dans le rubis — $\theta = {theta_2a}°$",
             fontsize=13)
ax.legend(fontsize=11)
ax.set_xlim(0, 6000)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.show()

print("\n=== Vecteurs propres à B0=0 G (theta=40°) ===")
vals0, vecs0 = diagonaliser(0, theta_2a)
base = ['|3/2>', '|1/2>', '|-1/2>', '|-3/2>']
for j in range(4):
    c = vecs0[:, j]
    expr = " + ".join([f"{c[k]:.3f}·{base[k]}" for k in range(4) if abs(c[k]) > 0.01])
    print(f"  |E{j+1}> = {expr}   (E={vals0[j]/1e9:.3f} GHz)")

print("\n=== Vecteurs propres à B0=6000 G (theta=40°) ===")
vals6, vecs6 = diagonaliser(6000, theta_2a)
for j in range(4):
    c = vecs6[:, j]
    expr = " + ".join([f"{c[k]:.3f}·{base[k]}" for k in range(4) if abs(c[k]) > 0.01])
    print(f"  |E{j+1}> = {expr}   (E={vals6[j]/1e9:.3f} GHz)")

print(f"\n=== Transitions à f_cav={f_cav/1e9} GHz pour theta={theta_2a}° ===")
n_trans = 0
for i in range(4):
    for j in range(i+1, 4):
        for B0 in B0_range:
            vals, _ = diagonaliser(B0, theta_2a)
            delta = abs(vals[j] - vals[i]) / 1e9
            if abs(delta - f_cav/1e9) < 0.05:
                n_trans += 1
                break
print(f"  -> {n_trans} transitions visibles")


theta_range = np.linspace(0, 190, 191)
B0_fin = np.linspace(0, 6000, 1200)
pts_B0    = []
pts_theta = []

for th in theta_range:
    for B0 in B0_fin:
        vals, _ = diagonaliser(B0, th)
        for i in range(4):
            for j in range(i+1, 4):
                ecart = (vals[j] - vals[i]) / 1e9
                if abs(ecart - f_cav/1e9) < 0.04:
                    pts_theta.append(th)
                    pts_B0.append(B0)

fig2, ax2 = plt.subplots(figsize=(9, 6))
ax2.scatter(pts_theta, pts_B0, s=1, color='tab:blue', alpha=0.6)

donnees_exp = {
    0:  [750, 1050, 3700, 3300],
    20: [1050, 2620, 4220, 5750],
    40: [460, 1500, 2600, 4400],
    60: [550, 2040, 3600],
    80: [900, 2140, 4820],
}
for th_exp, liste_B0 in donnees_exp.items():
    for B0_exp in liste_B0:
        ax2.scatter(th_exp, B0_exp, marker='x', s=80, color='red', zorder=5)

ax2.set_xlabel(r'Angle $\theta$ (°)', fontsize=13)
ax2.set_ylabel(r'$B_0$ (G)', fontsize=13)
ax2.set_title(r'Carte des transitions RPE — $f_{\rm cav}=9.19$ GHz', fontsize=13)
ax2.set_xlim(0, 190)
ax2.set_ylim(0, 6100)
from matplotlib.lines import Line2D
legende = [Line2D([0],[0], marker='o', color='w', markerfacecolor='tab:blue',
                  markersize=6, label='Simulation'),
           Line2D([0],[0], marker='x', color='red', markersize=8, lw=0,
                  label='Données expérimentales (2018)')]
ax2.legend(handles=legende, fontsize=11)
ax2.grid(alpha=0.3)
plt.tight_layout()
plt.show()


rho_rubis   = 3.98e3
V_ech       = 1e-9
M_Al2O3     = (2*26.98 + 3*16.) * 1e-3
N_A         = 6.022e23
n_Al2O3     = (rho_rubis / M_Al2O3) * N_A
n_Cr        = 0.001 * 2 * n_Al2O3
N_spins     = int(n_Cr * V_ech)
print(f"\n=== Nombre de spins Cr dans 1 mm³ : N = {N_spins:.3e} ===")

b1 = 100e-12
kB = k


def facteur_polarisation(Ea_Hz, Eb_Hz, niveaux_Hz, T):
    E = np.array(niveaux_Hz) * h
    E0 = np.min(E)
    Z = np.sum(np.exp(-(E - E0) / (kB * T)))
    pa = np.exp(-(Ea_Hz*h - E0) / (kB * T))
    pb = np.exp(-(Eb_Hz*h - E0) / (kB * T))
    return (pa - pb) / Z


def elem_matrice_Sy(vec_a, vec_b):
    return np.dot(np.conj(vec_a), Sy @ vec_b)


def calculer_omega(theta_deg, T_arr):
    B0_scan = np.linspace(0, 6000, 2400)
    paires  = []

    for idx in range(len(B0_scan) - 1):
        B0 = B0_scan[idx]
        vals, vecs = diagonaliser(B0, theta_deg)
        for i in range(4):
            for j in range(i+1, 4):
                ecart = (vals[j] - vals[i]) / 1e9
                if abs(ecart - f_cav/1e9) < 0.025:
                    paires.append((i, j, vals.copy(), vecs.copy()))
                    break

    resultats = []
    for (i, j, vals_res, vecs_res) in paires:
        me = elem_matrice_Sy(vecs_res[:, i], vecs_res[:, j])
        Omega_T = np.zeros(len(T_arr))
        for k_T, T in enumerate(T_arr):
            p = facteur_polarisation(vals_res[i], vals_res[j], vals_res, T)
            Omega_T[k_T] = (g_e * mu_B * b1 / (2*h)) * abs(me) * np.sqrt(N_spins * abs(p))
        resultats.append((i, j, Omega_T))
    return resultats


T_arr = np.logspace(np.log10(0.01), np.log10(30), 300)

fig3, ax3 = plt.subplots(figsize=(8, 5))
res_36 = calculer_omega(36, T_arr)
coul_trans = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
for nr, (i, j, Omega_T) in enumerate(res_36):
    ax3.plot(T_arr, Omega_T/1e6, color=coul_trans[nr % 5],
             label=fr'$|E_{i+1}\rangle \to |E_{j+1}\rangle$')
ax3.set_xscale('log')
ax3.set_xlabel('Température $T$ (K)', fontsize=13)
ax3.set_ylabel(r'$\Omega / 2\pi h$ (MHz)', fontsize=13)
ax3.set_title(r'Couplage spin-cavité $\Omega(T)$ — $\theta = 36°$', fontsize=13)
ax3.legend(fontsize=10)
ax3.grid(alpha=0.3, which='both')
plt.tight_layout()
plt.show()

fig4, ax4 = plt.subplots(figsize=(8, 5))
res_98 = calculer_omega(98, T_arr)
for nr, (i, j, Omega_T) in enumerate(res_98):
    ax4.plot(T_arr, Omega_T/1e6, color=coul_trans[nr % 5],
             label=fr'$|E_{i+1}\rangle \to |E_{j+1}\rangle$')
ax4.set_xscale('log')
ax4.set_xlabel('Température $T$ (K)', fontsize=13)
ax4.set_ylabel(r'$\Omega / 2\pi h$ (MHz)', fontsize=13)
ax4.set_title(r'Couplage spin-cavité $\Omega(T)$ — $\theta = 98°$', fontsize=13)
ax4.legend(fontsize=10)
ax4.grid(alpha=0.3, which='both')
plt.tight_layout()
plt.show()
