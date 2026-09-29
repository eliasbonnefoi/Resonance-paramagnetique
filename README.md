# TP numérique de Physique Quantique — RPE dans le rubis

**Résonance Paramagnétique Électronique de l'ion Cr³⁺ dans le rubis**

Elias Bonnefoi — ESPCI Paris PSL, 1ère année (2025-2026)

## Contenu

- [`rapport.pdf`](rapport.pdf) — rapport complet du TP (construction de l'hamiltonien de spin, diagonalisation numérique, carte des transitions RPE, couplage spin-cavité Ω(T))
- [`simulation_rpe_rubis.py`](simulation_rpe_rubis.py) — code Python de simulation

## Résumé

Le TP étudie la structure fine des niveaux d'énergie de l'ion Cr³⁺ (spin `S = 3/2`) dans une matrice de rubis, via la diagonalisation numérique de l'hamiltonien de spin combinant champ cristallin et effet Zeeman :

```
H/h = D·(Sz² − S²/3) + (g·μ_B·B0/h)·(cos θ·Sz + sin θ·Sx)
```

Trois volets sont traités :

1. **Niveaux d'énergie** en fonction de `B0` pour un angle donné, et identification des transitions résonantes avec la cavité micro-onde (`f_cav = 9.19 GHz`).
2. **Carte des transitions RPE** dans le plan `(θ, B0)`, comparée aux données expérimentales de 2018.
3. **Couplage spin-cavité Ω(T)** : dépendance en température du couplage collectif spin-photon, pour deux orientations du champ (`θ = 36°` et `θ = 98°`).

## Dépendances

```bash
pip install numpy matplotlib scipy
```

## Exécution

```bash
python simulation_rpe_rubis.py
```
