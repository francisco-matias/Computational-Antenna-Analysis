# Computational-Antenna-Analysis

Numerical solution of Hallén’s integral equation for the analysis of a folded dipole antenna using Method of Moments (MoM) discretization.

---

## Overview

This project implements a computational electromagnetic model of a folded dipole antenna based on the integral-equation formulation of Maxwell’s equations.

The current distributions on both arms of the antenna are obtained by solving the coupled Hallén integral equations using a pulse-basis (staircase) discretization approach.

---

## Mathematical Formulation

- Vector potential formulation of radiated fields
- Coupled integral equation system for current distributions
- Even-symmetry constraints
- Additional boundary conditions at the junction
- Matrix formulation of the discretized system

The resulting linear system is solved numerically to extract the unknown currents and constants.

---

## Numerical Method

- Method of Moments (MoM)
- Pulse basis expansion
- Quadrature-based kernel integration
- Dense matrix assembly
- Linear system solution for 2Ns + 2 unknowns

---

## Analysis Performed

- Current distribution convergence analysis (Ns variation)
- Input impedance evaluation
- Resonance identification
- Radiation pattern computation (polar plots)
- Half-Power Beamwidth (HPBW) extraction
- Directivity calculation
- Transmission-line analogy interpretation

---

## Technical Scope

- Computational Electromagnetics
- Integral Equation Methods
- Method of Moments
- Numerical Linear Algebra
- Antenna Theory
- Radiation Pattern Analysis

## Report

Full technical documentation available here:

[Antenna-Analysis.pdf](https://raw.githubusercontent.com/francisco-matias/Computational-Antenna-Analysis/cf2b5ce3a63e78e16cccfc164e8d69dee9083f2b/Antenna%20Analysis.pdf)
