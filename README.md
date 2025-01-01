# Proper Orthogonal Decomposition (POD) on the FitzHugh-Nagumo System

This repository contains a MATLAB implementation of Proper Orthogonal Decomposition (POD) applied to the FitzHugh-Nagumo system, a simplified model of excitable systems like neurons and cardiac cells. This ReadMe provides an overview of the project, the theoretical background, the implementation details, and instructions for usage.

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Background: FitzHugh-Nagumo System](#background-fitzhugh-nagumo-system)
3. [Proper Orthogonal Decomposition (POD)](#proper-orthogonal-decomposition-pod)
4. [Implementation](#implementation)
5. [Dependencies](#dependencies)
6. [Usage](#usage)
7. [Results](#results)
8. [Acknowledgments](#acknowledgments)

---

## Project Overview

This project demonstrates the application of POD for model order reduction on the FitzHugh-Nagumo system. By reducing the dimensionality of the system while retaining its essential dynamics, the computational cost of simulations can be significantly reduced.

---

## Background: FitzHugh-Nagumo System

The FitzHugh-Nagumo model is a two-variable simplification of the Hodgkin-Huxley model used to describe electrical impulses in excitable cells. The system of equations is given by:


`\varepsilon \frac{\partial v}{\partial t} = \varepsilon^2 \frac{\partial^2 v}{\partial x^2} - v(v-0.1)(1-v) - w + c, \\`
`\frac{\partial w}{\partial t} = b v - \gamma w + c,`


where:
- `v` represents the membrane potential,
- `w` is the recovery variable,
- `\varepsilon`, `\gamma`, `c`, and `b` are parameters.

This model captures the essential features of excitability and refractory periods.

---

## Proper Orthogonal Decomposition (POD)

POD is a technique for reducing the dimensionality of a system by projecting the dynamics onto a lower-dimensional subspace spanned by modes derived from the singular value decomposition (SVD) of snapshot data. This approach retains the dominant features of the system while discarding less significant components.

**Steps in POD:**

1. Generate a snapshot matrix by simulating the full system over time.
2. Compute the SVD of the snapshot matrix to extract the dominant modes.
3. Project the system onto these modes to obtain a reduced-order model.

---

## Implementation

This repository contains the following files:


- **`bc.m`**: Defines the boundary conditions of the FitzHugh-Nagumo system.
- **`nonlin.m`**: Defines the nonlinear behavior \(f=v(v-0.1)(1-v)\) of the system.
- **`genOperators.m`**: Generates the constant matrices \(c\) and \(A\), as well as defines some other constants, when computing the state vector \(y=[v,w]^T\).
- **`genTime.m`**: Generates the time constants \(tList\) and \(dt\), which is the list of time steps (e.g., \(1, 2, 3, \dots\)) and the size of the time steps \(\Delta t\).
- **`PODModes.m`**: Finds the first \(k\) modes of the given singular values using Relative Information Content (RIC).
- **`fitzFOM.m`**: Solves the full-order model (FOM) of the FitzHugh-Nagumo system, producing high-fidelity snapshots of the solution for both \(v\) and \(w\) over time. The generated snapshots are saved to a `.mat` file.
- **`fitzROM.m`**: Solves the reduced-order model (ROM) derived from Proper Orthogonal Decomposition (POD). Uses the precomputed POD basis to reconstruct approximate solutions for \(v\) and \(w\) with significantly reduced computational cost.
- **`script_fitz.m`**: Main file that runs `fitzROM.m` and plots the full-order solutions alongside the reduced-order solutions. Includes comments on how to set up and use the project effectively.

---

## Dependencies

This project requires MATLAB 2024b or newer. Older versions may work, but is not guaranteed.

---

## Usage

1. Press 'Code' and 'Download ZIP', or clone the repository:

```bash
git clone https://github.com/JonDamFlindt/FitzHugh-Nagumo.git
cd FitzHugh-Nagumo
```

2. Run the simulation:

```matlab
script_fitz
```

3. Customize parameters in `genOperators.m`, `genTime.m` and `script_fitz.m` to explore different settings of the FitzHugh-Nagumo system, and `fitzROM.m` for different error thresholds if using RIC for POD modes.

---

## Results

The output includes:

- Singular value spectrum for evaluating the effectiveness of the chosen truncation.
- Time evolution of \(v\) and \(w\) for the full and reduced-order systems.
- Comparison plots showing the accuracy of the reduced-order model.

---

## Acknowledgments

This project is part of a bachelor thesis on model order reduction at the University of Southern Denmark (SDU). Portions of the code were provided by my thesis supervisor, Professor Ralf Zimmermann. For more information on POD, its derivation and its applications, consult:

- Kutz, J. N., & Brunton, S. L. *Data-Driven Science and Engineering*
- Stewart, G. W. *Matrix Algorithms Volume 1: Basic Decompositions*
