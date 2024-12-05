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
- **`genOperators.m`**: Generates the constant matrices \(c\) and \(A\) when computing the state vector \(y=[v,w]^T\).
- **`PODModes.m`**: Finds the first \(k\) modes of the given singular values using Relative Information Content (RIC).

---

## Dependencies

This project requires MATLAB.

---

## Usage

1. Clone the repository:

```bash
git clone https://github.com/JonDamFlindt/FitzHugh-Nagumo.git
cd FitzHugh-Nagumo
```

2. Run the simulation:

```matlab
FitzROM
```

3. Customize parameters in `genOperators.m` to explore different settings of the FitzHugh-Nagumo system or POD truncation levels.

---

## Results

The output includes:

- Time evolution of \(v\) and \(w\) for the full and reduced-order systems.
- Comparison plots showing the accuracy of the reduced-order model.
- Singular value spectrum for evaluating the effectiveness of the chosen truncation.

---

## Acknowledgments

This project is part of a bachelor thesis on model order reduction at the University of Southern Denmark (SDU). Portions of the code were provided by my thesis supervisor, professor Ralf Zimmermann. For more information on POD, its derivation and its applications, consult:

- Kutz, J. N., & Brunton, S. L. *Data-Driven Science and Engineering*
- Stewart, G. W. *Matrix Algorithms Volume 1: Basic Decompositions*

