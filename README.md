# Decentralised Molecular Dynamics

A molecular dynamics simulation framework implemented in Rust that models various potential energy functions for diatomic molecules.

## Overview

This project provides a framework for simulating the dynamics of diatomic molecules using different potential energy models:

- **Harmonic Oscillator**: The classical spring model of molecular vibration
- **Morse Potential**: A more realistic anharmonic bond potential that accounts for dissociation
- **Lennard-Jones Potential**: A model for van der Waals interactions between atoms

The simulation uses atomic units (au) throughout its calculations.

## Features

- Simulation of H₂, Hg₂, and Ar₂ molecular dynamics
- Multiple potential energy models with accurate physical constants
- Verlet integration for numerical simulation
- Analytical solutions for harmonic oscillator to compare with numerical methods
- Comprehensive visualization of simulation results:
  - Displacement, force, acceleration, and velocity over time
  - Kinetic, potential, and total energy profiles
  - Comparison plots between different models
- Automatic data export to CSV files

## Getting Started

### Prerequisites

- Rust and Cargo (latest stable version recommended)
- A working Rust development environment

### Installation

1. Clone the repository:
   ```
   git clone https://github.com/yourusername/decentralised-dynamics.git
   cd decentralised-dynamics
   ```

2. Build the project:
   ```
   cargo build --release
   ```

3. Run the simulation:
   ```
   cargo run --release
   ```

## Simulation Parameters

You can modify various simulation parameters in `src/main.rs`:

- `DT_AU`: Time step in atomic units
- `EXP_LEN_AU`: Total simulation length
- `N_STEPS`: Number of simulation steps
- `PRINT_FREQ`: Frequency of console output
- `TEMP_K`: Temperature in Kelvin

Physical constants for each molecule (H₂, Hg₂, Ar₂) are also defined in the same file.

## Output

Simulation results are stored in the `output/` directory, organized by experiment name. Each experiment produces:

- CSV files with full simulation data
- PNG plots showing:
  - Displacement over time
  - Force over time
  - Acceleration over time
  - Velocity over time
  - Total energy over time
  - Combined energy plots (kinetic, potential, and total)

## Physics Background

This simulation implements three important potential energy models:

1. **Harmonic Oscillator Potential**: 
   V(r) = ½k(r-r₀)²
   
   A simple model that works well for small displacements from equilibrium.

2. **Morse Potential**:
   V(r) = D₀[1-e^(-α(r-r₀))]²
   
   Accounts for anharmonicity of real bonds and the possibility of dissociation.

3. **Lennard-Jones Potential**:
   V(r) = 4ε[(σ/r)¹² - (σ/r)⁶]
   
   Models van der Waals interactions with a balance of attractive and repulsive forces.

