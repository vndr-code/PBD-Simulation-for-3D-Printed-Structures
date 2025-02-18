# PBD Simulation for 3D Printed Structures

This project implements a Position-Based Dynamics (PBD) simulation tailored for 3D printed structures. The simulation models the behavior of a structure (e.g., a printed concrete object) using particles, spring constraints, and collision handling. The goal is to capture key physical phenomena such as layer stiffening (via age-dependent parameters) and collapse detection.

## Overview

The simulation is built on the following core ideas:

- **Particle Representation:**  
  Each point along a pre-generated path (loaded from a CSV) becomes a *Particle* with properties such as position, velocity, mass, radius, and an artificially assigned age (based on its layer).

- **Layer & Age Assignment:**  
  Particles are grouped into layers according to their z-coordinate. An artificial age is assigned based on the layer index (using a configurable age factor) to mimic the curing process—older layers become stiffer and less plastic.

- **Spring Constraints:**  
  Two types of springs connect the particles:
  - **Primary Springs:** Connect consecutive particles along the path.
  - **Secondary Springs:** Connect nearby particles (within a defined neighbor radius) to provide lateral support. For these springs, effective stiffness, plasticity, and yield ratio are adjusted based on particle age.

- **Dynamics and Constraint Projection:**  
  The simulation updates particle positions using a simple Euler integration with gravity and enforces boundary conditions. Collision detection prevents particles from overlapping. Spring constraints are projected (multiple times per timestep) to enforce distance constraints.

- **Collapse Detection:**  
  The simulation monitors total kinetic energy. If energy spikes above a configurable threshold, it is assumed that the structure has collapsed, and the simulation is frozen.

- **Output:**  
  Particle positions are recorded over time. The simulation can export an OBJ mesh and generate animated GIFs to visualize the structure and its progressive collapse.

## Project Structure

- **main.py**  
  Orchestrates the simulation: loads the path, creates particles, generates primary and secondary springs, runs the simulation, and calls visualization routines.

- **simulation.py**  
  Contains the core simulation code:
  - Particle class with properties such as position, velocity, mass, radius, age, and layer.
  - Functions to update particle dynamics (`simulate_particle()`), handle collisions (`handle_particle_collisions()`), assign layers and ages (`assign_layers()`), and create particles and springs.
  - `run_simulation()` loops over time steps, applies dynamics, projects spring constraints, detects collapse, and records positions.

- **path_loader.py**  
  Provides functions to load a pre-generated path from a CSV file and to compute rest lengths between consecutive points.

- **spring.py**  
  Defines the *Spring* class and its `project()` method.  
  For secondary springs, effective parameters (stiffness, plasticity, yield ratio) are adjusted based on the age of the connected particles.

- **visualization.py**  
  Contains functions to visualize simulation output:
  - `save_gif()` animates the particles with their connecting springs.
  - Variants such as `save_gif_path()`, `save_layer_gif_path()`, `save_gif_smooth()`, and `save_gif_smooth_gradient()` provide different visualization options (e.g., using smoothing or color gradients).

- **config.py**  
  Stores global constants used throughout the simulation (e.g., boundary limits, neighbor radius, age factors, and limits for spring parameters).

## Installation

### Dependencies

This project requires:
- Python 3.x
- [NumPy](https://numpy.org/)
- [Matplotlib](https://matplotlib.org/)
- [scikit-image](https://scikit-image.org/) (for marching cubes if using TPMS functionalities)

You can install these dependencies using pip:

```bash
pip install numpy matplotlib scikit-image

