# simulation.py
import numpy as np
#from config import INITIAL_LAYERS, BOUNDARIES
import config
from path_loader import load_path
from spring import Spring

class Particle:
    def __init__(self, pos, vel, mass=1.0, radius=1.0):
        self.pos = np.array(pos, dtype=float)
        self.vel = np.array(vel, dtype=float)
        self.mass = mass
        self.radius = radius
        self.inv_mass = 1.0 / mass if mass > 0 else 0.0
        self.age = 0

def simulate_particle(p, dt, boundaries=config.BOUNDARIES, damping=0.1):
    gravity = np.array([0, 0, -9.81])
    p.vel += gravity * dt
    p.pos += p.vel * dt

    # Check x boundaries
    if p.pos[0] < boundaries['x'][0]:
        p.pos[0] = boundaries['x'][0]
        p.vel[0] = -p.vel[0] * damping
    elif p.pos[0] > boundaries['x'][1]:
        p.pos[0] = boundaries['x'][1]
        p.vel[0] = -p.vel[0] * damping

    # Check y boundaries
    if p.pos[1] < boundaries['y'][0]:
        p.pos[1] = boundaries['y'][0]
        p.vel[1] = -p.vel[1] * damping
    elif p.pos[1] > boundaries['y'][1]:
        p.pos[1] = boundaries['y'][1]
        p.vel[1] = -p.vel[1] * damping

    # Check z boundaries
    if p.pos[2] < boundaries['z'][0]:
        p.pos[2] = boundaries['z'][0]
        p.vel[2] = -p.vel[2] * damping
    elif p.pos[2] > boundaries['z'][1]:
        p.pos[2] = boundaries['z'][1]
        p.vel[2] = -p.vel[2] * damping

def handle_particle_collisions(particles, damping=0.8):
    n = len(particles)
    for i in range(n):
        for j in range(i+1, n):
            p1, p2 = particles[i], particles[j]
            diff = p2.pos - p1.pos
            dist = np.linalg.norm(diff)
            min_dist = p1.radius + p2.radius
            if dist < min_dist and dist > 1e-6:
                penetration = min_dist - dist
                normal = diff / dist
                inv_mass_sum = p1.inv_mass + p2.inv_mass
                if inv_mass_sum == 0:
                    continue
                correction = normal * (penetration / inv_mass_sum)
                p1.pos -= p1.inv_mass * correction
                p2.pos += p2.inv_mass * correction
                v_rel = np.dot(p2.vel - p1.vel, normal)
                if v_rel < 0:
                    impulse = v_rel * damping
                    p1.vel += impulse * normal * p1.inv_mass
                    p2.vel -= impulse * normal * p2.inv_mass

def assign_layers(particles, threshold=0.1):
    """
    Assign a layer index to each particle based on its z-coordinate.
    Particles whose z-values differ by less than `threshold` are in the same layer.
    This function sets a new attribute `layer` on each Particle and updates config.INITIAL_LAYERS.
    """
    # Create list of (index, z) tuples.
    indexed_z = [(i, p.pos[2]) for i, p in enumerate(particles)]
    # Sort by z value.
    indexed_z.sort(key=lambda x: x[1])
    
    layer = 0
    layers = {}  # mapping from particle index to layer
    if indexed_z:
        layers[indexed_z[0][0]] = layer
        last_z = indexed_z[0][1]
        for idx, z in indexed_z[1:]:
            if abs(z - last_z) > threshold:
                layer += 1
                last_z = z
            layers[idx] = layer
    
    # Assign the computed layer to each particle.
    for i, p in enumerate(particles):
        p.layer = layers[i]
    
    # Update the global variable.
    #from config import INITIAL_LAYERS  # re-import in case
    config.INITIAL_LAYERS = [layers[i] for i in range(len(particles))]
    return particles

def create_particles_from_path(filename, threshold=0.1):
    points = load_path(filename)
    particles = [Particle(pos, [0, 0, 0], mass=1.0, radius=1.0) for pos in points]
    particles = assign_layers(particles, threshold)

    # this is only when the whole shape is made at once, so the age is artificially implemented
    for p in particles:
        p.age = p.layer * config.AGE_FACTOR

    return particles


# def create_particles_from_path(filename):
#     points = load_path(filename)
#     particles = [Particle(pos, [0, 0, 0], mass=1.0, radius=1.0) for pos in points]
#     return particles

def create_springs(particles, stiffness=0.9, plasticity=0.1, yield_ratio=0.8):
    """
    Create springs connecting consecutive particles.
    """
    springs = []
    for i in range(len(particles) - 1):
        p1 = particles[i].pos
        p2 = particles[i+1].pos
        rest_length = np.linalg.norm(p2 - p1)
        springs.append(Spring(i, i+1, rest_length, stiffness, plasticity, yield_ratio))
    return springs

def create_secondary_springs(particles, neighbor_radius, stiffness=0.1, plasticity=0.9, yield_ratio=0.8):
    """
    Create secondary springs connecting particles that are within a given neighbor_radius.
    """
    """
    Create secondary springs connecting particles that are within a given neighbor_radius.
    
    Parameters:
      particles: List of Particle objects.
      neighbor_radius: Maximum distance to consider for connecting two particles.
      stiffness, plasticity, yield_ratio: Spring parameters.
      
    Returns:
      A list of Spring objects connecting non-consecutive neighbors.
    """
    springs = []
    n = len(particles)
    for i in range(n):
        for j in range(i+2, n):  # Skip immediate consecutive particles.
            pos_i = particles[i].pos
            pos_j = particles[j].pos
            distance = np.linalg.norm(pos_j - pos_i)
            if distance <= neighbor_radius:
                # Mark these as secondary springs.
                springs.append(Spring(i, j, distance, stiffness, plasticity, yield_ratio, secondary=True))
    return springs


def fix_bottom_particles(particles, threshold=0.01):
    """
    Marks particles as fixed if their z coordinate is at or below a threshold.
    """
    for p in particles:
        if p.pos[2] <= threshold:
            p.inv_mass = 0  # Fix the particle by making it immovable.
            p.vel = np.zeros(3)
                                                                                                                        # 12 points about 50000 is ok, high res 36 points 100000
                                                                                                                        # also neighbours (radus) need to change, need to make it different
def run_simulation(num_steps=500, dt=0.02, path_filename='cylinder_path.csv', springs=None, energy_threshold=config.MAX_ENERGY): # energy threshold different for high res and low res
    particles = create_particles_from_path(path_filename)
    
    
    # Fix bottom particles (those at z=0, or nearly 0)
    fix_bottom_particles(particles, threshold=0.01)

    if springs is None:
        # Create primary springs automatically from particles
        springs = create_springs(particles)  
        # Create secondary springs if needed, e.g.:
        # springs += create_secondary_springs(particles, neighbor_radius=5.0)

    positions_over_time = []
    
    collapse_detected = False
    for step in range(num_steps):
        # Update particle positions
        for p in particles:
            simulate_particle(p, dt)
        
        handle_particle_collisions(particles)
        
        # Project spring constraints multiple times
        for _ in range(5):
            for s in springs:
                s.project(particles)
        
        # Record positions
        positions_over_time.append(np.array([p.pos.copy() for p in particles]))
        
        # Compute total kinetic energy as a spike indicator
        total_energy = sum(0.5 * p.mass * np.linalg.norm(p.vel)**2 for p in particles)
        if total_energy > energy_threshold:
            print(f"Collapse detected at step {step}: kinetic energy = {total_energy:.2f}")
            collapse_detected = True
            break  # Freeze simulation once collapse is detected
    
    if collapse_detected:
        # Optionally, fill the remaining time steps with the last recorded positions
        last_positions = positions_over_time[-1]
        for _ in range(num_steps - len(positions_over_time)):
            positions_over_time.append(last_positions)
    
    return np.array(positions_over_time)


