# main.py
#from config import INITIAL_LAYERS, BOUNDARIES
import config
from simulation import run_simulation, create_springs, create_particles_from_path, create_secondary_springs
from visualization import save_gif, save_gif_path, save_layer_gif_path, bezier_gif_path, save_gif_smooth, save_gif_smooth_gradient

def main():
    filename = 'convex_pyramid_path.csv'
    # Create particles from the path.
    particles = create_particles_from_path(filename)
    # Create primary (consecutive) springs:
    primary_springs = create_springs(particles, stiffness=0.9, plasticity=0.1, yield_ratio=0.8)

    # Create secondary springs:
    secondary_springs = create_secondary_springs(particles, neighbor_radius=config.NEIGHBOUR_RADIUS, 
                                                 stiffness=config.STIFFNESS_MIN, plasticity=config.PLASTICITY_MAX, yield_ratio=config.YIELD_RATIO_MIN)

    # Combine them if desired:
    all_springs = primary_springs + secondary_springs

    # Run the simulation (assume run_simulation uses these particles).
    positions = run_simulation(num_steps=100, dt=0.07, path_filename=filename, springs=all_springs )
    # Save the GIF, passing the springs list.
    #print(config.INITIAL_LAYERS)
    #save_layer_gif_path(positions, primary_springs , boundaries=config.BOUNDARIES,  filename="low_res_21.gif")
    #save_gif_smooth(positions, config.BOUNDARIES, filename="whole_smooth_path.gif", fps=30, line_width=5, s=0.5, num_curve_points=200)
    save_gif_smooth_gradient(positions, boundaries=config.BOUNDARIES, 
                             filename='convex_pyramid_path.gif', fps=30, 
                             s=0.5, num_curve_points=200, line_width=5)
    #bezier_gif_path(positions, primary_springs, boundaries=config.BOUNDARIES, filename="smoothed_path.gif", fps=30, line_width=5, curve_points=10, offset_magnitude=-10.5)

if __name__ == '__main__':
    main()
