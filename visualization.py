# visualization.py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
#from config import INITIAL_LAYERS, BOUNDARIES
import config
from scipy.interpolate import splprep, splev
from mpl_toolkits.mplot3d.art3d import Line3DCollection


def save_gif(positions, springs, boundaries=config.BOUNDARIES, filename="path_balls_with_springs.gif", fps=30):
    """
    Animates the particle positions and draws springs as red lines connecting the particles.
    'positions' is an array of shape (num_steps, num_particles, 3).
    'springs' is a list of Spring objects with attributes particle1 and particle2.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(boundaries['x'])
    ax.set_ylim(boundaries['y'])
    ax.set_zlim(boundaries['z'])
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(filename)
    
    # Create a scatter plot for particles.
    scatter = ax.scatter([], [], [], color='blue', s=100)
    
    # Create a line (artist) for each spring from the springs list.
    spring_lines = []
    for _ in springs:
        line, = ax.plot([], [], [], color='red', lw=2)
        spring_lines.append(line)
    
    def update(frame):
        pos = positions[frame]  # shape: (num_particles, 3)
        # Update particle positions.
        scatter._offsets3d = (pos[:, 0], pos[:, 1], pos[:, 2])
        
        # For each spring, update the corresponding line based on the connected particles.
        for i, spring in enumerate(springs):
            idx1 = spring.particle1
            idx2 = spring.particle2
            x_vals = [pos[idx1, 0], pos[idx2, 0]]
            y_vals = [pos[idx1, 1], pos[idx2, 1]]
            z_vals = [pos[idx1, 2], pos[idx2, 2]]
            spring_lines[i].set_data(x_vals, y_vals)
            spring_lines[i].set_3d_properties(z_vals)
        return [scatter] + spring_lines
    
    ani = FuncAnimation(fig, update, frames=len(positions), interval=20, blit=True)
    writer = PillowWriter(fps=fps)
    ani.save(filename, writer=writer)
    print(f"Saved animation as {filename}")

def save_gif_path(positions, springs, boundaries=config.BOUNDARIES, 
                  filename="path_with_color.gif", fps=30, line_width=5):
    """
    Animates the path (primary springs) in 3D and saves it as a GIF.
    
    Parameters:
      positions: numpy array of shape (num_steps, num_particles, 3)
      springs: list of Spring objects (each with attributes particle1 and particle2)
      boundaries: dict with limits for 'x', 'y', and 'z'
      filename: output filename
      fps: frames per second for the GIF
      line_width: thickness of the spring lines
      
    The color of each spring line is determined by the average z value of its endpoints,
    normalized between the bottom and top boundaries.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(boundaries['x'])
    ax.set_ylim(boundaries['y'])
    ax.set_zlim(boundaries['z'])
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(filename)
    
    # Pre-create one line (artist) for each spring.
    spring_lines = []
    for _ in springs:
        line, = ax.plot([], [], [], lw=line_width)
        spring_lines.append(line)
    
    def update(frame):
        pos = positions[frame]  # shape: (num_particles, 3)
        # For each spring, update the line based on the actual connection.
        for i, spring in enumerate(springs):
            idx1 = spring.particle1
            idx2 = spring.particle2
            x_vals = [pos[idx1, 0], pos[idx2, 0]]
            y_vals = [pos[idx1, 1], pos[idx2, 1]]
            z_vals = [pos[idx1, 2], pos[idx2, 2]]
            spring_lines[i].set_data(x_vals, y_vals)
            spring_lines[i].set_3d_properties(z_vals)
            
            # Compute the average z of the two endpoints.
            avg_z = (z_vals[0] + z_vals[1]) / 2.0
            # Normalize the average z with respect to the boundaries.
            norm = (avg_z - boundaries['z'][0]) / (boundaries['z'][1] - boundaries['z'][0])
            # Create a grayscale color: darker at the bottom (norm ~ 0), lighter at the top (norm ~ 1).
            norm = max(0, min(1, norm))  # Clamp norm to [0, 1]
            # Create a grayscale color.
            color = (norm, norm, norm)
            
            spring_lines[i].set_color(color)
        return spring_lines
    
    ani = FuncAnimation(fig, update, frames=len(positions), interval=20, blit=True)
    writer = PillowWriter(fps=fps)
    ani.save(filename, writer=writer)
    print(f"Saved animation as {filename}")


def save_layer_gif_path(positions, springs, boundaries=config.BOUNDARIES, 
                          filename="path_with_color.gif", fps=30, line_width=5):
    """
    Animates the path (springs) in 3D and saves it as a GIF.
    Each spring is colored based on the fixed (initial) layer of its endpoints.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(boundaries['x'])
    ax.set_ylim(boundaries['y'])
    ax.set_zlim(boundaries['z'])
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(filename)
    
    # Create a line (artist) for each spring.
    spring_lines = []
    for _ in springs:
        line, = ax.plot([], [], [], lw=line_width)
        spring_lines.append(line)
    
    # Use global INITIAL_LAYERS.
    max_layer = max(config.INITIAL_LAYERS) if config.INITIAL_LAYERS else 1
    
    def update(frame):
        pos = positions[frame]  # shape: (num_particles, 3)
        for i, spring in enumerate(springs):
            idx1 = spring.particle1
            idx2 = spring.particle2
            x_vals = [pos[idx1, 0], pos[idx2, 0]]
            y_vals = [pos[idx1, 1], pos[idx2, 1]]
            z_vals = [pos[idx1, 2], pos[idx2, 2]]
            spring_lines[i].set_data(x_vals, y_vals)
            spring_lines[i].set_3d_properties(z_vals)
            
            # Use the fixed layer info from the global variable.
            avg_layer = (config.INITIAL_LAYERS[idx1] + config.INITIAL_LAYERS[idx2]) / 2.0
            norm = avg_layer / max_layer
            norm = max(0, min(1, norm))
            color = (norm, 0, norm)
            spring_lines[i].set_color(color)
        return spring_lines
    
    ani = FuncAnimation(fig, update, frames=len(positions), interval=20, blit=True)
    writer = PillowWriter(fps=fps)
    ani.save(filename, writer=writer)
    print(f"Saved animation as {filename}")


def bezier_curve(p0, p1, cp, num_points=10):
    """
    Computes a quadratic Bézier curve between p0 and p1 with control point cp.
    Returns an array of shape (num_points, 3).
    """
    t_vals = np.linspace(0, 1, num_points)[:, None]  # column vector
    curve = (1 - t_vals)**2 * p0 + 2 * (1 - t_vals) * t_vals * cp + t_vals**2 * p1
    return curve

def smooth_polyline(points, s=0.5, num_points=200):
    """
    Given an array of points (shape (N,3)), fit a spline and generate a smooth curve.
    Parameters:
      points: (N,3) array of (x,y,z) coordinates.
      s: smoothing factor (0 means interpolation through all points).
      num_points: number of points in the resulting smooth curve.
    Returns:
      A (num_points,3) array representing the smooth polyline.
    """
    points = np.array(points)
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    tck, u = splprep([x, y, z], s=s)
    u_fine = np.linspace(0, 1, num_points)
    x_new, y_new, z_new = splev(u_fine, tck)
    return np.column_stack((x_new, y_new, z_new))

def save_gif_smooth_gradient(positions, boundaries=config.BOUNDARIES, 
                             filename="smooth_gradient_path.gif", fps=30, 
                             s=0.5, num_curve_points=200, line_width=5):
    """
    Animates a smooth polyline for each frame using all particle positions,
    and applies a gradient (from dark to light) along the curve.
    
    Parameters:
      positions: numpy array of shape (num_steps, num_particles, 3)
      boundaries: dict with limits for 'x', 'y', and 'z'
      filename: output filename
      fps: frames per second for the GIF
      s: smoothing factor for the spline
      num_curve_points: number of points in the smooth curve
      line_width: thickness of the drawn curve
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(boundaries['x'])
    ax.set_ylim(boundaries['y'])
    ax.set_zlim(boundaries['z'])
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(filename)
    
    # We'll update a Line3DCollection each frame.
    # Create an initial empty collection.
    lc = Line3DCollection([], linewidths=line_width)
    ax.add_collection3d(lc)
    
    cmap = plt.get_cmap('spring')  # Use a grayscale colormap.
    
    def update(frame):
        pos = positions[frame]  # shape: (num_particles, 3)
        # Generate a smooth curve through these points.
        smooth_line = smooth_polyline(pos, s=s, num_points=num_curve_points)
        # Create segments: each segment is a pair of points.
        segments = [smooth_line[i:i+2] for i in range(len(smooth_line)-1)]
        # Create a normalized parameter for each segment.
        t_vals = np.linspace(0, 1, len(segments))
        # Use the colormap to get colors (darker at t=0, lighter at t=1).
        colors = cmap(t_vals)
        # Update the Line3DCollection.
        lc.set_segments(segments)
        lc.set_color(colors)
        return [lc]
    
    ani = FuncAnimation(fig, update, frames=len(positions), interval=20, blit=True)
    writer = PillowWriter(fps=fps)
    ani.save(filename, writer=writer)
    print(f"Saved animation as {filename}")

def smooth_polyline(points, s=0.5, num_points=200):
    """
    Given an array of points (shape (N,3)), fit a spline and generate a smooth curve.
    
    Parameters:
      points: (N,3) array of (x,y,z) coordinates.
      s: smoothing factor (0 means interpolation through all points).
      num_points: number of points in the resulting smooth curve.
      
    Returns:
      (num_points, 3) array representing the smooth polyline.
    """
    points = np.array(points)
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    # Fit a spline to the points.
    tck, u = splprep([x, y, z], s=s)
    u_fine = np.linspace(0, 1, num_points)
    x_new, y_new, z_new = splev(u_fine, tck)
    return np.column_stack((x_new, y_new, z_new))

def save_gif_smooth(positions, boundaries=config.BOUNDARIES, filename="smooth_path.gif", fps=30, line_width=5, s=0.5, num_curve_points=200):
    """
    Animates a smooth polyline for each frame using all particle positions.
    
    Instead of drawing each spring individually, this function takes the complete
    set of particle positions (assumed to be in order along the primary path), applies
    spline interpolation to create a smooth curve, and then draws that curve.
    
    Parameters:
      positions: numpy array of shape (num_steps, num_particles, 3)
      boundaries: dict with limits for 'x', 'y', and 'z'
      filename: output filename
      fps: frames per second for the GIF
      line_width: thickness of the drawn curve
      s: smoothing factor for the spline (adjust for desired smoothness)
      num_curve_points: number of points for the smooth curve
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(boundaries['x'])
    ax.set_ylim(boundaries['y'])
    ax.set_zlim(boundaries['z'])
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Smoothed Path Visualization")
    
    # Create a single line artist for the smoothed path.
    line, = ax.plot([], [], [], color='red', lw=line_width)
    
    def update(frame):
        pos = positions[frame]  # shape: (num_particles, 3)
        # Generate a smooth curve through all particle positions.
        smooth_line = smooth_polyline(pos, s=s, num_points=num_curve_points)
        line.set_data(smooth_line[:, 0], smooth_line[:, 1])
        line.set_3d_properties(smooth_line[:, 2])
        return [line]
    
    ani = FuncAnimation(fig, update, frames=len(positions), interval=20, blit=True)
    writer = PillowWriter(fps=fps)
    ani.save(filename, writer=writer)
    print(f"Saved animation as {filename}")

def bezier_gif_path(positions, springs, boundaries, filename="smoothed_path.gif", fps=30, line_width=5, curve_points=10, offset_magnitude=0.5):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(boundaries['x'])
    ax.set_ylim(boundaries['y'])
    ax.set_zlim(boundaries['z'])
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Path Simulation with Smoothed Springs")
    
    # Create a line for each spring.
    spring_lines = []
    for _ in springs:
        line, = ax.plot([], [], [], lw=line_width)
        spring_lines.append(line)
    
    def update(frame):
        pos = positions[frame]  # shape: (num_particles, 3)
        # Update each spring using a Bézier curve.
        for i, spring in enumerate(springs):
            idx1 = spring.particle1
            idx2 = spring.particle2
            p0 = pos[idx1]
            p1 = pos[idx2]
            # Compute midpoint.
            midpoint = (p0 + p1) / 2.0
            # Choose an offset (for example, a vertical offset):
            offset = np.array([0, 0, offset_magnitude])
            cp = midpoint + offset
            curve = bezier_curve(p0, p1, cp, num_points=curve_points)
            spring_lines[i].set_data(curve[:, 0], curve[:, 1])
            spring_lines[i].set_3d_properties(curve[:, 2])
        return spring_lines
    
    from matplotlib.animation import FuncAnimation, PillowWriter
    ani = FuncAnimation(fig, update, frames=len(positions), interval=20, blit=True)
    writer = PillowWriter(fps=fps)
    ani.save(filename, writer=writer)
    print(f"Saved animation as {filename}")

