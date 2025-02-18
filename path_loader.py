# path_loader.py
import csv
import numpy as np

def load_path(filename):
    points = []
    with open(filename, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            points.append((float(row['x']), float(row['y']), float(row['z'])))
    return points

def compute_rest_lengths(points):
    rest_lengths = []
    for i in range(len(points) - 1):
        p0 = np.array(points[i])
        p1 = np.array(points[i + 1])
        d = np.linalg.norm(p1 - p0)
        rest_lengths.append(d)
    return rest_lengths
