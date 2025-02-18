# config.py
#higher subdivisions need to change the subsettings
BOUNDARIES = {'x': (-20, 20), 'y': (-20, 20), 'z': (0, 40)}
INITIAL_LAYERS = None  # This will be updated by assign_layers.
NEIGHBOUR_RADIUS = 10#10 # affects secondary springs, between which neighbours, low res is 12 points per circle, but might need to change depending on the geometry

AGE_FACTOR= 7 #7

STIFFNESS_AGE_FACTOR = 0.1
PLASTICITY_AGE_FACTOR = 0.05
YIELD_AGE_FACTOR = 0.1

STIFFNESS_MIN = 0.5#0.05#
STIFFNESS_MAX = 2.0

PLASTICITY_MIN = 0.05
PLASTICITY_MAX = 0.9 #2.0 #0.9

YIELD_RATIO_MIN = 0.05#0.5
YIELD_RATIO_MAX = 1.5 #2.0 #1.5

MAX_ENERGY = 30000#50000

