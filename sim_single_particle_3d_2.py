import concurrent.futures
import csv
import multiprocessing as mp
import multiprocessing.pool
from functools import partial

import numpy as np
import matplotlib.pyplot
import math

# 1) Vectors need to be accurately described
# 2) F=qv x B needs to be incorporated
# 3) Data visualization
# 4) More particles (with particle class)

def integrate_b_field(field, x_range, y_range, z_range, step_size = 0.1/100):
    # For each discrete square/cube, the line integral of the perimeter is equal to the change in the B_field within the same area
    # ∇xE = -dB/dt
    # E = integral(-dB/dt)
    # We need to integrate in every direction (6 squares facing each direction in a cube)
    dz_0 = 0
    dz_1 = 0
    for x in np.linspace(x_range[0], x_range[1], math.ceil((x_range[1] - x_range[0]) / step_size)):
        for y in np.linspace(y_range[0], y_range[1], math.ceil((y_range[1] - y_range[0]) / step_size)):
            # X-Y Face (iterate through x/y,x+1/y,x+1/y+1, x/y+1) - Looking at Z-axis of dB/dt @z=0
            dz_0 += -field(x, y, z_range[0])[2]
            # X-Y Face (iterate through x/y,x+1/y,x+1/y+1, x/y+1) - Looking at Z-axis of dB/dt @z=1
            dz_1 += -field(x, y, z_range[1])[2]

    dx_0 = 0
    dx_1 = 0
    for y in np.linspace(y_range[0], y_range[1], math.ceil((y_range[1] - y_range[0]) / step_size)):
        for z in np.linspace(z_range[0], z_range[1], math.ceil((z_range[1] - z_range[0]) / step_size)):
            # Y-Z Face (iterate through y/z, y+1/z, y+1/z+1, y/z+1) - Looking at X-axis of dB/dt @x=0
            dx_0 += -field(x_range[0], y, z)[0]
            # Y-Z Face (iterate through y/z, y+1/z, y+1/z+1, y/z+1) - Looking at X-axis of dB/dt @x=1
            dx_1 += -field(x_range[1], y, z)[0]

    dy_0 = 0
    dy_1 = 0
    for x in np.linspace(x_range[0], x_range[1], math.ceil((x_range[1] - x_range[0])/step_size)):
        for z in np.linspace(z_range[0], z_range[1], math.ceil((z_range[1] - z_range[0])/step_size)):
            # X-Z Face (iterate through x/z, x+1/z, x+1/z+1, x/z+1) - Looking at Y-axis of dB/dt @y=0
            dy_0 += -field(x, y_range[0], z)[1]
            # X-Z Face (iterate through x/z, x+1/z, x+1/z+1, x/z+1) - Looking at Y-axis of dB/dt @y=1
            dy_1 += -field(x, y_range[1], z)[1]

    return #TODO: HOW TO COMBINE VECTORS?

PI = 3.141592

# 0,0,0 is the bottom-left corner, closest to the screen
# X (WIDTH), Y (HEIGHT), Z (DEPTH)
# Units in meters
WIDTH = 1
HEIGHT = 1
DEPTH = 1

### PARTICLE DEFINITION
# TODO: Turn this into a class
# Particle will start in the center of the simulated area on the left hand side of the screen and move towards the right
particle_position = np.array([0.5, 0.5, 0.5])
particle_velocity = 0
particle_mass = 1.6726e-27      # kg
particle_charge = 1.602e-19     # Coulombs

### B FIELD DEFINITION - position (meters), field frequency (Hz) time (s)
# Magnetic field will be defined via continuous equations
FIELD_AMPLITUDE = 1e-6#5         # Amps/Meter (H), Magnetic Field flux (B) is related to (H) via ∇B = u∇H
MU_0 = 4e-7*PI              # H/m == T*m/A
FIELD_FREQUENCY = 200000      # Hz
B_STEP_SIZE = 0.1/100    # Must be significantly smaller than E_STEP_SIZE
def B_field(x, y, z, w, t):
    return math.sin(2 * PI * t * w) * FIELD_AMPLITUDE * x * MU_0 * np.array([[1], [0], [0]])  # Teslas (V*s/m^2)

def dt_B_field(x, y, z, w, t):
    return 2*PI*w*math.cos(2*PI*t*w) * FIELD_AMPLITUDE * x * MU_0 * np.array([[1], [0], [0]])      # T/s (V/m^2)

### E FIELD DEFINITION
# Electric field will be discretized
E_STEP_SIZE = 0.1
E_field = np.ndarray([math.ceil(WIDTH / E_STEP_SIZE), math.ceil(HEIGHT / E_STEP_SIZE), math.ceil(DEPTH / E_STEP_SIZE)], dtype=np.ndarray)

if __name__ == "__main__":
    # Data File
    f = open("sim_output_firstattempt.csv", 'w', newline='')
    writer = csv.DictWriter(f, fieldnames=['time', 'force', 'accel', 'velocity', 'position'])
    writer.writeheader()
    writer.writerow({'time': 0, 'force': 0, 'accel': 0, 'velocity': particle_velocity,
                     'position': particle_position[0]})

    ### SIMULATION LOOP
    NUM_STEPS = 100
    TIME_STEP_SIZE = 0.0000001   # seconds
    time = 0
    for i in range(NUM_STEPS):

        test = []
        future = np.ndarray([math.ceil(WIDTH / E_STEP_SIZE), math.ceil(HEIGHT / E_STEP_SIZE), math.ceil(DEPTH / E_STEP_SIZE)], dtype=multiprocessing.pool.ApplyResult)
        with mp.Pool() as pool:
            for x, y, z in np.ndindex(E_field.shape):
                future[x][y][z] = pool.apply_async(integrate_b_field, [partial(dt_B_field, w=FIELD_FREQUENCY, t=time), (x, x + E_STEP_SIZE), (y, y + E_STEP_SIZE), (z, z + E_STEP_SIZE), B_STEP_SIZE])#, callback=cb)
            for x, y, z in np.ndindex(E_field.shape):
                future[x][y][z].wait()
                E_field[x][y][z] = future[x][y][z].get()


                #E_field[x][y][z] = integrate_b_field(partial(dt_B_field, w=FIELD_FREQUENCY, t=time), (x,x+E_STEP_SIZE), (y, y+E_STEP_SIZE), (z, z+E_STEP_SIZE), step_size=B_STEP_SIZE)
        time += TIME_STEP_SIZE
        #print(E_field)
        print(f"Timestep: {time}")
        # TODO: Add vectorized math
        force = particle_charge*E_field[math.floor(particle_position[0]*WIDTH/E_STEP_SIZE)][math.floor(particle_position[1]*HEIGHT/E_STEP_SIZE)][math.floor(particle_position[2]*DEPTH/E_STEP_SIZE)]
        accel = force/particle_mass
        particle_velocity += accel*TIME_STEP_SIZE
        particle_position[0] += particle_velocity*TIME_STEP_SIZE
        print(f"Particle Position: {particle_position}")
        print(f"Particle Accel: {accel}")
        print(f"Particle Velocity: {particle_velocity}")
        writer.writerow({'time': (i+1)*TIME_STEP_SIZE, 'force': force, 'accel': accel, 'velocity': particle_velocity,
                         'position': particle_position[0]})

        # F = qvB (q known, v known, B known) - TODO
        # F = qE
        # F = ma
        # a = dv/dt
        # v = dPosition/dt

