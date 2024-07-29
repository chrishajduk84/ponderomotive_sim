import math
import random

from constants import ELEMENTARY_CHARGE, SOLAR_WIND_DENSITY, SOLAR_WIND_SPEED, PROTON_MASS, ELECTRON_MASS, \
    VACUUM_PERMITIVITY, PI
from particle import Particle2D, FixedParticle2D
from matplotlib import pyplot as plt
import numpy as np

def particle_cleaner(all_particles, x_limits, y_limits):
    to_delete = []
    for i, _ in enumerate(all_particles):
        p = all_particles[i]
        if p.x > x_limits[1] or p.x < x_limits[0] or p.y > y_limits[1] or p.y < y_limits[0]:
            to_delete.append(i)
    to_delete = reversed(to_delete)
    for i in to_delete:
        del all_particles[i]

def particle_generator(all_particles, x_limits, y_limits, time_step, velocity=SOLAR_WIND_SPEED, density=SOLAR_WIND_DENSITY):
    # For every time_step a particle travels a certain distance.
    minimum_distance = time_step * velocity
    # The minimum distance travelled by a particle with `velocity` represents the area within which new particles can be generated
    generated_area = minimum_distance*(y_limits[1] - y_limits[0])   # In m^2
    # Within this generated area we expect a certain number of protons according to the density
    num_protons = generated_area * density

    # Each of the protons will be randomly distributed within  the generated_area
    # Each of the protons will also have a corresponding electron because the plasma must be quasi-neutral
    for i in range(math.ceil(num_protons)):
        x = random.random() * minimum_distance
        y = random.random() * (y_limits[1] - y_limits[0]) + y_limits[0]
        all_particles.append(Particle2D(x, y, velocity, 0, ELEMENTARY_CHARGE, PROTON_MASS))
        all_particles.append(Particle2D(x, y, velocity, 0, -ELEMENTARY_CHARGE, ELECTRON_MASS))


# Particles coming from -X -> +X
# Initialize limits
x_limits = [0, 0.1]
y_limits = [0, 1]

# Initialize Simulation Objects
fixed_particles = []
sim_particles = []

fig = plt.figure()
ax = fig.gca()

# Initialize electrostatic grid - Pattern: every 25cm
y_spacing = np.linspace(y_limits[0], y_limits[1], 100)
for i in y_spacing:
    fixed_particles.append(FixedParticle2D((x_limits[1] - x_limits[0])/2, i, ELEMENTARY_CHARGE))
for i in y_spacing:
    fixed_particles.append(FixedParticle2D((x_limits[1] - x_limits[0])/2+0.01, i, -ELEMENTARY_CHARGE))

NUM_TIME_STEPS = 10000
TIME_STEP = (x_limits[1] - x_limits[0])/SOLAR_WIND_SPEED/NUM_TIME_STEPS
for t in range(NUM_TIME_STEPS):
    if t < 2:
        particle_generator(sim_particles, x_limits, y_limits, TIME_STEP)

    if t % 20 == 0:
        for f in fixed_particles:
            if f.charge != 0:
                f.charge = 0
            else:
                if f.x < 0.055:
                    f.charge = ELEMENTARY_CHARGE
                else:
                    f.charge = -ELEMENTARY_CHARGE

    for p in sim_particles:
        force_x = 0.0
        force_y = 0.0
        for f in fixed_particles:
            rx = p.x - f.x
            ry = p.y - f.y
            force_x += 1/(4*PI*VACUUM_PERMITIVITY)*p.charge*f.charge/(rx**2) * (-1 if rx < 0 else 1)
            force_y += 1/(4*PI*VACUUM_PERMITIVITY)*p.charge*f.charge/(ry**2) * (-1 if ry < 0 else 1)
        accel_x = force_x/p.mass
        accel_y = force_y/p.mass
        p.vx += accel_x*TIME_STEP
        p.vy += accel_y*TIME_STEP
        p.x += p.vx*TIME_STEP
        p.y += p.vy*TIME_STEP
        # for q in sim_particles:
        #     if p == q:
        #         continue
        #     else:
        #         rx = p.x - q.x
        #         ry = p.y - q.y
        #         if rx != 0:
        #             force_x = 1 / (4 * PI * VACUUM_PERMITIVITY) * p.charge * q.charge / (rx ** 2) * (-1 if rx < 0 else 1)
        #             accel_x = force_x / p.mass
        #             p.vx += accel_x * TIME_STEP
        #             p.x += p.vx * TIME_STEP
        #
        #             #q
        #             accel_x = -force_x / q.mass
        #             q.vx += accel_x * TIME_STEP
        #             q.x += q.vx * TIME_STEP
        #
        #         if ry != 0:
        #             force_y = 1 / (4 * PI * VACUUM_PERMITIVITY) * p.charge * q.charge / (ry ** 2) * (-1 if ry < 0 else 1)
        #             accel_y = force_y / p.mass
        #             p.vy += accel_y * TIME_STEP
        #             p.y += p.vy * TIME_STEP
        #
        #             # For Q
        #             accel_y = -force_y / q.mass
        #             q.vy += accel_y * TIME_STEP
        #             q.y += q.vy * TIME_STEP

    particle_cleaner(sim_particles, x_limits, y_limits)

    plt.cla()
    ax.scatter([p.x for p in fixed_particles], [p.y for p in fixed_particles], color='y')
    sim_particles_positive = []
    sim_particles_negative = []
    for p in sim_particles:
        print(f"{p.x}, {p.vx}")
        if p.charge > 0:
            sim_particles_positive.append(p)
        else:
            sim_particles_negative.append(p)
    ax.scatter([p.x for p in sim_particles_positive], [p.y for p in sim_particles_positive], color='r')
    ax.scatter([p.x for p in sim_particles_negative], [p.y for p in sim_particles_negative], color='b')
    ax.set_xlim([0, 0.1])
    plt.pause(0.01)
    plt.show(block=False)
    print(f"TIME STEP: {t}")
