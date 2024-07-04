import math
import random
import matplotlib.pyplot as plt

from particle import Particle1D

# Start with a random continuous stream of particles travelling at a constant uniform velocity
v = 400*1000            # 400km/s => 400 m/s
#p = 5 / (0.01 ^ 3)   # 5 particles per cm^3 (Protons) or 5/(0.01^3) particles per m^3
linear_p = 5/0.00001          # For 1d case lets put all the particles that would normally fill a volume into a 1d-line
# Since the solar wind/plasma is quasi-neutral, the protons are in charge pairs

spatial_resolution = 0.00001                   # meters
temporal_resolution = spatial_resolution/v/10  # seconds

SPATIAL_BOUNDARY = (0, 1)
particles_in_boundaries = []    # Probably better to convert to linked list later

f_field_1 = lambda q1, x: 1/(4*3.14159256*8.854e-12)*q1*1.602176e-19/((x-0.01)**2) * (-1 if x < 0.01 else 1)   #V/m
f_field_2 = lambda q1, x: 1/(4*3.14159256*8.854e-12)*q1*-1.602176e-19/((x-0.02)**2) * (-1 if x < 0.02 else 1)   #V/m

fig = plt.figure()
ax = fig.gca()

TIME_STEPS = 100000
for t_step in range(TIME_STEPS):
    print(f"Time stamp: {t_step}")
    if t_step < 100:
        # Generate particles flying in from 0m -> 1m on the x-axis
        for i in range(math.floor(linear_p*spatial_resolution)):
            p = Particle1D(x=random.random()*spatial_resolution, v=v, charge=1.602176e-19, mass=1.672e-27)
            e = Particle1D(x=random.random()*spatial_resolution, v=v, charge=-1.602176e-19, mass=9.109e-31)
            particles_in_boundaries.append(p)
            particles_in_boundaries.append(e)

    # Simulate change in speed and position under an electrostatic field
    indices_to_remove = []
    to_plot_x = [part.x for part in particles_in_boundaries]
    to_plot_y = [part.charge for part in particles_in_boundaries]
    ax.scatter(to_plot_x, to_plot_y)
    for i in range(len(particles_in_boundaries)):
        # F = qE
        F = f_field_1(particles_in_boundaries[i].charge, particles_in_boundaries[i].x) + \
            f_field_2(particles_in_boundaries[i].charge, particles_in_boundaries[i].x)
        # Update acceleration
        a = F/particles_in_boundaries[i].mass
        # Update velocity
        particles_in_boundaries[i].v += a*temporal_resolution
        # Update position
        particles_in_boundaries[i].x += particles_in_boundaries[i].v*temporal_resolution

        # If they have exited the simulation boundary, lets mark them for removal (can't remove while we are iterating)
        if particles_in_boundaries[i].x > SPATIAL_BOUNDARY[1]:
            indices_to_remove.append(i)
        to_plot_x.append(particles_in_boundaries[i].x)
        to_plot_y.append(particles_in_boundaries[i].charge)
    plt.cla()
    ax.scatter(to_plot_x, to_plot_y)
    ax.set_xlim([0, 0.025])
    plt.pause(0.01)
    plt.show(block=False)

    for i in reversed(indices_to_remove):
        particles_in_boundaries.pop(i)




if __name__ == "__main__":
    pass