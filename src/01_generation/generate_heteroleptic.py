import math
import random
import numpy as np
import gsd.hoomd

# parameters
N_star = 1000
f = 6
n_m = 16
L = 750.0

bond_length = 1.2
terminal_bond_length = 0.45
min_dist = 0.85

sticky_flag = 1
terminal_diameter = 0.1
normal_diameter = 1.0

frac_C = 0.3

seed = 12345
max_trials_bead = 10000

random.seed(seed)
np.random.seed(seed)

N_m = 1 + f * n_m
N_particles = N_star * N_m
outfile = "Initial_config_N_" + str(N_particles) + "_L_" + str(int(L)) + ".gsd"

N_C_star = int(round(frac_C * N_star))
N_A_star = N_star - N_C_star


def random_unit_vector():
    u = 2.0 * random.random() - 1.0
    phi = 2.0 * math.pi * random.random()
    s = math.sqrt(1.0 - u * u)
    return np.array([s * math.cos(phi), s * math.sin(phi), u])


def too_close(pos, star_positions):
    for i in range(len(star_positions) - 1):
        dr = pos - star_positions[i]
        if np.dot(dr, dr) < min_dist * min_dist:
            return True
    return False


# safe sphere radius around each star center
R_max = (n_m - 1) * bond_length + terminal_bond_length
R_safe = R_max + 1.0

# cubic lattice of cell centers, one star per cell
ncell_1d = math.ceil(N_star ** (1.0 / 3.0))
cell_size = L / ncell_1d

cell_centers = []
for ix in range(ncell_1d):
    for iy in range(ncell_1d):
        for iz in range(ncell_1d):
            cx = -0.5 * L + (ix + 0.5) * cell_size
            cy = -0.5 * L + (iy + 0.5) * cell_size
            cz = -0.5 * L + (iz + 0.5) * cell_size
            cell_centers.append(np.array([cx, cy, cz]))

random.shuffle(cell_centers)

positions = np.zeros((N_particles, 3))
typeid = np.zeros(N_particles, dtype=np.uint32)
diameter = np.full(N_particles, normal_diameter)

bonds_group = []
bonds_typeid = []
angles_group = []
angles_typeid = []

particle_types = ['A', 'B', 'C', 'D']
bond_types = ['A-A', 'A-B', 'C-C', 'C-D']
angle_types = ['A-A-A', 'A-O-A', 'C-C-C', 'C-O-C']

# 0 = A-star, 1 = C-star
star_labels = [0] * N_A_star + [1] * N_C_star
random.shuffle(star_labels)

global_index = 0

for s in range(N_star):
    if star_labels[s] == 1:
        backbone_type = 2
        sticky_type = 3
        backbone_bond_type = 2
        sticky_bond_type = 3
        chain_angle_type = 2
        core_angle_type = 3
    else:
        backbone_type = 0
        sticky_type = 1
        backbone_bond_type = 0
        sticky_bond_type = 1
        chain_angle_type = 0
        core_angle_type = 1

    center = cell_centers[s]

    core_idx = global_index
    positions[core_idx] = center
    typeid[core_idx] = backbone_type
    diameter[core_idx] = normal_diameter
    global_index = global_index + 1

    star_positions = [center]
    arm_first_indices = []

    for arm in range(f):
        prev_idx = core_idx
        prev_pos = positions[prev_idx]
        first_idx_this_arm = None

        for m in range(n_m):
            bead_idx = global_index
            is_terminal = (m == n_m - 1)

            if is_terminal:
                step_len = terminal_bond_length
            else:
                step_len = bond_length

            for trial in range(max_trials_bead):
                direction = random_unit_vector()
                trial_pos = prev_pos + step_len * direction

                if np.linalg.norm(trial_pos - center) > R_safe:
                    continue

                if not is_terminal and too_close(trial_pos, star_positions):
                    continue

                positions[bead_idx] = trial_pos

                if sticky_flag == 1 and is_terminal:
                    typeid[bead_idx] = sticky_type
                    diameter[bead_idx] = terminal_diameter
                    bonds_typeid.append(sticky_bond_type)
                else:
                    typeid[bead_idx] = backbone_type
                    diameter[bead_idx] = normal_diameter
                    bonds_typeid.append(backbone_bond_type)

                bonds_group.append((prev_idx, bead_idx))

                if m == 0:
                    first_idx_this_arm = bead_idx

                if m == 1:
                    angles_group.append((core_idx, bead_idx - 1, bead_idx))
                    angles_typeid.append(chain_angle_type)
                elif m >= 2:
                    angles_group.append((bead_idx - 2, bead_idx - 1, bead_idx))
                    angles_typeid.append(chain_angle_type)

                prev_idx = bead_idx
                prev_pos = trial_pos
                star_positions.append(trial_pos.copy())
                global_index = global_index + 1
                break

        arm_first_indices.append(first_idx_this_arm)

    for arm in range(f):
        a = arm_first_indices[arm]
        b = core_idx
        c = arm_first_indices[(arm + 1) % f]
        angles_group.append((a, b, c))
        angles_typeid.append(core_angle_type)

frame = gsd.hoomd.Frame()
frame.particles.N = N_particles
frame.particles.position = positions
frame.particles.types = particle_types
frame.particles.typeid = typeid
frame.particles.diameter = diameter
frame.configuration.box = [L, L, L, 0, 0, 0]

frame.bonds.N = len(bonds_group)
frame.bonds.types = bond_types
frame.bonds.typeid = np.array(bonds_typeid, dtype=np.uint32)
frame.bonds.group = np.array(bonds_group, dtype=np.int32)

frame.angles.N = len(angles_group)
frame.angles.types = angle_types
frame.angles.typeid = np.array(angles_typeid, dtype=np.uint32)
frame.angles.group = np.array(angles_group, dtype=np.int32)

with gsd.hoomd.open(name=outfile, mode='x') as traj:
    traj.append(frame)

print("Wrote", outfile)
print("N_particles =", N_particles)
print("N_star =", N_star, "N_m =", N_m)
print("N_A_star =", N_A_star, "N_C_star =", N_C_star)
