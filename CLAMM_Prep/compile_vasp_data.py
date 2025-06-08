import numpy as np
from numpy.linalg import norm
import os
import copy


def cart_to_frac(cart_coord, basis):
    """
    transform from Cartesian coordinate to direct/fraction coordinate
    :param cart_coord: list in the format of [0, 0, 0]
    :param basis: list in the format of [[0.0, 0.0, 3.6], [0.0, 3.6, 0.0], [3.6, 0.0, 0.0]]
    :return: frac_coord
    """
    a1 = np.array(basis[0])
    a2 = np.array(basis[1])
    a3 = np.array(basis[2])
    cart_coord = np.array(cart_coord)
    trans_matr = np.vstack([a1, a2, a3]).T
    inv_matr = np.linalg.inv(trans_matr)
    frac_coord = np.matmul(inv_matr, cart_coord.T).T

    return list(frac_coord)


def apply_pbc(coord):
    """
    Apply periodic boundary conditions to a point
    :param coord: fraction coordinate for a given point
    :return: fraction coordinates using PBCs
    """
    pbc_coord = copy.deepcopy(coord)
    for i in range(3):
        pbc_coord[i] = np.around(pbc_coord[i], decimals=15) % 1
        if abs(pbc_coord[i] - 1) < 0.1:
            pbc_coord[i] = 0

    return pbc_coord


def read_lat(contcar_lines):
    """
    read lattice constant and lattice vector from contcar/poscar file
    :param contcar_lines: read lines of contcar/poscar
    :return: lattice info from contcar/poscar
    """
    scale = float(contcar_lines[1])
    lat_const = []
    lat_vec = []
    lat_ang = []
    for i in range(2, 5):
        vec = np.array([float(x) for x in contcar_lines[i].split()]) * scale
        lat_vec.append(vec)
        lat_const.append(norm(vec))
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        cos_ = np.dot(lat_vec[j], lat_vec[k]) / (lat_const[j] * lat_const[k])
        lat_ang.append(np.arccos(np.clip(cos_, -1, 1)))

    return lat_const, lat_vec, lat_ang


def read_poscar_coords(poscar_lines):
    """
    read species and composition from atat-formatted files
    :param poscar_lines: read lines of atat-formatted poscar files, i.e. no species information
    :return: frac coords from poscar
    """
    scale = float(poscar_lines[1])
    lat_vec = []
    for i in range(2, 5):
        vec = np.array([float(x) for x in poscar_lines[i].split()]) * scale
        lat_vec.append(vec)
    pos_start = 0
    for i in range(len(poscar_lines)):
        if 'artesian' in poscar_lines[i] or 'irect' in poscar_lines[i]:
            pos_start = i + 1
    atom_sum = len(poscar_lines) - pos_start
    pos_list = [[]] * atom_sum
    for i in range(atom_sum):
        pnt = poscar_lines[pos_start + i].split()
        pos = []
        for j in range(3):
            pos.append(float(pnt[j]))
        if 'artesian' in poscar_lines[pos_start - 1]:
            pos = cart_to_frac(pos, lat_vec)
        pos = apply_pbc(pos)
        pos_list[i] = pos

    return pos_list


def read_pos(outcar_lines, contcar_lines, poscar_lines, species):
    """
    read species and composition from contcar and poscar file
    :param species: users-defined species sequence
    :param outcar_lines: read lines of outcar
    :param contcar_lines: read lines of contcar
    :param poscar_lines: read lines of poscar
    :return: energy, species and position info in seq from contcar with unrelaxed coords from poscar
                pos_list format [species, position, spin]
    """
    # read enrg
    enrg = None
    for i in range(len(outcar_lines)):
        if "TOTEN" in outcar_lines[i]:
            enrg = outcar_lines[i].split()
            enrg = float(enrg[4])
    # read coords and comps
    poscar_coords = read_poscar_coords(poscar_lines)
    scale = float(contcar_lines[1])
    lat_vec = []
    for i in range(2, 5):
        vec = np.array([float(x) for x in contcar_lines[i].split()]) * scale
        lat_vec.append(vec)
    spec = contcar_lines[5].split()
    comp = contcar_lines[6].split()
    atom_sum = int(np.sum(np.array(comp, dtype=np.float64)))
    spec_list = [[]] * len(species)
    pos_list = [[]] * atom_sum
    seq = []
    for i in range(len(spec)):
        for j in range(int(comp[i])):
            seq.append(spec[i])
    for i in range(atom_sum):
        pnt = contcar_lines[8 + i].split()
        pos = []
        for j in range(3):
            pos.append(float(pnt[j]))
        if 'artesian' in contcar_lines[7]:
            pos = cart_to_frac(pos, lat_vec)
        pos = apply_pbc(pos)
        dist = 0.3
        for j in range(len(poscar_coords)):
            new_dist = np.sum(np.abs(np.subtract(pos, poscar_coords[j])))
            if new_dist < dist:
                dist = new_dist
                pos_list[i] = [seq[i], poscar_coords[j], 0]  # 0 means non-mag
        if dist == 0.3:
            print('no match for certain coordinates! possible over-relaxed structure!')
    spec_dict = {i: seq.count(i) for i in seq}
    for i in range(len(species)):
        for keys in spec_dict:
            if species[i] == keys:  # if element from user-specified species is in contcar
                spec_list[i] = [spec_dict[keys], species[i]]
                break
            else:                   # if element from user-specified species is not in contcar
                spec_list[i] = [0, species[i]]

    return enrg, spec_list, pos_list


def read_mag_pos(outcar_lines, contcar_lines, poscar_lines, species):
    """
    read energy, position, and magnetism from outcar and contcar
    :param outcar_lines: read lines of outcar
    :param poscar_lines: read lines of poscar
    :param contcar_lines: read lines of contcar
    :param species: users-defined species sequence
    :return: energy, position and magnetism in seq
    """
    enrg, spec_list, mag_pos_list = read_pos(outcar_lines, contcar_lines, poscar_lines, species)
    atom_num = len(mag_pos_list)
    for i in range(len(outcar_lines)):
        if "magnetization (x)" in outcar_lines[i]:
            for j in range(atom_num):
                mag = outcar_lines[i + j + 4].split()
                mag_pos_list[j][2] = float(mag[4])

    return enrg, spec_list, mag_pos_list


def compile_vasp(root_dir, output_dir, species, read_mag, use_spin_tol, spin_tol):
    """
    generate metadata from vasp calculation
    :param root_dir: the top directory of your vasp data
    :param output_dir: path and name of your output metadata file
    :param species: species and sequence
    :param read_mag: choose whether to read magnetism
    :param use_spin_tol: choose whether to use spin tolerance
    :param spin_tol: spin tolerance here in the same sequence of species
    :return: metadata of vasp results
    """
    output = open(output_dir, 'w')
    for subdir, dirs, files in os.walk(root_dir):
        poscar_lines = []
        outcar_lines = []
        contcar_lines = []
        name = subdir.replace(root_dir, "")
        name = name.replace("/", "")
        flag = 0
        for file in files:
            if 'POSCAR' in files and 'CONTCAR' in files and 'OUTCAR' in files:
                if file == 'CONTCAR':
                    contcar = open(subdir + '/' + file, 'r')
                    contcar_lines = contcar.readlines()
                    contcar.close()
                if file == "POSCAR":
                    poscar = open(subdir + '/' + file, 'r')
                    poscar_lines = poscar.readlines()
                    poscar.close()
                    if len(poscar_lines) == 0:
                        print(subdir, ': empty structure file!')
                if file == "OUTCAR":
                    outcar = open(subdir + '/' + file, 'r')
                    outcar_lines = outcar.readlines()
                    outcar_len = len(outcar_lines)
                    for i in range(outcar_len):
                        if 'Elapsed time (sec):' in outcar_lines[i]:
                            flag = 1
                    if flag == 0:
                        print(subdir, ': unfinished job!')
                    outcar.close()
            else:
                print(subdir, ': no vasp files!')
        if len(poscar_lines) > 0 and len(outcar_lines) > 0 and len(contcar_lines) > 0 and flag == 1:
            lat_const, lat_vec, lat_ang = read_lat(poscar_lines)
            if not read_mag:
                enrg, spec_list, pos_list = read_pos(outcar_lines, contcar_lines, poscar_lines, species)
            else:
                enrg, spec_list, pos_list = read_mag_pos(outcar_lines, contcar_lines, poscar_lines, species)
                if use_spin_tol:
                    for mag in pos_list:
                        if len(spin_tol[species.index(mag[0])]) == 1:
                            if spin_tol[species.index(mag[0])] == 0:
                                mag[2] = 0
                            elif np.abs(mag[2]) < spin_tol[species.index(mag[0])]:
                                mag[2] = 0
                            else:
                                mag[2] = 1 * np.sign(mag[2])
                        elif len(spin_tol[species.index(mag[0])]) == 2:
                            if np.abs(mag[2]) < min(spin_tol[species.index(mag[0])]):
                                mag[2] = 0
                            elif np.abs(mag[2]) < max(spin_tol[species.index(mag[0])]):
                                mag[2] = 1 * np.sign(mag[2])
                            else:
                                mag[2] = 2 * np.sign(mag[2])
                        else:
                            print('Error: check the spin tolerance setting!')
            pos_list.sort(key=lambda x: species.index(x[0]))
            # start writing outputs
            # print(name)
            output.write("# ")
            for i in range(len(species)):
                output.write(str(species[i]) + " ")
            output.write('\n')
            for i in range(len(spec_list)):
                output.write(str(spec_list[i][0]) + "\t")
            output_line = name + "\t" + str(enrg) + "\t" \
                + str(lat_const[0]) + "\t" + str(lat_const[1]) + "\t" + str(lat_const[2]) + "\t" \
                + str(lat_ang[0]) + "\t" + str(lat_ang[1]) + "\t" + str(lat_ang[2]) + "\n"
            output.write(output_line)
            output_line = str(lat_vec[0][0]) + "\t" + str(lat_vec[0][1]) + "\t" + str(lat_vec[0][2]) + "\n"
            output.write(output_line)
            output_line = str(lat_vec[1][0]) + "\t" + str(lat_vec[1][1]) + "\t" + str(lat_vec[1][2]) + "\n"
            output.write(output_line)
            output_line = str(lat_vec[2][0]) + "\t" + str(lat_vec[2][1]) + "\t" + str(lat_vec[2][2]) + "\n"
            output.write(output_line)
            for i in range(len(pos_list)):
                output_line = "\t" + str(i) + "\t" + str(species.index(pos_list[i][0])) + "\t" \
                              + str(pos_list[i][2]) + "\t" + str(pos_list[i][1][0]) + "\t" \
                              + str(pos_list[i][1][1]) + "\t" + str(pos_list[i][1][2]) + "\n"
                output.write(output_line)
    output.close()



# Stand in code for ATAT input files
# def read_atat_pos(poscar_lines, species, lat_vec):
#     """
#     read species and composition from atat-formatted files
#     :param poscar_lines: read lines of atat-formatted poscar files, i.e. no species information
#     :param species: users-defined species sequence
#     :param lat_vec: lattice vector
#     :return: species and position info from atat-formatted files
#     """
#     pos_start = 0
#     for i in range(len(poscar_lines)):
#         if 'artesian' in poscar_lines[i] or 'irect' in poscar_lines[i]:
#             pos_start = i + 1
#     atom_sum = len(poscar_lines) - pos_start
#     spec_list = [[]] * len(species)
#     pos_list = [[]] * atom_sum
#     seq = []
#     for i in range(atom_sum):
#         pnt = poscar_lines[pos_start + i].split()
#         pos = []
#         for j in range(3):
#             pos.append(float(pnt[j]))
#         if 'artesian' in poscar_lines[pos_start - 1]:
#             pos = cart_to_frac(pos, lat_vec)
#         pos = apply_pbc(pos)
#         pos_list[i] = [pnt[3], pos]
#         seq.append(pnt[3])
#     pos_list.sort(key=lambda x: species.index(x[0]))
#     spec_dict = {i: seq.count(i) for i in seq}
#     for i in range(len(species)):
#         for keys in spec_dict:
#             if species[i] == keys:
#                 spec_list[i] = [spec_dict[keys], species[i]]
#                 break
#             else:
#                 spec_list[i] = [0, species[i]]
#
#     return spec_list, pos_list
#
#
# def read_pos(contcar_lines, species, lat_vec):
#     """
#     read species and composition from contcar/poscar file
#     :param species: users-defined species sequence
#     :param contcar_lines: read lines of contcar
#     :param lat_vec: lattice vector
#     :return: species and position info from contcar/poscar
#     """
#     spec = contcar_lines[5].split()
#     comp = contcar_lines[6].split()
#     atom_sum = int(np.sum(np.array(comp, dtype=np.float64)))
#     spec_list = [[]] * len(species)
#     pos_list = [[]] * atom_sum
#     seq = []
#     for i in range(len(spec)):
#         for j in range(int(comp[i])):
#             seq.append(spec[i])
#     for i in range(atom_sum):
#         pnt = contcar_lines[8 + i].split()
#         pos = []
#         for j in range(3):
#             pos.append(float(pnt[j]))
#         if 'artesian' in contcar_lines[7]:
#             pos = cart_to_frac(pos, lat_vec)
#         pos = apply_pbc(pos)
#         pos_list[i] = [seq[i], pos]
#     pos_list.sort(key=lambda x: species.index(x[0]))
#     spec_dict = {i: seq.count(i) for i in seq}
#     for i in range(len(species)):
#         for keys in spec_dict:
#             if species[i] == keys:  # if element from user-specified species is in contcar
#                 spec_list[i] = [spec_dict[keys], species[i]]
#                 break
#             else:                   # if element from user-specified species is not in contcar
#                 spec_list[i] = [0, species[i]]
#
#     return spec_list, pos_list


# def read_contcar_seq(contcar_lines):
#     """
#     read vasp output species sequence
#     :param contcar_lines: read lines of contcar
#     :return: atomic sequence in contcar and outcar for one-to-one correspondence
#     """
#     spec = contcar_lines[5].split()
#     comp = contcar_lines[6].split()
#     seq = []
#     for i in range(len(spec)):
#         for j in range(int(comp[i])):
#             seq.append(spec[i])
#
#     return seq
