import os
import json
import copy
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt


def frac_to_cart(frac_coord, basis):
    """
    Transform from fraction/direct coord to Cartesian coord for a given structure
    :param frac_coord: list in the format of [0,0,0]
    :param basis: list in the format of [[0.0, 0.0, 3.6], [0.0, 3.6, 0.0], [3.6, 0.0, 0.0]]
    :return: cart_coord
    """
    frac_coord = np.array(frac_coord)
    trans_matr = np.vstack(basis).T
    # inv_matr = np.linalg.inv(trans_matr)
    cart_coord = np.matmul(trans_matr, frac_coord.T).T

    return list(cart_coord)


def cart_to_frac(cart_coord, basis):
    """
    transform from Cartesian coordinate to direct/fraction coordinate
    :param cart_coord: list in the format of [0, 0, 0]
    :param basis: list in the format of [[0.0, 0.0, 3.6], [0.0, 3.6, 0.0], [3.6, 0.0, 0.0]]
    :return: frac_coord
    """
    cart_coord = np.array(cart_coord)
    trans_matr = np.vstack(basis).T
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


def read_str(contcar_lines, outcar_lines):
    """
    read species and composition from contcar/poscar file
    :param outcar_lines: read ines of outcar
    :param contcar_lines: read lines of contcar
    :return: species and position info in terms of a dictionary from contcar/poscar
    """
    str_dict = {}
    # read lat_vec
    scale = float(contcar_lines[1])
    lat_vec = []
    for i in range(2, 5):
        vec = np.array([float(x) for x in contcar_lines[i].split()]) * scale
        lat_vec.append(vec)
    # read lat_pnt
    spec = contcar_lines[5].split()
    comp = contcar_lines[6].split()
    seq = []
    for i in range(len(spec)):
        for j in range(int(comp[i])):
            seq.append(spec[i])
    atom_sum = int(np.sum(np.array(comp, dtype=np.float64)))
    pos_list = []
    for i in range(atom_sum):
        pnt = contcar_lines[8 + i].split()
        pos = []
        for j in range(3):
            pos.append(float(pnt[j]))
        if 'artesian' in contcar_lines[7]:
            pos = cart_to_frac(pos, lat_vec)
        pos = apply_pbc(pos)
        pos_list.append(pos)
    # read mag
    mag_list = [0.0] * atom_sum
    for i in range(len(outcar_lines)):
        if "magnetization (x)" in outcar_lines[i]:
            for j in range(atom_sum):
                mag = outcar_lines[i + j + 4].split()
                mag_list[j] = float(mag[4])
    str_dict['Spec'] = seq
    str_dict['LatVec'] = lat_vec
    str_dict['LatPnt'] = pos_list
    str_dict['Spin'] = mag_list

    return str_dict


def read_pair(str_dict, target, pair, clust):
    """
    find the number of certain pairs in a structure
    :param str_dict: read str from contcar and outcar
    :param target: species of interest
    :param pair: species to be paired with target
    :param clust: in cart coord
    :return: dict of # pairs and corresponding mag on target atoms
    """
    mag_list = []
    new_clust = copy.deepcopy(clust)
    if target in str_dict['Spec'] and pair in str_dict['Spec']:
        for i in range(len(str_dict['Spec'])):
            if str_dict['Spec'][i] == target:
                count = 0
                for j in range(len(new_clust)):
                    new_clust[j] = cart_to_frac(clust[j], str_dict['LatVec'])
                    coord = copy.deepcopy(str_dict['LatPnt'][i])
                    coord = np.sum([new_clust[j], coord], axis=0)
                    coord = apply_pbc(coord)
                    spec = 'unknown'
                    dist = 0.8
                    for k in range(len(str_dict['LatPnt'])):
                        new_dist = np.sum(np.abs(np.subtract(coord, str_dict['LatPnt'][k])))
                        if new_dist < dist:
                            dist = new_dist
                            spec = str_dict['Spec'][k]
                    if spec == pair:
                        count += 1
                    elif spec == 'unknown':
                        print('no match!')
                mag_list.append([count, str_dict['Spin'][i]])

    return mag_list


def find_mag_distrib(root_dir, output_dir, species, hist_bins=20):
    """
    get magnetic distribution of each species from the data set
    :param root_dir: the top directory of your vasp data
    :param output_dir: path and name of your output metadata file
    :param species: species and sequence
    :return: magnetic distribution list
    """
    distrib_list = [[] for _ in range(len(species))]  # distribution list for each element
    for subdir, dirs, files in os.walk(root_dir):
        contcar_lines = []
        outcar_lines = []
        seq_lines = []
        flag = 0
        for file in files:
            if 'POSCAR' in files and 'CONTCAR' in files and 'OUTCAR' in files:
                if file == 'CONTCAR':
                    seq_file = open(subdir + '/' + file, 'r')
                    seq_lines = seq_file.readlines()
                    seq_file.close()
                if file == "POSCAR":
                    contcar = open(subdir + '/' + file, 'r')
                    contcar_lines = contcar.readlines()
                    contcar.close()
                    if len(contcar_lines) == 0:
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
        if len(contcar_lines) > 0 and len(outcar_lines) > 0 and len(seq_lines) > 0 and flag == 1:
            str_dict = read_str(seq_lines, outcar_lines)
            for i in range(len(species)):
                for j in range(len(str_dict['Spec'])):
                    if str_dict['Spec'][j] == species[i]:
                        distrib_list[i].append(str_dict['Spin'][j])

    for i in range(len(distrib_list)):
        dist_i = np.array(distrib_list[i])
        # plot histogram
        plt.hist(dist_i, bins=hist_bins, color='blue', edgecolor='black', alpha=0.7)
        plt.title('Magnetic Distribution of ' + species[i])
        plt.xlabel('Magnetic Moment')
        plt.ylabel('Frequency')
        plt.savefig('mag_distrib_' + species[i] + '.png')
        plt.close()


    with open(output_dir, 'w') as filehandler:
        json.dump(distrib_list, filehandler)


def find_pair_mag(root_dir, output_dir, target, pair, clust):
    mag_dict = {}
    for subdir, dirs, files in os.walk(root_dir):
        contcar_lines = []
        outcar_lines = []
        seq_lines = []
        flag = 0
        for file in files:
            if 'POSCAR' in files and 'CONTCAR' in files and 'OUTCAR' in files:
                if file == 'CONTCAR':
                    seq_file = open(subdir + '/' + file, 'r')
                    seq_lines = seq_file.readlines()
                    seq_file.close()
                if file == "POSCAR":
                    contcar = open(subdir + '/' + file, 'r')
                    contcar_lines = contcar.readlines()
                    contcar.close()
                    if len(contcar_lines) == 0:
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
        if len(contcar_lines) > 0 and len(outcar_lines) > 0 and len(seq_lines) > 0 and flag == 1:
            str_dict = read_str(seq_lines, outcar_lines)
            mag_list = read_pair(str_dict, target, pair, clust)
            for i in range(len(mag_list)):
                count = mag_list[i][0]
                mag = mag_list[i][1]
                if count not in mag_dict.keys():
                    mag_dict[count] = [mag]
                else:
                    mag_dict[count].append(mag)

    with open(output_dir, 'w') as filehandler:
        json.dump(mag_dict, filehandler)
