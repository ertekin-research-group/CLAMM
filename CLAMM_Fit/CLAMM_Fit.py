import parse
import count
import symop
import cefit
from pymatgen.core import Molecule
from pymatgen.symmetry import analyzer as syman
import json
import time
import yaml


start_time = time.time()
with open('param_in', 'r') as filehandle:
    param = yaml.safe_load(filehandle)
# to-do: write functions to read in params after setting the default values
# keys = ['do_count', 'do_fit', 'use_avg_enrg']
# for key in keys:
#     for k, v in param.items():
#         if key in k:
#             print(v)
do_count = param.get('do_count')
if do_count is None:
    do_count = True  # default value
do_fit = param.get('do_fit')
if do_fit is None:
    do_fit = True  # default value
use_avg = param.get('use_avg_enrg')
if use_avg is None:
    use_avg = False  # default value
rescale_enrg = param.get('rescale_enrg')
if rescale_enrg is None:
    rescale_enrg = False  # default value
ep_comp = param.get('ep_comp')
ep_enrg = param.get('ep_enrg')
fit_lasso = param.get('fit_lasso')
if fit_lasso is None:
    fit_lasso = False  # default value
fit_ridge = param.get('fit_ridge')
if fit_ridge is None:
    fit_ridge = False  # default value
fit_eln = param.get('fit_eln')
if fit_eln is None:
    fit_eln = False  # default value
lat_in = param.get('lat_in')
data_file = param.get('data_file')
clust_in = param.get('clust_in')
species = param.get('species')

str_list = parse.parse_str(data_file)
if use_avg:
    str_list = parse.find_avg_str(str_list)
else:
    str_list = parse.find_uniq_str(str_list)
print('# of unique structures', len(str_list), flush=True)
with open('str_out', 'w') as filehandle:
    json.dump(str_list, filehandle)
clust_list = parse.parse_clust(clust_in)
sym_list = symop.find_sym(lat_in)
symeq_clust_list = []
for orig_clust in clust_list:
    orig_clust = count.scale_clust(orig_clust)  # transform to scaled Cartesian coord
    symeq_clust_list.append(symop.find_eq_clust(sym_list, orig_clust))
with open('symeq_clust_out', 'w') as filehandle:
    json.dump(symeq_clust_list, filehandle)

# write MC rule file
# count, deco_list = parse.parse_count('count_out')
# spec_pntsym_list = []
# for symeq_clust in symeq_clust_list:
#     if len(symeq_clust[0][0]) <= 2:
#         spec_pntsym_list.append([[None]] * len(symeq_clust))
#     else:
#         pntsym_list = []
#         for clust in symeq_clust:
#             coords = clust[0]
#             hypo_molec = Molecule(['H'] * len(coords), coords)
#             pntsym = syman.PointGroupAnalyzer(hypo_molec).get_symmetry_operations()
#             pntsym_list.append(pntsym)
#         spec_pntsym_list.append(pntsym_list)
# with open('eci_out_weighted', 'r') as filehandle:
#     eci_list = json.load(filehandle)
# cefit.write_eci('CLUSTERS_weighted', symeq_clust_list, deco_list, eci_list, spec_pntsym_list, species)

if do_count:
    spec_pntsym_list = []
    new_clust_list = []
    for symeq_clust in symeq_clust_list:
        new_clust_list.append(symeq_clust[0])
        if len(symeq_clust[0][0]) <= 2:
            spec_pntsym_list.append([[None]] * len(symeq_clust))
        else:
            pntsym_list = []
            for clust in symeq_clust:
                coords = clust[0]
                hypo_molec = Molecule(['H'] * len(coords), coords)
                pntsym = syman.PointGroupAnalyzer(hypo_molec).get_symmetry_operations()
                pntsym_list.append(pntsym)
            spec_pntsym_list.append(pntsym_list)
    count_list = count.count_singlelattice(symeq_clust_list, spec_pntsym_list, str_list, new_clust_list)
    with open('count_out', 'w') as filehandle:
        json.dump(count_list, filehandle)
    # count spin pairs
    # count_spin_list = count.count_spin_pair(symeq_clust_list, spec_pntsym_list, str_list, new_clust_list)
    # with open('count_spin_out', 'w') as filehandle:
    #     json.dump(count_spin_list, filehandle)

if do_fit:
    count, deco_list = parse.parse_count('count_out')
    if rescale_enrg:
        enrg = parse.parse_scaled_enrg('str_out', ep_comp, ep_enrg)
    else:
        enrg = parse.parse_enrg('str_out')
    all_eci = cefit.all_data_lasso(count, enrg)
    with open('eci_out', 'w') as filehandle:
        json.dump(all_eci, filehandle)
    # write MC rule file
    spec_pntsym_list = []
    for symeq_clust in symeq_clust_list:
        if len(symeq_clust[0][0]) <= 2:
            spec_pntsym_list.append([[None]] * len(symeq_clust))
        else:
            pntsym_list = []
            for clust in symeq_clust:
                coords = clust[0]
                hypo_molec = Molecule(['H'] * len(coords), coords)
                pntsym = syman.PointGroupAnalyzer(hypo_molec).get_symmetry_operations()
                pntsym_list.append(pntsym)
            spec_pntsym_list.append(pntsym_list)
    with open('eci_out', 'r') as filehandle:
        eci_list = json.load(filehandle)
    cefit.write_eci('CLUSTERS', symeq_clust_list, deco_list, eci_list, spec_pntsym_list, species)

    if fit_lasso:
        lasso_eci = cefit.lasso_fit(count, enrg)
        with open('eci_out_lasso', 'w') as filehandle:
            json.dump(lasso_eci.tolist(), filehandle)
        # write MC rules after fitting
        with open('symeq_clust_out', 'r') as filehandle:
            symeq_clust_list = json.load(filehandle)
        with open('eci_out_lasso', 'r') as filehandle:
            eci_list = json.load(filehandle)
        cefit.write_eci('CLUSTERS_lasso', symeq_clust_list, deco_list, eci_list, spec_pntsym_list, species)

    if fit_ridge:
        ridge_eci = cefit.ridge_fit(count, enrg)
        with open('eci_out_ridge', 'w') as filehandle:
            json.dump(ridge_eci.tolist(), filehandle)
        # write MC rules after fitting
        with open('symeq_clust_out', 'r') as filehandle:
            symeq_clust_list = json.load(filehandle)
        with open('eci_out_ridge', 'r') as filehandle:
            eci_list = json.load(filehandle)
        cefit.write_eci('CLUSTERS_ridge', symeq_clust_list, deco_list, eci_list, spec_pntsym_list, species)

    if fit_eln:
        eln_eci = cefit.eln_fit(count, enrg)
        with open('eci_out_eln', 'w') as filehandle:
            json.dump(eln_eci.tolist(), filehandle)
        # write MC rules after fitting
        with open('symeq_clust_out', 'r') as filehandle:
            symeq_clust_list = json.load(filehandle)
        with open('eci_out_eln', 'r') as filehandle:
            eci_list = json.load(filehandle)
        cefit.write_eci('CLUSTERS_eln', symeq_clust_list, deco_list, eci_list, spec_pntsym_list, species)

print("--- %s seconds ---" % (time.time() - start_time))
