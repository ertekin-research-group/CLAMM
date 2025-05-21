import symop
import json
import yaml
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LassoCV
from sklearn.linear_model import RidgeCV
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import normalize
from sklearn.utils import resample

with open('param_in', 'r') as stream:
    param = yaml.safe_load(stream)
n = int(param['sample_times'])
sample_ratio = float(param['sample_ratio'])
kfold = int(param['kfold'])
alpha_range = param['alpha_range']
l1 = param['l1_ratio']
coef_tol = float(param['convergence'])
kf = KFold(n_splits=kfold, shuffle=True, random_state=123456)
alpha_cv = np.logspace(alpha_range[0], alpha_range[1], num=100)


def all_data_norm(count, enrg):
    """
    Normalized lasso fitting to all data
    :param enrg: energy list
    :param count: count list containing clusters, decorations, and counts
    :return: list of ECIs
    """
    coef_list = []
    count_norm = normalize(count)
    model = LassoCV(alphas=alpha_cv, cv=kf, max_iter=1000000000, tol=coef_tol)
    model.fit(count_norm, enrg)
    coef_list.append(model.intercept_)
    coef_list.extend(model.coef_.tolist())
    rmse = np.sqrt(mean_squared_error(model.predict(count_norm), enrg))
    score = model.score(count_norm, enrg)
    coef_num = np.sum(model.coef_ != 0)
    print(model.alpha_, rmse, score, coef_num, flush=True)

    return coef_list


def all_data_lasso(count, enrg):
    """
    Lasso fitting to all data
    :param enrg: energy list
    :param count: count list containing clusters, decorations, and counts
    :return: list of ECIs
    """
    coef_list = []
    model = LassoCV(alphas=alpha_cv, cv=kf, max_iter=1000000000, tol=coef_tol)
    model.fit(count, enrg)
    coef_list.append(model.intercept_)
    coef_list.extend(model.coef_.tolist())
    rmse = np.sqrt(mean_squared_error(model.predict(count), enrg))
    score = model.score(count, enrg)
    cvs = model.mse_path_.mean(axis=1).min()
    coef_num = np.sum(model.coef_ != 0)
    print(model.alpha_, rmse, score, coef_num, cvs, flush=True)

    return coef_list


def all_data_loocv_ridge(count, enrg):
    """
    Ridge LOO CV fitting to all data
    :param enrg: energy list
    :param count: count list containing clusters, decorations, and counts
    :return: list of ECIs
    """
    coef_list = []
    model = RidgeCV(alphas=alpha_cv)
    model.fit(count, enrg)
    coef_list.append(model.intercept_)
    coef_list.extend(model.coef_.tolist())
    rmse = np.sqrt(mean_squared_error(model.predict(count), enrg))
    score = model.score(count, enrg)
    cvs = model.mse_path_.mean(axis=1).min()
    coef_num = np.sum(model.coef_ != 0)
    print(model.alpha_, rmse, score, coef_num, cvs, flush=True)

    return coef_list


def ridge_fit(count, enrg):
    """
    Ridge fitting to energy per atom
    :param enrg: energy list
    :param count: count list containing clusters, decorations, and counts
    :return: list of ECIs
    """
    # bootstrap ridge
    sample_size = round(sample_ratio * len(enrg))
    attr_list = [[] for _ in range(n)]
    coef_list = [[] for _ in range(n)]
    print('Ridge: alpha, rmse, score, coef_num, cvs', flush=True)
    for i in range(n):
        x, y = resample(count, enrg, n_samples=sample_size, random_state=i)
        model = RidgeCV(alphas=alpha_cv, cv=kf)
        model.fit(x, y)
        coef_list[i].append(model.intercept_)
        coef_list[i].extend(model.coef_.tolist())
        rmse = np.sqrt(mean_squared_error(model.predict(count), enrg))
        score = model.score(count, enrg)
        cvs = model.mse_path_.mean(axis=1).min()
        coef_num = np.sum(model.coef_ != 0)
        print(model.alpha_, rmse, score, coef_num, cvs, flush=True)
        attr_list[i].extend([model.alpha_, rmse, score, float(coef_num), cvs])
    coef_mean = np.mean(coef_list, axis=0)
    with open('fit_coef_ridge', 'w') as filehandle:
        json.dump(coef_list, filehandle)
    with open('fit_attr_ridge', 'w') as filehandle:
        json.dump(attr_list, filehandle)

    return coef_mean


def lasso_fit(count, enrg):
    """
    Lasso fitting to energy per atom
    :param enrg: energy list
    :param count: count list containing clusters, decorations, and counts
    :return: list of ECIs
    """
    # bootstrap lasso
    sample_size = round(sample_ratio * len(enrg))
    coef_list = [[] for _ in range(n)]
    attr_list = [[] for _ in range(n)]
    print('alpha, rmse, score, coef_num, cvs', flush=True)
    for i in range(n):
        x, y = resample(count, enrg, n_samples=sample_size, random_state=i)
        model = LassoCV(alphas=alpha_cv, cv=kf, max_iter=1000000000, tol=coef_tol)
        model.fit(x, y)
        coef_list[i].append(model.intercept_)
        coef_list[i].extend(model.coef_.tolist())
        rmse = np.sqrt(mean_squared_error(model.predict(count), enrg))
        score = model.score(count, enrg)
        cvs = model.mse_path_.mean(axis=1).min()
        coef_num = np.sum(model.coef_ != 0)
        print(model.alpha_, rmse, score, coef_num, cvs, flush=True)
        attr_list[i].extend([model.alpha_, rmse, score, float(coef_num), cvs])
    coef_mean = np.mean(coef_list, axis=0)
    print('# of lasso selected features', np.sum(coef_mean != 0), flush=True)
    with open('fit_coef_lasso', 'w') as filehandle:
        json.dump(coef_list, filehandle)
    with open('fit_attr_lasso', 'w') as filehandle:
        json.dump(attr_list, filehandle)

    return coef_mean


def eln_fit(count, enrg):
    """
    ElasticNet fitting to energy per atom
    :param enrg: energy list
    :param count: count list containing clusters, decorations, and counts
    :return: list of ECIs
    """
    # bootstrap elasticnet
    sample_size = round(sample_ratio * len(enrg))
    attr_list = [[] for _ in range(n)]
    coef_list = [[] for _ in range(n)]
    print('Eln: alpha, l1_ratio, rmse, score, coef_num', flush=True)
    for i in range(n):
        x, y = resample(count, enrg, n_samples=sample_size, random_state=i)
        model = ElasticNetCV(alphas=alpha_cv, cv=kf, max_iter=1000000000, tol=coef_tol, l1_ratio=l1)
        model.fit(x, y)
        coef_list[i].append(model.intercept_)
        coef_list[i].extend(model.coef_.tolist())
        rmse = np.sqrt(mean_squared_error(model.predict(count), enrg))
        score = model.score(count, enrg)
        cvs = model.mse_path_.mean(axis=1).min()
        coef_num = np.sum(model.coef_ != 0)
        print(model.alpha_, model.l1_ratio_, rmse, score, coef_num, cvs, flush=True)
        attr_list[i].extend([model.alpha_, model.l1_ratio_, rmse, score, float(coef_num), cvs])
    coef_mean = np.mean(coef_list, axis=0)
    print('# of eln selected features', np.sum(coef_mean != 0), flush=True)
    with open('fit_coef_eln', 'w') as filehandle:
        json.dump(coef_list, filehandle)
    with open('fit_attr_eln', 'w') as filehandle:
        json.dump(attr_list, filehandle)

    return coef_mean


def write_eci(name, symeq_clust_list, deco_list, eci_list, pntsym_list, spec_seq):
    """
    write the clusters and ecis as a rule file for the magnetic MC simulation
    :param name: name of output file
    :param symeq_clust_list: all symmetry equivalent clusters
    :param deco_list: decoration list from count list
    :param eci_list: eci list from fitting
    :param pntsym_list: all point symmetry operations for the given clusters
    :param spec_seq: species order like ['Fe', 'Ni', 'Cr']
    :return: a file like this
            #
            Type : 1
            Motif :
            0, 0, 0 : 0, 1, 0
            0, 0, 0 : 0, -1, 0
            Deco : 0, 0 : 1, 1 : 0, 1 : 2, 1
            Enrg : -0.003 : 0.02 : -0.02 : -0.01
            #
            Type : 0
            Motif :
            0, 0, 0 : 1, 0, 0 : 0, 1, 0
            0, 0, 0 : -1, 0, 0 : 0, -1, 0
            Deco : 0, 0, 0 : 1, 1, 1 : 0, 1, 1 : 1, 0, 0 : 2, 2, 1
            Enrg : -0.002 : 0.01 : -0.025 : -0.012 : 1.1
    """
    output = open(name, 'w')
    output.write('Motif : intercept \n')
    output.write('Enrg : ' + str(eci_list[0]) + '\n#')
    for i in range(len(symeq_clust_list)):
        start = 0
        for k in range(0, i):
            start = start + len(deco_list[k])
        # multiplicity = len(symeq_clust_list[i])
        deco = deco_list[i]
        spin = symeq_clust_list[i][0][2][0]
        enrg_list = []
        output.write('\nType : ' + str(spin))
        output.write('\nMotif :\n')
        for j in range(len(symeq_clust_list[i])):
            clust = symeq_clust_list[i][j]
            motif = clust[0]
            size = len(motif)
            for k in range(size-1):
                output.write(', '.join(map(str, motif[k])) + ' : ')
            output.write(', '.join(map(str, motif[size-1])) + '\n')
        # output.write('\nMultiplicity : ' + str(multiplicity))
        output.write('Deco')
        for k in range(len(deco)):
            spec_list = symop.find_eq_spec_list(deco[k], symeq_clust_list[i][0], pntsym_list[i][0], spec_seq)
            enrg = [eci_list[start + k + 1]] * len(spec_list)
            enrg_list.extend(enrg)
            for m in range(len(spec_list)):
                output.write(' : ' + str(spec_list[m]))
        output.write('\nEnrg')
        for k in range(len(enrg_list)):
            output.write(' : ' + str(enrg_list[k]))
        output.write('\n#')

    output.close()
