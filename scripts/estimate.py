import argparse
import os
import math
import pandas as pd
import numpy as np
import json

from typing import Dict, Tuple, List
from pathlib import Path


def list_to_dic(l: List[Tuple[int, int]]) -> Dict[int, int]:
    """
    Convert a list of tuples to a dictionary.
    """
    list_d = {}
    for e in l: list_d[e[0]] = e[1]
    return list_d


def compute_rhas(distro_a: Dict[int, int], distro_b: Dict[int, int], deg_err: float,
                 add_err: float = 0.005) -> Dict[int, float]:
    """
    Compute the relative Hausdorff with additive slack (RHAS) distance (pointwise) between two distributions.
    Both distributions are given as dictionaries mapping degree to count.
    """
    distros = [distro_a, distro_b]
    freq_errs = [[], []]

    for i in [0, 1]:  # compare distros[i] against distros[1-i]
        A = distros[i]
        B = distros[1 - i]

        for deg, freq in A.items():
            if deg == 0 or freq == 0: continue

            # Search degrees d' in the multiplicative window around deg (inclusive)
            lo = int(math.ceil((1.0 - deg_err) * deg))
            hi = int(math.floor((1.0 + deg_err) * deg))

            min_err = float("inf")
            for d_ in range(lo, hi + 1):
                other_freq = B.get(d_)
                if other_freq is None: continue

                diff = abs(freq - other_freq)
                residual = max(0.0, diff - add_err)
                rel_err = residual / freq
                if rel_err < min_err: min_err = rel_err

            # "cap at 100" if nothing matched in the window
            if min_err == float("inf"): min_err = 100.0

            freq_errs[i].append([deg, min_err])

    rh_max = {}
    first = list_to_dic(freq_errs[0])
    second = list_to_dic(freq_errs[1])
    max_deg = max(max(first), max(second))
    for deg in range(max_deg + 1):
        if deg in first and deg in second:
            rh_max[deg] = max(first[deg], second[deg])
        elif deg in first:
            rh_max[deg] = first[deg]
        elif deg in second:
            rh_max[deg] = second[deg]
    return rh_max


def fetch_exact(root_path: str, dataset_name: str) -> pd.DataFrame:
    """
    Fetch (ground truth) exact data from os files
    """
    print(f'Retrieving exact dataframe...')
    dataset_root = Path(root_path) / Path(dataset_name)
    # -- read exact df
    exact_path = Path(f'{dataset_root}') / Path(f'{dataset_name}_exact.txt')
    exact_df = pd.read_csv(exact_path, sep=',')
    print(f'> Fetched exact dataframe with shape: {exact_df.shape}')
    return exact_df


def fetch_data(trial_root_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Fetch data from experiments, from os files -> build head and tail dataframes (for a single trial)
    """

    # -- read params (json file)
    params_path = trial_root_path / Path('params_info.json')
    with open(params_path, 'r') as f: params = json.load(f)

    # -- read degree estimates S_h, S_t
    cur_degree_head_df = pd.read_csv(trial_root_path / Path(f'exp_head_node_degrees.txt'), sep=',',
                                     names=['node_id', 'est_degree'])
    cur_degree_tail_df = pd.read_csv(trial_root_path / Path(f'exp_tail_node_degrees.txt'), sep=',',
                                     names=['node_id', 'est_degree'])
    print(f'[Degrees] Head Shape: {cur_degree_head_df.shape} | Tail Shape: {cur_degree_tail_df.shape}')

    # -- read triangles estimates T_h, T_t
    cur_triangles_head_df = pd.read_csv(trial_root_path / Path('exp_head_node_triangles.txt'), sep=',',
                                        names=['node_id', 'est_local_triangles'])
    cur_triangles_tail_df = pd.read_csv(trial_root_path / Path('exp_tail_node_triangles.txt'), sep=',',
                                        names=['node_id', 'est_local_triangles'])
    # group by node_id and take mean (average over trials for this experiment)
    print(f'[Triangles] Head Shape: {cur_triangles_head_df.shape} | Tail Shape: {cur_triangles_tail_df.shape}')

    # -- merge degree and triangle estimates
    head_df = pd.merge(cur_degree_head_df, cur_triangles_head_df, on='node_id', how='left')
    tail_df = pd.merge(cur_degree_tail_df, cur_triangles_tail_df, on='node_id', how='left')
    # print(f'[Merged] head dataframe with shape: {head_df.shape} | tail dataframe with shape: {tail_df.shape}')
    # -- read info.csv data
    info_df = pd.read_csv(trial_root_path / Path('exp_info.csv'),
                          names=['sample_size', 'unique_sample_size', 'aux_sample_size', 'head_budget',
                                 'tail_budget', 'deg_thresh', 'time'])
    unique_nodes_head = cur_degree_head_df['node_id'].nunique()
    unique_nodes_tail = cur_degree_tail_df['node_id'].nunique()
    info_df['unique_nodes_head'] = unique_nodes_head
    info_df['unique_nodes_tail'] = unique_nodes_tail

    return head_df, tail_df, info_df


def compute_cc_hist(hist: Dict[int, int]) -> Dict[int, int]:
    """
    Compute complementary cumulative distribution from histogram
    """
    cc_hist = {}
    max_d = int(np.ceil(max(hist.keys())))
    for d in range(max_d + 1):
        if d in hist: cc_hist[d] = sum([hist[k] for k in hist if k >= d])

    return cc_hist


def bin_exact_df(df_exact: pd.DataFrame) -> pd.DataFrame:
    df_exact_grouped = df_exact.groupby('degree').agg(
        {'triangles': 'sum', 'node_id': 'count'}).reset_index().rename(columns={'node_id': 'node_count'})
    print(f'> Exact grouped shape: {df_exact_grouped.shape}')

    # compute s_lcc (S_LCC({d}): triangles / degree choose 2
    df_exact_grouped['s_lcc'] = df_exact_grouped.apply(
        lambda row: (row['triangles'] / ((row['degree'] * (row['degree'] - 1)) / 2)) if row['degree'] >= 2 else 0,
        axis=1)

    # compute wedges (W({d}): degree choose 2 * node count
    df_exact_grouped['deg_choose_2'] = df_exact_grouped['degree'].apply(lambda d: (d * (d - 1)) / 2 if d >= 2 else 1)
    df_exact_grouped['wedges'] = df_exact_grouped['deg_choose_2'] * df_exact_grouped['node_count']

    return df_exact_grouped


# -- Retrieve Node and Wedge Averaged Binned Degree-wise Clustering Coefficient
def bin_estimated_df(df_head: pd.DataFrame, df_tail: pd.DataFrame, deg_thresh: int) -> pd.DataFrame:
    # Group by degree and compute T({d}) and C({d}) for exact, head and tail dataframes
    df_head_grouped = df_head.groupby('est_degree').agg(
        {'est_local_triangles': 'sum', 'node_id': 'count'}).reset_index().rename(columns={'node_id': 'est_node_count'})
    df_tail_grouped = df_tail.groupby('est_degree').agg(
        {'est_local_triangles': 'sum', 'node_id': 'count'}).reset_index().rename(columns={'node_id': 'est_node_count'})
    print(f'Head grouped shape: {df_head_grouped.shape} | Tail grouped shape: {df_tail_grouped.shape}')

    # compute s_lcc (S_LCC({d}): triangles / degree choose 2
    df_head_grouped['s_lcc'] = df_head_grouped.apply(
        lambda row: (row['est_local_triangles'] / ((row['est_degree'] * (row['est_degree'] - 1)) / 2)) if row[
                                                                                                              'est_degree'] >= 2 else 0,
        axis=1)
    df_tail_grouped['s_lcc'] = df_tail_grouped.apply(
        lambda row: (row['est_local_triangles'] / ((row['est_degree'] * (row['est_degree'] - 1)) / 2)) if row[
                                                                                                              'est_degree'] >= 2 else 0,
        axis=1)

    # compute wedges (W({d}): degree choose 2 * node count
    df_head_grouped['deg_choose_2'] = df_head_grouped['est_degree'].apply(lambda d: (d * (d - 1)) / 2 if d >= 2 else 1)
    df_head_grouped['est_wedges'] = df_head_grouped['deg_choose_2'] * df_head_grouped['est_node_count']
    df_tail_grouped['deg_choose_2'] = df_tail_grouped['est_degree'].apply(lambda d: (d * (d - 1)) / 2 if d >= 2 else 1)
    df_tail_grouped['est_wedges'] = df_tail_grouped['deg_choose_2'] * df_tail_grouped['est_node_count']

    # combine apx df by combining head and tail estimator given degree threshold \tau
    max_degree = max(int(df_head_grouped['est_degree'].max()), int(df_tail_grouped['est_degree'].max())) + 1
    combined_scc, combined_deg, combined_node_count, combined_triangles, combined_wedges = [], [], [], [], []
    for deg in range(max_degree + 1):
        if deg <= deg_thresh:
            # use head estimator (df_head)
            row = df_head_grouped[df_head_grouped['est_degree'] == deg]
        else:
            # use tail estimator (df_tail)
            row = df_tail_grouped[df_tail_grouped['est_degree'] == deg]
        if not row.empty:
            combined_deg.append(deg)
            combined_scc.append(row['s_lcc'].values[0])
            combined_node_count.append(row['est_node_count'].values[0])
            combined_triangles.append(row['est_local_triangles'].values[0])
            combined_wedges.append(row['est_wedges'].values[0])

    df_apx_combined = pd.DataFrame(
        {'est_degree': combined_deg, 'est_s_lcc': combined_scc, 'est_node_count': combined_node_count,
         'est_local_triangles': combined_triangles, 'est_wedges': combined_wedges})

    print(f'> Combined head and tail estimator | df shape: {df_apx_combined.shape}')

    return df_apx_combined


# -- main function for generating plots
def retrieve_estimates_df(dataset_name: str, root_path: str, target_params: Dict[str, float]) -> Tuple[
    pd.DataFrame, pd.DataFrame, int, pd.DataFrame]:
    exact_df = fetch_exact(root_path, dataset_name)

    print(f'Binning exact df per single degree...')
    binned_exact_df = bin_exact_df(exact_df)

    N_TRIALS = 10
    # cumulative dfs over trials
    cum_apx_list, cum_info_list = [], []
    p_h, p_t, h_p, t_p, a_p = target_params['p_sample_head'], target_params['p_sample_tail'], target_params[
        'head_memory_perc'], target_params['tail_memory_perc'], target_params['aux_memory_perc'],

    exp_root = Path(root_path) / Path(dataset_name) / Path(f'node_sample_h{p_h}_t{p_t}') / Path(
        f'aux_b{a_p}_hp{h_p}_tp{t_p}')
    for idx_trial in range(1, N_TRIALS + 1):
        print(f'>>> Retrieving estimates for trial {idx_trial}/{N_TRIALS}...')
        trial_root = exp_root / Path(f'trial_{idx_trial:02d}')

        head_df, tail_df, info_df = fetch_data(trial_root)
        deg_thresh = info_df['deg_thresh'].values[0]

        # Algorithm Estimate: compute NDCC and WDCC exact and estimated distributions for each single degree {d}
        print(f'Binning estimated df per single degree...')
        apx_df_per_deg = bin_estimated_df(head_df, tail_df, deg_thresh)
        # add trial id
        apx_df_per_deg['trial'] = idx_trial
        # append to cumulative df
        cum_apx_list.append(apx_df_per_deg)
        cum_info_list.append(info_df)

    binned_apx_df = pd.concat(cum_apx_list, ignore_index=True)
    cum_info_df = pd.concat(cum_info_list, ignore_index=True)
    avg_deg_thresh = int(cum_info_df['deg_thresh'].mean())

    return binned_apx_df, binned_exact_df, avg_deg_thresh, cum_info_df


def bin_df(df_input: pd.DataFrame, bins: List[float], is_exact: bool) -> pd.DataFrame:
    # -- bin dataframe with pd.cut
    df_input = df_input.copy()
    col_to_cut = 'degree' if is_exact else 'est_degree'
    df_input['degree_bin'] = pd.cut(df_input[col_to_cut], bins=bins, right=False)
    if is_exact:
        cols_to_sum = ['s_lcc', 'triangles', 'wedges', 'node_count']
    else:
        cols_to_sum = ['est_s_lcc', 'est_local_triangles', 'est_wedges', 'est_node_count']
    agg_dict = {col: 'sum' for col in cols_to_sum}
    df_input_binned = df_input.groupby('degree_bin', observed=False).agg(agg_dict).reset_index()
    df_input_binned['binned_id'] = range(len(df_input_binned))
    return df_input_binned


def retrieve_distros(df_exact: pd.DataFrame, df_apx: pd.DataFrame, bin_size: float) -> Tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    # params for RHAS (fixed)
    delta, eta = 0.1, 0.005

    max_degree = df_exact['degree'].max()
    # -- compute degree intervals
    bins = [bin_size ** i for i in range(int(math.ceil(math.log(max_degree, bin_size))) + 1)]
    df_exact_binned = bin_df(df_exact, bins, is_exact=True)
    df_exact_binned['ndcc'] = np.where(df_exact_binned['node_count'] > 0,
                                       df_exact_binned['s_lcc'] / df_exact_binned['node_count'], 0)
    df_exact_binned['wdcc'] = np.where(df_exact_binned['wedges'] > 0,
                                       df_exact_binned['triangles'] / df_exact_binned['wedges'], 0)

    cc_exact_ndcc = dict(zip(df_exact_binned['binned_id'], df_exact_binned['ndcc']))
    cc_exact_wdcc = dict(zip(df_exact_binned['binned_id'], df_exact_binned['wdcc']))

    # bin over trials
    trial_range = sorted(df_apx['trial'].unique())
    apx_list, rhas_ndcc_list, rhas_wdcc_list = [], [], []
    for trial in trial_range:
        df_sub = df_apx[df_apx['trial'] == trial]
        df_apx_binned = bin_df(df_sub, bins, is_exact=False)
        # compute metric ndcc and wdcc
        df_apx_binned['est_ndcc'] = np.where(df_apx_binned['est_node_count'] > 0,
                                             df_apx_binned['est_s_lcc'] / df_apx_binned['est_node_count'], 0)
        df_apx_binned['est_wdcc'] = np.where(df_apx_binned['est_wedges'] > 0,
                                             df_apx_binned['est_local_triangles'] / df_apx_binned['est_wedges'], 0)

        # compute RHAS with custom params

        cc_apx_ndcc = dict(zip(df_apx_binned['binned_id'], df_apx_binned['est_ndcc']))
        cc_apx_wdcc = dict(zip(df_apx_binned['binned_id'], df_apx_binned['est_wdcc']))
        rhas_distance_ndcc = compute_rhas(cc_exact_ndcc, cc_apx_ndcc, delta, eta)
        rhas_distance_wdcc = compute_rhas(cc_exact_wdcc, cc_apx_wdcc, delta, eta)

        x_axis_rh = sorted(rhas_distance_ndcc.keys())
        df_rh_ndcc = pd.DataFrame({'binned_id': x_axis_rh, 'rh_distance': [rhas_distance_ndcc[d] for d in x_axis_rh]})
        df_rh_wdcc = pd.DataFrame({'binned_id': x_axis_rh, 'rh_distance': [rhas_distance_wdcc[d] for d in x_axis_rh]})

        apx_list.append(df_apx_binned)
        rhas_ndcc_list.append(df_rh_ndcc)
        rhas_wdcc_list.append(df_rh_wdcc)

    return df_exact_binned, pd.concat(apx_list, ignore_index=True), pd.concat(rhas_ndcc_list,
                                                                              ignore_index=True), pd.concat(
        rhas_wdcc_list, ignore_index=True)


def main(args: argparse.Namespace):
    with open(f'./utils/target_params.json', 'r') as f: target_params = json.load(f)[args.dataset_name]

    binned_apx_df, binned_exact_df, avg_deg_thresh, cum_info_df = retrieve_estimates_df(args.dataset_name,
                                                                                        args.root_path, target_params)
    exact_distros, est_distros, rhas_ndcc, rhas_wdcc = retrieve_distros(binned_exact_df, binned_apx_df, args.bin_size)

    # save to output
    exact_distros.to_csv(Path(args.output_folder) / Path(f'{args.dataset_name}_exact_distros.csv'), index=False)
    est_distros.to_csv(Path(args.output_folder) / Path(f'{args.dataset_name}_est_distros.csv'), index=False)
    rhas_ndcc.to_csv(Path(args.output_folder) / Path(f'{args.dataset_name}_rhas_ndcc.csv'), index=False)
    rhas_wdcc.to_csv(Path(args.output_folder) / Path(f'{args.dataset_name}_rhas_wdcc.csv'), index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Script for generating plots")
    parser.add_argument('--dataset_name', '-n', type=str, help='Name of the dataset to process')
    parser.add_argument('--root_path', '-r', type=str, help='Root path for experiment results')
    parser.add_argument('--bin-size', '-b', type=float, help='Bin size for degree intervals, i.e., D_i = [b^i, b^{i+1})')
    parser.add_argument('--output_folder', '-o', type=str,
                        help='Folder for saving NDCC and WDCC distro, and RHAS distances')
    args = parser.parse_args()
    main(args)
