import argparse
import os
import json

from pathlib import Path
from typing import Dict

# -- global variables
RUN_ONE_PASS_ALGO_EXECUTABLE_PATH = 'code/build/BOLIDE'
EXACT_EXECUTABLE_PATH = 'code/build/RunExact'
N_TRIALS = 10

def write_sample_to_file(sample: Dict[int, int], file_path: Path) -> None:
    """Writes the sampled nodes and their degrees to a file.

    Args:
        sample (Dict[int, int]): A dictionary mapping node IDs to their sampled degrees.
        file_path (Path): The path to the output file.
    """
    with open(file_path, 'w') as f:
        for node, deg in sample.items():
            f.write(f'{node} {deg}\n')


def run_one_pass_algo(args: argparse.Namespace) -> None:
    # -- read dataset info
    with open('utils/datasets_map.json', 'r') as f: dataset_info = json.load(f)[args.dataset_name]

    num_nodes = dataset_info['nodes']
    num_edges = dataset_info['edges']

    with open(f'./utils/target_params.json', 'r') as f: params = json.load(f)[args.dataset_name]

    p_sample_head = params['p_sample_head']
    p_sample_tail = params['p_sample_tail']
    head_memory_perc = params['head_memory_perc']
    tail_memory_perc = params['tail_memory_perc']
    aux_memory_perc = params['aux_memory_perc']

    # -- assertions
    assert (0 <= head_memory_perc <= 1 and 0 <= tail_memory_perc <= 1), "Error! hp and ht must be in [0, 1]"
    assert (0 <= p_sample_head <= 1 and 0 <= p_sample_tail <= 1), "Error! p_sample_head and p_sample_tail must be in [0, 1]"
    assert (0 <= aux_memory_perc <= 1), "Error! aux_memory_perc must be in [0, 1]"


    output_root = Path(args.output_dir) / Path(args.dataset_name)
    os.makedirs(output_root, exist_ok=True)

    """
    Run Exact Algorithm for computing ground truth
    """
    # -- Run Exact
    exact_path = output_root / Path(f'{args.dataset_name}_exact.txt')
    command = (f"{EXACT_EXECUTABLE_PATH} "
               f"{args.dataset_path} "
               f"{exact_path} "
               )
    os.system(command)


    """
    Run one pass algorithm: NodeSampler + TriangleCounter
    """
    # -- set params: num_edges_head, num_edges_tail are not known, so pick percentages of total num_edges
    head_memory_budget = int(head_memory_perc * num_edges)
    tail_memory_budget = int(tail_memory_perc * num_edges)

    sample_root = (output_root / Path(f'node_sample_h{p_sample_head}_t{p_sample_tail}'))
    random_seed_sample_seq = [int.from_bytes(os.urandom(4), byteorder="big") for _ in range(N_TRIALS)]
    # fix random_seed sequence for all aux budgets and all trials
    random_seed_tc_seq = [int.from_bytes(os.urandom(4), byteorder="big") for _ in range(N_TRIALS)]

    for trial in range(1, N_TRIALS+1):

        aux_memory_budget = int(aux_memory_perc * num_edges)

        exp_output_root = sample_root / Path(f'aux_b{aux_memory_perc}_hp{head_memory_perc}_tp{tail_memory_perc}') / Path(f'trial_{trial:02d}')
        os.makedirs(exp_output_root, exist_ok=True)

        # -- write json info within exp_output_root with params
        params_info = {
            'idx_trial': trial,
            'dataset_name': args.dataset_name,
            'dataset_path': args.dataset_path,
            'num_nodes': num_nodes,
            'num_edges': num_edges,
            'p_sample_head': p_sample_head,
            'p_sample_tail': p_sample_tail,
            'aux_memory_perc': aux_memory_perc,
            'head_memory_perc': head_memory_perc,
            'tail_memory_perc': tail_memory_perc,
            'aux_memory_budget': aux_memory_budget,
            'head_memory_budget': head_memory_budget,
            'tail_memory_budget': tail_memory_budget,
        }

        # dump to json
        with open(exp_output_root / Path('params_info.json'), 'w') as f: json.dump(params_info, f, indent=4)

        output_path_exp = exp_output_root / Path(f'exp.txt')
        print(f"\n\n>>>>> Running Trial: {trial}/{N_TRIALS} <<<<<")
        random_seed_sample = random_seed_sample_seq[trial-1]
        random_seed = random_seed_tc_seq[trial-1]
        command = (f"{RUN_ONE_PASS_ALGO_EXECUTABLE_PATH} "
                   f"{args.dataset_path} " # graph stream path
                   f"{p_sample_head} " # p_h
                   f"{p_sample_tail} " # p_t
                   f"{0.5} " # set eps to 0.5 to recover tai roughly 10 as in the paper
                   f"{output_path_exp} " # output path
                   f"{random_seed_sample} " # random seed for node sampling
                   f"{random_seed} " # random seed for edge sampling
                   f"{aux_memory_budget} " # B_{M_h}
                   f"{head_memory_budget} " # B_{M_t}
                   f"{tail_memory_budget}" # B_{A_h} + B_{A_t}
                   )
        os.system(command)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Script for running one-pass streaming algorithm.")
    parser.add_argument('--dataset_path', '-d', type=str, required=True, help='Path to the dataset file.')
    parser.add_argument('--dataset_name', '-n', type=str, required=True, help='Name of the dataset being processed.')
    parser.add_argument('--output_dir', '-o', type=str, required=True, help='Directory to save the output results.')
    args = parser.parse_args()
    run_one_pass_algo(args)