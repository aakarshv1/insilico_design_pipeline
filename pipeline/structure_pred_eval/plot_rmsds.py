import matplotlib.pyplot as plt
import pandas as pd
import os
import argparse
import numpy as np

def get_best_results(casp_file):
    with open(casp_file, 'r') as file:
        lines = file.readlines()

    header = lines[0].strip(' ').split()
    data_lines = [line.strip().split() for line in lines[2:-1]]
    df = pd.DataFrame(data_lines, columns=header)
    return df[['RMS_CA', 'RMSD[L]', 'RMSD[D]']].min().to_numpy()

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('results_dir', type=str,
                        help='Directory containing RMSD scores')
    
    args = parser.parse_args()
    results_dir = args.results_dir
    graph_dir = results_dir + '/' + 'rmsd_stats'
    casp_dir = "pipeline/structure_pred_eval/casp15_results"

    if not os.path.isdir(graph_dir):
        os.mkdir(graph_dir)
    
    problems = [x for x in os.listdir(results_dir) if 'problem' in x and len(pd.read_csv(results_dir + '/' + x + '/' + 'single_motif_scores.csv')['motif_ca_rmsd']) >= 48]
    num_subplots = len(problems)

    rows = int(np.ceil(np.sqrt(num_subplots)))
    cols = int(np.ceil(num_subplots / rows))

    fig, axs = plt.subplots(rows, cols, figsize=(rows*4,cols*4))
    
    stats = []

    for problem, ax in zip(problems, axs.flatten()):
        name = problem.split('=')[1]
        rmsds = pd.read_csv(results_dir + '/' + problem + '/' + 'single_motif_scores.csv')['motif_ca_rmsd'].to_numpy()
        ax.hist(rmsds)
        # plt.title(f"Alpha Carbon RMSD for Full-Seq Conditioned Structures, {name} ({len(rmsds)} Samples)")
        ax.set_title(f"{name} ({len(rmsds)} Samples)")
        ax.set_xlabel('RMSD')
        ax.set_ylabel('Counts')

        best_casp = get_best_results(casp_dir + f"/{name}.txt")

        stats.append(
            {"Name": name, 
            "Samples": len(rmsds),
            "Genie_CA_RMSD_avg": np.mean(rmsds),
            "Genie_CA_RMSD_std": np.std(rmsds),
            "Best_RMS_CA": best_casp[0],
            "Best_RMSD[L]": best_casp[1],
            "Best_RMSD[D]": best_casp[2]
            }
        )
    pd.DataFrame(stats).to_csv(graph_dir + '/' + 'casp_comparison.csv', index=False)
    fig.suptitle("Alpha Carbon RMSD for Full-Seq Conditioned Structures, CASP15")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(graph_dir + '/' + 'rmsd.png')

    
