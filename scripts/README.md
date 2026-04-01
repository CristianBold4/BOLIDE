# Running Scripts 

Instructions for reproducing results obtained by running BOLIDE algorithm and computing estimates of NDCC and WDCC distributions.

Requirements: Python, numpy, pandas, matplotlib, seaborn, Jupyter Notebook. 

The scripts refer to the dataset configurations in `../utils/target_params.json` folders. Please adjust parameters accordingly. 
Moreover, the scripts access `../utils/datasets_map.json` to retrieve the number $m$ of edges of the dataset, in order to set properly the memory budget parameters. 

## Running BOLIDE algorithm

Script `run_one_pass_algo.py` takes as input the dataset path, the dataset name (used for output), and the output directory. 
First, the script runs the exact algorithm and saves the node-degree-triangle map in output. 
Then, the script runs BOLIDE algorithm for *N_TRIALS* independent repetitions (*N_TRIALS*=10 by default). 
Then, for each trial, it creates a subfolder (given the target parameters) in which it saves the node-degree maps $S_h, S_t$ and the node-triangle maps $\hat{T}_h, \hat{T}_t$. 
Moreover, it saves additional statistics such as the degree threshold $\tau$, the number of unique edges, the runtime, etc. 

## Computing binned degree-wise clustering coefficients distributions

Script `estimate.py` takes as input the dataset name, the folder produced by the above script, the bin size (for producing intervals) and the output folder in which storing the final distributions. 
The script considers the exact and estimated outputs and group them by singleton degrees, i.e., it computes $S_{LCC}(\{d\}), C(\{d\}), T(\{d\})$ and $W(\{d\})$. 
Then, it creates the degree intervals $D_i = [b^i, b^{i+1})$ given the input bin size $b$. 
Finally, it computes exact distributions $\{NDCC(D_i)\}_i$ and $\{WDCC(D_i)\}_i$, and estimated distributions $\{\widehat{NDCC}(D_i)\}_i$ and $\{\widehat{WDCC}(D_i)\}_i$, together with the pointwise RHAS distance between them. 
The parameters for RHAS distance are harcoded to $\delta = 0.1$ and $\delta = 0.005$. 


## Plotting Results

The Jupyter notebook `../notebook/plot_figures.ipynb` reads from input the exact and estimated distribution, together with the RHAS pointwise distances (please, adjust input paths accordingly), and plots the figure as in the paper. 


