# Single-Pass Estimation of the Clustering Coefficient Distribution in Graph Streams

Cristian Boldrin and C. Seshadhri, "Single-Pass Estimation of the Clustering Coefficient Distribution in Graph Streams", under review for VLDB 2027.

**Supplementary Material** is located in the root folder (open pdf [here](https://github.com/CristianBold4/BOLIDE/blob/main/supp_material.pdf))

---

## Installation and Usage
Here are the instructions for running *BOLIDE* algorithm: *B*inned degree *O*ne-pass c*L*ustering coeffic*I*ents *D*istribution *E*stimation.
Code is deployed in *C++ 22* under *gcc 13.3.0* compiler. Additionally, *CMake 3.16.3+* is required.

Within the `code` folder:
1. Compile the code
   <br><br>
   `bash compile.sh`
   <br><br>
   The binaries will be generated inside the `code/build` folder.
   <br><br>

2. Preprocess the raw dataset
   <br><br>
   `./build/DataPreprocessing <dataset_path> <delimiter> <skip> <output_path>`
   <br><br>
   where *dataset_path* is the filepath to the dataset to be preprocessed,
   *delimiter* is the character used to separate the rows in the dataset,
   *skip* is the number of lines to skip before starting to read the dataset, and *output_path* is the
   path where the preprocessed dataset will be saved.
   <br><br>

3. Run the Exact Algorithm (if needed):
   <br><br>
   `./build/RunExact <graph_path> <output_path>`
   <br><br>
   where *graph_path* is the path to the preprocessed dataset at point (2) and *output_path* is the path where the output will be saved.
   <br><br>

4. Run *BOLIDE* Algorithm:
   <br><br>
   `./build/BOLIDE  <graph_path> <p_head> <p_tail> <eps> <output_path> <random_seed_sample> <random_seed_triangle_counter> <aux_sample_size> <head_budget> <tail_budget>`
   <br><br>
   where
   -  *graph_path* is the path to the preprocessed dataset at point (2);
   -   *p_head* is the probability of sampling head nodes ($p_h$ in the paper);
   -   *p_tail* is the probability for sampling tail nodes ($p_t$ in the paper) via sample-and-hold;
   -   *eps* is the parameter inducing degree threshold $\tau$, computed as $\tau = 3 \log( 1 / \varepsilon) / \varepsilon^2$ (see *headtail* paper). In our experiments, we set *eps* = 0.48 so that $\tau \approx 10$;
   -   *output_path* is the path where the output will be saved ($S_h, S_t, \hat{T}_h, \hat{T}_t$ and information containing runtime, $\tau$, etc;
   -   *random_seed_sample* is the random seed used for sampling nodes;
   -   *random_seed_triangle_counter* is the random seed used for sampling edges;
   -   *aux_sample_size* is the size of the auxiliary samples, i.e., $B_{A_h} + B_{A_t}$ assuming that $B_{A_h}| = |B_{A_t}$;
   -   *head_budget* is the memory budget for head main sample, i.e,. $B_{M_h}$;
   -   *tail_budget* is the memory budget for tail main sample, i.e,. $B_{M_t}$.
   <br><br>

## Datasets

Here are the links to the datasets we used to perform the experiments.

| Dataset                                                              | Nodes | Edges | Triangles |
|----------------------------------------------------------------------|-------|-------|-----------|
| [as-skitter](https://snap.stanford.edu/data/as-Skitter.html)         | 1.7M  | 11.1M | 28.7M     |
| [livejournal](https://snap.stanford.edu/data/soc-LiveJournal1.html)  | 4.8M  | 42.8M | 285.7M    |
| [com-orkut](https://snap.stanford.edu/data/com-Orkut.html)           | 3.1M  | 117M  | 627M      |
| [twitter](https://anlab-kaist.github.io/traces/WWW2010)              | 41M   | 1.2B  | 34.8B     |
| [com-friendster](https://snap.stanford.edu/data/com-Friendster.html) | 65M   | 1.8B  | 4.1B      |




## Experiments and Results

To reproduce experiments and results of the paper, please refer to the `scripts` folder.



