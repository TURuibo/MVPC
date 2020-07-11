# MVPC
In this repository, we provide the implementation of Missing Value PC (MVPC) for both linear Gaussian and binary cases. MVPC is a framework of causal discovery in the presence of different missingness mechanisms, including missing completely at random (MCAR), missing at random (MAR), and missing not at random (MNAR). MVPC is based on the PC algorithm and contains two methods for correcting wrong results produced by missing value issues, *Permutation-based Correction (PermC)* and *Density Ratio Weighted (DRW) correction method*. More details can be found in [the paper [1]](https://arxiv.org/abs/1807.04010). The implementation is based on the R package [pcalg](https://cran.r-project.org/web/packages/pcalg/index.html).  

## Installation

Step 1: Install [R](https://www.r-project.org).

Step 2: Intall attached packages.  
	mipfp_3.2.1, numDeriv_2016.8-1, Rsolnp_1.16, cmm_0.12, DescTools_0.99.28, e1071_1.7-2, ks_1.11.6, weights_1.0, mice_3.4.0, gdata_2.18.0, Hmisc_4.2-0, ggplot2_3.1.1, Formula_1.2-3, survival_2.44-1.1, lattice_0.20-38, mvtnorm_1.0-11, pcalg_2.6-2, graph_1.60.0, BiocGenerics_0.28.0  

## Structure of Repository

The directory structure is:

* src: The source code for MVPC.  
	* MissingValuePC.R: Implementation of the  MVPC framework.  
	* CITest.R: Implementation of conditional independence tests.  
	* SyntheticDataGeneration.R: Generation of sythetic data.  
	* Evaluation.R: Evaluation metrics, such as Structural Hamming distance, recall and precision of the causal skeleton results.  
* exps: The experiments on the binary and continuous variable datasets.  
	* bins: Binary variable experiments with different sample sizes in MAR/MNAR datasets.  
	* conti: Linear Gaussain experiments.  
		* Mul_cause.R: Experiments for the case with multiple parents of missingness indicators.  
		* NoS_mar.R: Experiments with different sample sizes in MAR datasets.  
		* NoS_mnar.R: Experiments with different sample sizes in MNAR datasets.  
		* NoV: Experiments with different number of variables in MAR datasets.  
	* results: Some results in the experiments.  
* data: Used datasets for binary experiments.

* demo.R: An example of applying MVPC to the linear Gaussian case with missing values.  

## Reference
[1] Causal Discovery in the Presence of Missing Data, AISTATS 2019, Ruibo Tu\*, Cheng Zhang\*, Paul Ackermann, Karthika Mohan, Clark Glymour, Hedvig Kjellström, Kun Zhang\*
