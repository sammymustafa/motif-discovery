# Motif Discovery in Biological Sequences

This repository showcases advanced computational techniques applied to motif discovery in biological sequences, with a particular focus on DNA sequences and gene expression data in breast cancer patients.

## Key Highlights 
### 1. Gibbs Sampling for Motif Discovery
* Implementation of Gibbs Sampler: A stochastic approach to discover sequence motifs in DNA sequences, utilizing the Gibbs sampling algorithm to estimate position weight matrices (PWMs).
* Test Datasets Analysis: Application of the Gibbs sampler on synthetic datasets (with identical and degenerate motifs) and real datasets (yeast transcription factor binding sites for ACE2 and MBP1).
* PWM Identification: Determining the most consistently occurring motifs by analyzing PWMs across various runs due to the stochastic nature of Gibbs sampling.

### 2. K-means Clustering in Gene Expression Data
* K-means Algorithm Development: Custom implementation of k-means clustering to categorize gene expression profiles in breast cancer patient data, aiming to distinguish different cancer subtypes.
* Analysis on Multiple Tissue Samples: Evaluating the clustering algorithm on two distinct tissue data sets, with an emphasis on understanding the behavior of the algorithm in varying data distributions.

### 3. Viterbi Algorithm Application
* Sequence Decoding: Employing the Viterbi algorithm for decoding sequences in the context of hidden Markov models, providing insights into the underlying biological processes.

### 4. Construction of Emission and Transmission Matrices
* Matrix Formulation: Developing emission and transmission matrices as part of hidden Markov models to understand the probabilistic behavior of sequence motifs and gene expression patterns.

### 5. Convolutional Neural Network (CNN) Utilization
* Deep Learning in Motif Discovery: Leveraging CNNs to capture complex patterns in biological sequences, enhancing the capability to identify significant motifs.

### 6. Expectation-Maximization Algorithm in Motif Analysis
* Algorithmic Implementation: Applying the expectation-maximization algorithm to refine the motif discovery process, particularly in the context of incomplete or noisy data.
