# Characterizing Side-Channel Leakage in ILP-Based Genomic Analytics

## Abstract

Genomic data is among the most sensitive data an individual can generate [1]. Because DNA is shared across biological relatives, it also encodes information about one’s parents, children, and extended family who never consented to its analysis, and it cannot be changed or reissued the way a compromised password or credit card number can.

Exposure of this information can reveal health conditions, ancestry, and familial relationships without consent of all parties involved. Existing approaches to securing genomic computation fall into two broad categories: cryptographic techniques and trusted execution environments (TEEs).

Cryptographic techniques, such as homomorphic encryption [2] and secure multi-party computation [3, 4], offer strong guarantees but remain too slow for many real analysis pipelines. TEEs offer fast but incomplete data protection during cloud-based processing: they isolate data from the host system, but observable behaviors, such as memory-access patterns, control flow, and execution time, may still leak information.

This project investigates such leakage in genomic analysis tools that formulate their computational problems using integer linear programming (ILP). We will begin with Aldy and CITUP, two representative genomic tools that use ILP for pharmacogenomic genotyping and tumor phylogeny reconstruction, respectively.

Our initial goal is to understand how private genomic inputs affect the behavior of the underlying ILP solver, including its branch-and-bound traversal, node selection, pruning decisions, memory accesses, and runtime. We will then evaluate whether these execution differences reveal properties of the private input and identify which components would need to become data-oblivious.

Ultimately, this work aims to inform the design of side-channel-resistant execution techniques that could apply across multiple ILP-based genomic analysis workloads running inside TEEs.

## References

[1] Richard Annan, Justin Noland, Kamaria Perkins, Xiaohong Yuan, Kaushik Roy, and Letu Qingge.  
"Genomic privacy and security in the era of artificial intelligence and quantum computing."  
*Discover Computing*, 28, June 2025.

[2] Marcelo Blatt, Alexander Gusev, Yuriy Polyakov, and Shafi Goldwasser.  
"Optimized homomorphic encryption solution for secure genome-wide association studies."  
*BMC Medical Genomics*, 13(7):83, 2020.

[3] Hyunghoon Cho, David J. Wu, and Bonnie Berger.  
"Secure genome-wide association analysis using multiparty computation."  
*Nature Biotechnology*, 36(6):547-551, 2018.

[4] Karthik A. Jagadeesh, David J. Wu, Johannes A. Birgmeier, Dan Boneh, and Gill Bejerano.  
"Deriving genomic diagnoses without revealing patient genomes."  
*Science*, 357(6352):692-695, 2017.