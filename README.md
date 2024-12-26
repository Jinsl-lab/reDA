# reDA: differential abundance testing on scATAC-seq data using random walk with restart
 reDA is a differential abundance test framework based on random walk with restart, for scATAC-seq data. To better measure the abundance of cells under different conditions, reDA introduces a random walk with restart, which can better capture the local and global information of the shared nearest neighbor (SNN) graph, mitigate the effects of information loss and information redundancy, and thus reliably identify condition-specific cell subpopulations based on association tests.
## installation
### Python Dependencies
reDA depends on the following Python packages:
```
numpy 1.26.2
scipy 1.11.3
pandas 1.3.5
argparse 1.4.0
scanpy 1.9.5
packaging 23.2
anndata 0.10.2
rpy2 3.4.5
multianndata 0.0.4
```
