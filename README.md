# AFDB90v4

This repository contains all analysis code, data and metadata generated for the current submission of our manuscript "What is hidden in the darkness? Deep-learning assisted large-scale protein family curation uncovers novel protein families and folds".

## Repo organisation

The code is organised in python notebooks (for major data analysis), python scripts (for large-scale data generation and processing) and bash scripts. The Notebooks are divided into four main analysis tasks, and describe which scripts were used to generate and analyse the data (which is also provided precomputed in https://zenodo.org/record/7759303#.ZBr5BezMIeY). Outlier detec

## How to use this repo

Just download the code in the repo and the data in [Zenodo](https://zenodo.org/record/7759303#.ZBr5BezMIeY) and follow the Jupyter Notebooks.

## Dependencies

The code was written in Python 3.6.

For the *analysis of the data*, common, standard python modules were used. Extra modules required are:
- networkx
- scipy
- seaborn
- pandas
- datashader
- geometricus

For the *generation* of the data, we used the `dbuilder` package, which is part of the ProteinUniverseAtlas project and can be found in https://github.com/ProteinUniverseAtlas/dbuilder. Shape-mers were generated using the trained `ShapemerLearn` model from geometricus, for which model code can be found in https://github.com/TurtleTools/geometricus.

