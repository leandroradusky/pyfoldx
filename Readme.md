# pyFoldX: Python Bindings for FoldX

Official bindings of [FoldX](http://foldxsuite.crg.eu/) for Python programming language.

## Dependencies

- Python >= 3.6
- Linux or MacOS
- [Biopython](https://biopython.org/), [pandas](https://pandas.pydata.org/), [numpy](https://numpy.org/)

## Installing

0) Before start, a FoldX executable is needed by pyFoldX to function. You can obtain it for free upon registration [here](http://foldxsuite.crg.eu/). Once you have the  executable in your filesystem, add the following line to your .bashrc:

    export FOLDX_BINARY=/your/path/to/foldx

1) Download the git codebase or clone this repository executing:

	git clone https://github.com/leandroradusky/pyFoldX.git

2) (optional) Create a virtual environment to install the tool:

    vitualenv pyfoldx-env
    . pyfoldx-env/bin/activate

3) Install pyFoldx and its requirements:

    cd pyfoldx
    python setup.py install

## Basic tutorials

- [Working with PDB Structures](https://github.com/leandroradusky/pyFoldX/blob/master/notebooks/StructureUsage.ipynb)
- [Working with Ensembles of structures](https://github.com/leandroradusky/pyFoldX/blob/master/notebooks/EnsembleStability.ipynb)
- [Mapping mutations along Ensembles](https://github.com/leandroradusky/pyFoldX/blob/master/notebooks/EnsembleMutation.ipynb)
- TODO: [Visualizing generated molecules in Jupyter notebooks]()

## Advanced tutorials

- TODO: [Assesing mutation accuracy]()
- TODO: [Including other structural and nonstructural features to improve energy predictions]()
- [Parameterizing a molecule on the fly](https://github.com/leandroradusky/pyFoldX/blob/master/notebooks/paramX_AtomNames.ipynb)
- TODO: [Generating a dataset of minimized structures]()

## Features under development

- Molecule visualization using NGL