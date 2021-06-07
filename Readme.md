# pyfoldx: Python Bindings for FoldX

Official bindings of [FoldX](http://foldxsuite.crg.eu/) for Python programming language.

## Dependencies

- Python >= 3.8
- Linux or MacOS
  
## Installing

1) Before start, a FoldX executable is needed by pyFoldX to function. You can obtain it for free upon registration [here](http://foldxsuite.crg.eu/). Once you have the  executable in your filesystem, add the following line to your .bashrc:

    ```
    export FOLDX_BINARY=/your/path/to/foldx
    ```

2) (optional) Create a virtual environment to install the tool:

    ```
    virtualenv pyfoldx-env
    . pyfoldx-env/bin/activate
    ```

3) Install pyFoldx and its requirements:

    ```
    pip install pyfoldx
    ```

Note: if you have several python versions within your system or environment, you should run the proper pip (i.e. pip3.8).

## Basic tutorials

- [Working with PDB Structures](https://github.com/leandroradusky/pyFoldX/blob/master/notebooks/StructureUsage.ipynb)
- [Working with Ensembles of structures](https://github.com/leandroradusky/pyFoldX/blob/master/notebooks/EnsembleStability.ipynb)
- [Mapping mutations along Ensembles](https://github.com/leandroradusky/pyFoldX/blob/master/notebooks/EnsembleMutation.ipynb)

## Advanced tutorials

- TODO: [Assesing mutation accuracy]()
- TODO: [Including other structural and nonstructural features to improve energy predictions]()
- [Parameterizing a molecule on the fly](https://github.com/leandroradusky/pyFoldX/blob/master/notebooks/paramX_AtomNames.ipynb)
- TODO: [Generating a dataset of minimized structures]()

## Features under development

- Molecule visualization using NGL