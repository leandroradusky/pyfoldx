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

## Tutorials

- [Working with PDB Structures](https://github.com/leandroradusky/pyFoldX/blob/master/notebooks/StructureUsage.ipynb)
- [Working with Ensembles of structures](https://github.com/leandroradusky/pyFoldX/blob/master/notebooks/EnsembleStability.ipynb)
- (temporarily down) [Working with FoldX repaired structures.]()
- [Mapping mutations along Ensembles](https://github.com/leandroradusky/pyFoldX/blob/master/notebooks/EnsembleMutation.ipynb)
- [Parameterizing a molecule on the fly](https://github.com/leandroradusky/pyFoldX/blob/master/notebooks/paramX_AtomNames.ipynb)

## Datasets & Resources

- [Energetic and structural features computed on Missense3D-DB proposed best structure.](https://github.com/leandroradusky/pyfoldx/blob/master/notebooks/data/missense3d-benchmarking_PDBeKB_foldx.csv)
- [Energetic and structural features computed on Missense3D-DB along protein ensembles.](https://github.com/leandroradusky/pyfoldx/blob/master/notebooks/data/ensembles_mutations_foldx.csv)
- [New parameterized molecules.](https://github.com/leandroradusky/pyfoldx/tree/master/molecules)

