# ThermodynamicTablesFANTOM

This repository contains a [Jupyter notebooks](#jupyter-notebooks) and [datasets](#datasets) to display and to build binary files of thermodynamic tables used in FANTOM-2D.

This repo comes as a complement of the following research paper and related supplementary material:

- Theunissen, T., Huismans, R.S., Lu, G. and Riel, N. Relative continent/mid-ocean ridge elevation: a reference case for isostasy in geodynamics. Earth-Science Reviews. [https://doi.org/10.1016/j.earscirev.2022.104153](https://doi.org/10.1016/j.earscirev.2022.104153)

- SM_ESR_isostasy repository on [github](https://github.com/tth030/SM_ESR_isostasy)

***Please cite the source when using these data.***

This Repository allows:

- Displaying and downloading thermodynamic solutions including input files, raw data and grids of density, melt fraction,...

While you can run these notebooks on Binder we encourage you to clone this repo and to run jupyter-lab locally.

## Content

- [How to run the notebooks?](#how-to-run-the-notebooks)
    - [Run in the cloud (Binder)](#run-in-the-cloud-binder)
    - [Install and run locally (Conda)](#install-and-run-locally-conda)
- [Jupyter notebooks](#jupyter-notebooks)
- [Datasets](#datasets)

## How to run the notebooks?

### Run in the cloud (Binder)

You can run the notebooks in your browser without installing anything thanks to
[binder](https://mybinder.org/). Just follow the link below and it will launch 
remotely a new notebook server for you:

- [Run on binder](https://mybinder.org/v2/gh/tth030/ThermodynamicTablesFANTOM/main?labpath=start.ipynb)

### Install and run locally (Conda)

Assuming that you have `git` and [conda](https://conda.io/docs/index.html)
installed, you can install all the packages required to run the notebooks in a
new conda environment using the following commands:

```bash
$ git clone https://github.com/tth030/ThermodynamicTablesFANTOM.git
$ cd ThermodynamicTablesFANTOM
$ conda env create -f environment.yml
$ conda activate ThermodynamicTablesFANTOM
```

You also need to install a few Jupyterlab extensions with the following command
(this step won't be required anymore with Jupyterlab >= 3.x):

```bash
$ jupyter labextension install \
    @jupyter-widgets/jupyterlab-manager \
    @pyviz/jupyterlab_pyviz \
    dask-labextension \
    ipygany
```

Finally run the command below to start the Jupyterlab application. It should
open a new tab in your browser.

```bash
$ jupyter-lab start.ipynb
```

## Jupyter notebooks

- `start.ipynb`: general introduction and disclaimers.
- `thermo.ipynb`: data display and building grid from results of thermodynamic calculations.

## Datasets

### Disclaimer:

***Please cite this paper when using these data.***

Theunissen, T., Huismans, R.S., Lu, G. and Riel, N. Relative continent/mid-ocean ridge elevation: a reference case for isostasy in geodynamics. Earth-Science Reviews. [https://doi.org/10.1016/j.earscirev.        2022.104153](https://doi.org/10.1016/j.earscirev.2022.104153)


## License

The source code and data assets are under the following licenses:

    script: MIT. (LICENSE)
    data:
        thermodyn: CC-BY-4.0 (LICENSE)

For a full description, please visit the README and LICENSE files of each package in the corresponding package folders.
