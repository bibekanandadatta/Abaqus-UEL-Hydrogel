# Abaqus-UEL-Hydrogel

Abaqus/Standard user element subroutines for coupled chemo-mechanics of hydrogel. The documentation provide the necessary details of the coupled theory, constitutive model, and finite element implementation procedure.

> [!WARNING]
> If you are new in the realm of finite element programming, especially developing user element codes in Abaqus, simpler examples of Abaqus UEL subroutine implementation which can also be a good starting point are: [Abaqus-UEL-Elasticity](https://github.com/bibekananda-datta/Abaqus-UEL-Elasticity) and [Abaqus-UEL-Hyperelasticity](https://github.com/bibekananda-datta/Abaqus-UEL-Hyperelasticity).



## Obtaining the file

If you have `git` installed, you can clone the repository to your local machine using
```bash
git clone https://github.com/bibekananda-datta/Abaqus-UEL-Hydrogel.git
```

You can also `fork` the repository and sync as updates are deployed, test and develop your code by creating a separate branch.

Alternatively, you can download the files in a `zip` folder in this repository using the `code` drop-down menu on the top right corner. In this approach, you will not receive any bug fixes and updates.

> [!NOTE]
> Compiling the source code requires the LAPACK library from the Intel oneMKL package. See below for the details.



## Configuring Abaqus and executing the subroutine

To run user subroutines in Abaqus, you will need to install Intel Visual Studio and Intel oneAPI package and link them with Abaqus. Follow [this blog tutorial](https://www.bibekanandadatta.com/blog/2021/link-intel-and-vs-abaqus-2020/) if you have not done it before. Additionally, [see this blog post](https://www.bibekanandadatta.com/blog/2024/lapack-Intel-Fortran-Abaqus/) to learn how to link and use the LAPACK library from the Intel oneMKL package to Abaqus user subroutines.


Make sure that the user subroutine and input file are in the same directory. Using the `Abaqus command line terminal` or `cmd terminal` or `PowerShell terminal`, you can execute the following command from the directory to execute the subroutine.

```bash
abaqus interactive double analysis ask_delete=off job=<your_job_name> input=<input_file_name.inp> user=../src/uel_hydrogel.for
```
Specify the variable names (inside < >) in the above command as needed. For additional information on executing user subroutines, check the Abaqus user manual.



## Citation

If you use this repository (documentation or source code), please consider citing this from the following:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15725220.svg)](https://doi.org/10.5281/zenodo.15725220)

APA format:
```
Datta, B., & Nguyen, Thao D. (2025, June 23). A finite element model and Abaqus user element (UEL) implementation of hydrogel chemo-mechanics. Zenodo. https://doi.org/10.5281/zenodo.15725220.
```

BibTeX:
``` bibtex
@misc{dattaFiniteElementModel2025,
  author       = {Datta, Bibekananda and Nguyen, Thao D.},
  title        = {{A finite element model and Abaqus user element (UEL) implementation of hydrogel chemo-mechanics}},
  month        = jun,
  year         = 2025,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.15725220},
  url          = {https://doi.org/10.5281/zenodo.15725220}
}
```