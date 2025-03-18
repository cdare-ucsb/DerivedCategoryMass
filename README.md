# DerivedCategoryMass

DerivedCategoryMass is a project stemming from my Ph.D. work on the mass evolution of certain objects in the catagories relevant to the _B_-Model of string theory; however, this question also provides context about the relation between birational models of certain Calabi-Yau varieties obtained from considering the moduli-space of the skyscraper sheaf under wall-crossing.



## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [License](#license)

## Introduction

In physics we typically contextualize the setting in which large-scale events occur as 4D-spacetime, which is often modeled using some Lorentzian manifold. Field theories (such as QFTs, CFTs) additionally describe the data of quantum / local fields over the spacetime manifold as sections of tensor bundles (e.g. scalar fields = Higgs bundle, spinor fields = electrons, quarks), along with some symmetry data. Motivated by the Kaluza-Klein model — which predicts that at low energies a dimension unifying electromagnetism and gravity should be undetectable — we assume that the fundamental forces arise from dimensions too small to observe directly, and thus must be mathematially _compact_. This is the first geometrically significant property we need to extract before moving on.

Next, it turns out that in current models there exist notable issues caused by the fact that the Higgs boson directly couples to every non-massless particle (this ultimately leading to large loop corrections). In theory, this should cause the Higgs boson's mass to grow as larger as the Plank mass; however, experimentally the Higgs boson's mass stays far below the Plank-scale thus causing the <a href="https://en.wikipedia.org/wiki/Hierarchy_problem">Hierarchy Problem</a>. One of the currently solutions is an added "graded" symmetry on the compactified space, known as supersymmetry (SUSY). In particular, the assumption that our theory survives at least __N__=1 supersymmetry in low energies requires that the compactified dimensions admit a covariantly constant spinor, and thus has a _SU(n)_ holonomy — in other words, our compactified dimensions are Calabi-Yau.

In the gauge theory setting of this compactified space, not all gauge fields in our field configuration should be physically significant. In order to admit a vacuum state solution, there must exist a global minimum of the Yang-Mills action; for Gauge fields, the occurs if and only if the connection on the vector bundle is a Hermite-Einstein connection. By a theorem of <a href="https://en.wikipedia.org/wiki/Kobayashi%E2%80%93Hitchin_correspondence">Donaldson-Ulhenbeck-Yau</a>, this is mathematically equivalent to the bundle being slope-stable with respect to some Kähler class. In string theory, objects depending on some choice of Kähler class are often modified by world-sheet instanton corrections. In this context, the vacuaa are described as BPS vacua and correspond to Bridgeland-stable objects in the derived category of our compactified dimensions (i.e. the Calabi-Yau threefold). The problem of finding mass for BPS states is trivial since they satisfy the BPS bound _M = |Z|_ (where Z is the central charge of a Bridgeland stability condition); however, finding the mass of other particles turns out to be a difficult problem regarding how they decompose. 


Most current supersymmetric models require 3 compactified complex dimensions for anomaly cancellation to mathematically work; unfortunately, only one example of a complex compact Kahler threefold that admits geometric Bridgeland stability conditions is known [Li (2019)](https://link.springer.com/article/10.1007/s00222-019-00888-z). There are two ways to simplify the above problem: one one hand, we can look at fewer compactified dimensions than are currently expected. In particular, if we assume that there are only two complex dimensions (as opposed to three in current supersymmetric models), our problem reduces to predicting mass on a __projective K3 surface__ — since a very general projective K3 surface has Picard rank 1, we add this to the assumption in the program to make interactions between divisors easier. On the other hand, we can also look at "local models" of our compactified Calabi-Yau manifolds, which are non-compact quasi-projective Calabi-Yau's whose derived category behaves like the derived category of our original Calabi-Yau in an open neighborhhood. Under the assumption that our 3 compactified dimensions contain a projective plane, this leads to the __"Local P2"__ model. We also provide a __"Local P1"__ model which behaves as an even simpler approximation of our projective K3 surfaces.

-------------------------------------------


## Installation


### Installation of Python
If you are interested in trying out the current implementation, the first step is making sure that you have a recent installation of Python 3 on your system. 

#### Windows
 this requires obtaining the most recent [Python Installation](https://www.python.org/downloads/) from the site, opening the `.exe` file, and checking the box for "Add Python to PATH" (important!). Next, verify the installation by running either

```
python --version
```
or
```
python3 --version
```
in your command prompt. Whichever one produces a valid output will affect the first command / term you use in the `Usage` step.

#### MacOS
Simply run
```
brew install python
```
again, make sure to verify the installation by running either 
```
python --version
```
or
```
python3 --version
```
to determine how to run the executable.

#### Linux (Ubuntu/Debian)

Run
```
sudo apt update
sudo apt install python3
```
and again verify the installation.


### Clone the repository from GitHub (or Download)

The user may notice that the code is readily available on this webpage, but not on their local machine. There are effectively two options to move the code to your local machine to run the program

1. Download the director as a `.zip` file by finding the green "<> Code" button on the repository webpage and selecting "Download ZIP". Next, unzip the directory and move it to a location you can easily navigate to from the command line.

2. Clone via `git` using the command line. This requires that the `git` package be installed on the local machine — one can follow similar steps to the installation of python to install the `git` package (e.g. for MacOS, one would simply run `brew install git` and verify installation by running `git --version`). Again, move to a local directory that is easy to navigate to from the command line and run
```
cd path/of/your/choosing && git clone https://github.com/cdare-ucsb/DerivedCategoryMass.git
```

### Install Python dependencies

While there is a decent amount of standalone code present, this program also heavily utilizes other libraries to handle front-end tasks such as translation of user inputs (via `Flask`) and training of neural networks (via `PyTorch`). Thus, in order to run this program the user must also have the same dependencies installed — to do this, run
```sh
pip install -r requirements.txt
```
(for some users who already have these packages installed, it may be judicious to simply install these requirements in a virtual Python environment so as to not adjust current installations).

-------------------------------------------

## Usage

After the installation step has been completed and Flask / Plotly / etc. have been added to the current Python installation, the GUI for the program can be run from a Flask server using

```python
python run.py
```

This should open up the default browser to a webpage where the user can choose from 3 different geometric models for the compactified dimensions relevant to superstring theory — each model simplifies the standard complex projective threefold assumption in its own unique way, which is described on the relevant page. By navigating to the bottom of the page, the current utility is that the user can explore the mass asymptotics of D-branes under different charges

![Example image of usage](/app/static/images/github_README.png)


## Project Structure 

The project is currently separated into the front-end (handled by flask) in the `app` directory, the back-end (currently python) in the `src` directory, and the unit testing in the `tests/` directory. 

```
__pycache__/
.pytest_cache/
.vscode/
app/
    __init__.py
    models/
    routes.py
    static/
    templates/
config.py
data/
docs/
README.md
LICENSE
requirements.txt
run.py
src/
    ChainComplex.py
    ChernCharacter.py
    CoherentSheaf.py
    DerivedCategoryObject.py
    DistinguishedTriangle.py
    LocalP2.py
    MassPlot.py
    model.py
    ProjectiveCY.py
    SphericalTwist.py
tests/
```

- `app`: Contains the application code for Flask to serve the correct files
- `config.py`: Configuration files for the project
- `data`: Directory for storing data files
- `docs`: Documentation files
- `src`: The main source code for the backend mass computations
- `tests`: Unit tests for the project

Currently the majority of documentatinon is given for the back-end code in the `src/` directory. The documentation modus used for this project is Doxygen, and the documentation files can be easily accessed from the root directory of this project (i.e. `/DerivedCategoryMass`) by running 

```
firefox docs/html/index.html  # Linux
open docs/html/index.html      # macOS
start docs\html\index.html     # Windows
```

![Example of Doxygen documentation](/app/static/images/github_README_documentation.png)


## Future Improvements

K3 surfaces are possibly the most interesting example of the three Calabi-Yau geometries to choose from since their cohomology is 
well-enough behaved to allow one to compute successive twists of line bundles. Currently, the program only allows for one to compute
up to two twists around line bundles — theoretically we can compute as many as we want, though implementing the mass function for this will generally require a large number of cases since there are many possible permutations of the order of O(a), O(b), O(c), O(d)... . In particular, implementing a n-fold spherical twist will generally require (n+1)! cases.

If larger numbers of spherical twists were to be considered for K3 surfaces, it would also be wise to create a class for an arbitrary number of spherical twists; currently using two classes for SphericalTwist and DoubleSphericalTwist works fine, but
it would be a bit cumbersome to keep on creating a new class for each number of twists.

It is difficult mathematically to compute higher numbers of spherical twists for Local P1 and Local P2 since the derived RHom of a line bundle with a single spherical twist leads to a long exact sequence that is incredibly difficult to resolve. Ultimately, this requires homological algebra that is beyond me.

Another potential improvement would to be implementing interactions between divisors; for starters, this would allow us to consider twists for higher-rank K3 surfaces. However, the main benefit is this would allow us to consider much more interesting geometries such as local P1 x P1, and other del Pezzo surfaces.

## License 

This project is licensed under the MIT License. See the LICENSE file for details.