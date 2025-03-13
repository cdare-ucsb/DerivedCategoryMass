# DerivedCategoryMass

DerivedCategoryMass is a project stemming from my Ph.D. work on the mass evolution of certain objects in the catagories relevant to the _B_-Model of string theory; however, this question also provides context about the relation between birational models of certain Calabi-Yau varieties obtained from considering the moduli-space of the skyscraper sheaf under wall-crossing.



## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [License](#license)

## Introduction

For a general picture of applications where this is useful, in physics we typically contextualize the setting in which large-scale events occur as 4D-spacetime (typically modeled using some Lorentzian manifold). Field theories (such as QFTs, CFTs) additionally describe the data of quantum / local fields over the spacetime manifold as sections of tensor bundles (e.g. scalar fields = Higgs bundle, spinor fields = electrons, quarks), along with some symmetry data. Motivated by the Kaluza-Klein model — which predicts that at low energies a dimension unifying electromagnetism and gravity should be undetectable — we assume that the fundamental forces arise from dimensions too small to observe directly, and thus must be mathematially _compact_. This is the first geometrically significant property we need to extract before moving on.

Next, it turns out that in current models there are some notable issues since the Higgs boson directly couples to every non-massless particle, leading to large loop corrections. However, experimentally the Higgs boson's mass stays far below the Plank-scale, leading to the <a href="https://en.wikipedia.org/wiki/Hierarchy_problem">Hierarchy Problem</a>. One of the currently solutions is an added graded symmetry on the compactified space, known as supersymmetry (SUSY). It turns out that the assumption that our theory survives at least __N__=1 supersymmetry in low energies requires that the compactified dimensions admit a covariantly constant spinor, and thus has a _SU(n)_ holonomy — in other words, our compactified dimensions are Calabi-Yau.

Backtracking to the gauge theory setting, not all gauge fields in our field configuration should be physically significant. In order to admit a vacuum state solution, there must exist a global minimum of the Yang-Mills action; for Gauge fields, the occurs if and only if the connection on the vector bundle is a Hermite-Einstein connection. By a theorem of <a href="https://en.wikipedia.org/wiki/Kobayashi%E2%80%93Hitchin_correspondence">Donaldson-Ulhenbeck-Yau</a>, this is mathematically equivalent to the bundle being slope-stable with respect to some Kähler class. In string theory, objects depending on some choice of Kähler class are often modified by world-sheet instanton corrections. In this context, the vacuaa are described as BPS vacua and correspond to Bridgeland-stable objects in the derived category of our compactified dimensions (i.e. the Calabi-Yau threefold). The problem of finding mass for BPS states is trivial since they satisfy the BPS bound _M = |Z|_ (where Z is the central charge of a Bridgeland stability condition); however, finding the mass of other particles turns out to be a difficult problem regarding how they decompose. 


Most current supersymmetric models require 3 compactified complex dimensions for anomaly cancellation to mathematically work; unfortunately, only one example of a complex compact Kahler threefold that admits geometric Bridgeland stability conditions is known [Li (2019)](https://link.springer.com/article/10.1007/s00222-019-00888-z). There are two ways to simplify the above problem: one one hand, we can look at fewer compactified dimensions than are currently expected. In particular, if we assume that there are only two complex dimensions (as opposed to three in current supersymmetric models), our problem reduces to predicting mass on a __projective K3 surface__ — since a very general projective K3 surface has Picard rank 1, we add this to the assumption in the program to make interactions between divisors easier. On the other hand, we can also look at "local models" of our compactified Calabi-Yau manifolds, which are non-compact quasi-projective Calabi-Yau's whose derived category behaves like the derived category of our original Calabi-Yau in an open neighborhhood. Under the assumption that our 3 compactified dimensions contain a projective plane, this leads to the __"Local P2"__ model. We also provide a __"Local P1"__ model which behaves as an even simpler approximation of our projective K3 surfaces.

## Installation

If you are interested in trying out the current implementation, you first want to move local directory of your choosing and clone the repository via

```
cd path/of/your/choosing && git clone https://github.com/cdare-ucsb/DerivedCategoryMass.git
```

Once the files are installed locally, there are some prerequisite Python libraries that are needed. To install the required dependencies, run:

```sh
pip install -r requirements.txt
```

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
    LocalP1.py
    LocalP2.py
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

## License 

This project is licensed under the MIT License. See the LICENSE file for details.