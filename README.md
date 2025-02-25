# LocalP2Mass

DeepGanModel is a project for generating deep learning models using Generative Adversarial Networks (GANs).

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [License](#license)

## Introduction

DeepGanModel is designed to simplify the process of creating and training GANs. It provides a set of tools and utilities to help you build, train, and evaluate GAN models.

## Installation

To install the required dependencies, run:

```sh
pip install -r requirements.txt
```

## Usage

To run the project, execute:

```python
python run.py
```

## Project Structure 

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
requirements.txt
run.py
src/
tests/
```

- `app`: Contains the application code for Flask to serve the correct files
- `config.py`: Configuration files for the project
- `data`: Directory for storing data files
- `docs`: Documentation files
- `src`: The main source code for the backend mass computations
- `tests`: Unit tests for the project

## License 

This project is licensed under the MIT License. See the LICENSE file for details.