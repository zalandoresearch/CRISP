# CRISP

<details><summary>Table of Contents</summary><p>

- [CRISP](#crisp)
  - [What is CRISP](#what-is-crisp)
  - [Build Instructions](#build-instructions)
    - [Tool Setup](#tool-setup)
    - [Dependencies](#dependencies)
    - [Compilation](#compilation)
    - [Test](#test)
  - [Contributing](#contributing)
  - [Citing CRISP](#citing-crisp)
  - [License](#license)
</p></details><p></p>


## What is CRISP

CRISP (*C*OVID-19 *RI*sk *S*core *P*rediction) is a probabilistic graphical model for COVID-19 infection spread through a population based on the SEIR model extended by (1) mutual contacts between pairs of individuals across time across various channels (e.g., Bluetooth contact traces), as well as (2) test outcomes at given times for infection, exposure and immunity tests. The micro-level model keeps track of the infection state for each individual at every point in time, ranging from susceptible, exposed, infectious to recovered. Our algorithm uses block-Gibbs sampling to draw samples of the latent infection status of each individual over the entire time period of analysis, given the latent infection status of all contacts and test outcome data. Our implementation enables full inference at the scale of millions of contacts between thousands of individuals. To the best of our knowledge, this is a first model for infection spread based on this detailed micro-level data; most infection models are macro-level models that reason over entire populations.

## Build Instructions

### Tool Setup
This code uses CMake. If you have homebrew installed, try first
```brew install cmake```

Otherwise, see if you can download a pre-compiled binary at https://cmake.org/download/.

### Dependencies 
We are using PyBind submodules. In order to get these dependencies installed, run

```git submodule update --init --recursive```

### Compilation
In order to compile the C++ Gibbs sampler as a Python library, run the following step

```
cd code
cmake CMakeLists.txt
make clean
make
```

### Test
If you want to test the code, just run the following command in the ``code`` folder

```
python3 test_gibbs_sampler.py
```
You can also experiment with the option ``--setup 1`` or ``--setup 2``. 

## Contributing

Thanks for your interest in contributing! There are many ways to get involved; start with our [contributor guidelines](/CONTRIBUTING.md) and then check these [open issues](https://github.com/zalandoresearch/CRISP/issues) for specific tasks.

## Citing CRISP
If you use CRISP in a scientific publication, we would appreciate references to the following paper:

**[CRISP: A Probabilistic Model for Individual-Level COVID-19 Infection Risk Estimation Based on Contact Data. Ralf Herbrich, Rajeev Rastogi, Roland Vollgraf.](https://github.com/zalandoresearch/CRISP/raw/master/crisp.pdf)**

Biblatex entry:
```latex
@online{crisp2020,
  author       = {Ralf Herbrich and Rajeev Rastogi and Roland Vollgraf},
  title        = {CRISP: A Probabilistic Model for Individual-Level COVID-19 Infection Risk Estimation Based on Contact Data},
  date         = {2020-06-05},
  year         = {2020},
  eprintclass  = {cs.LG},
  eprinttype   = {arXiv},
  eprint       = {cs.LG/2006.XXXXX},
}
```

## License

The MIT License (MIT) Copyright © [2020] Zalando SE, https://tech.zalando.com

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
