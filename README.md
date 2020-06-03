# COVID-19 Infection

We present CRISP (*C*OVID-19 *RI*sk *S*core *P*rediction), a probabilistic graphical model for COVID-19 infection spread through a population based on the SEIR model extended by (1) mutual contacts between pairs of individuals across time across various channels (e.g., Bluetooth contact traces), as well as (2) test outcomes at given times for infection, exposure and immunity tests. The micro-level model keeps track of the infection state for each individual at every point in time, ranging from susceptible, exposed, infectious to recovered. We develop a Monte Carlo EM algorithm to infer contact-channel specific infection transmission probabilities. Our algorithm uses Gibbs sampling to draw samples of the latent infection status of each individual over the entire time period of analysis, given the latent infection status of all contacts and test outcome data. We provide implementation details to fully parallelize the inference algorithms on GPUs and enable full inference at the scale of millions of contacts between thousands of individuals. To the best of our knowledge, this is a first model for infection spread based on this detailed micro-level data; most infection models are macro-level models that reason over entire populations.

## Build Instructions

This code uses CMake. If you have homebrew installed, try first

```brew install cmake```

Also, we are using PyBind submodules. In order to get these dependencies installed, run

```git submodule update --init --recursive```

Then, the following code work

```
cmake CMakeLists.txt -DCMAKE_BUILD_TYPE=Release
make clean
make
```

If you want to test the code, just run

```
python -m test_gibbs_sampler
```

## Python COVID C++ Library Build

```g++ -O3 -std=c++11 -I/Library/Developer/CommandLineTools/Library/Frameworks/Python3.framework/Versions/3.7/Headers/ -L/Library/Developer/CommandLineTools/Library/Frameworks/Python3.framework/Versions/3.7/lib/ -l python3.7 -m64 -shared -o covid.so python_bindings.cpp```
