<a href="https://github.com/chi-feng/LMIS-OED/actions?query=workflow%3A%22C%2B%2B+build+and+tests%22">![C++ build and tests](https://github.com/chi-feng/LMIS-OED/workflows/C++%20build%20and%20tests/badge.svg)</a>

# A layered multiple importance sampling scheme for focused optimal Bayesian experimental design

https://arxiv.org/abs/1903.11187

## Build and Install
```
$ git clone --recursive https://github.com/chi-feng/LMIS-OED.git
```
Eigen sources are included as a submodule. If you cloned the repository without the submodules, you can get them with the command
```
$ git submodule update --init
```
Prerequisites to build (on Ubuntu/Debian)
```
$ sudo apt install build-essential cmake
```
Run `build.sh` to produce build artifacts in `build/`

### Example
```
./build/Driver -experiment mossbauer -dim 3 -index 0 -poi 1 -sigeps 0.4 -design -1.3,0,1.3 -N 100 -M1 100 -M2 100 -useMIS 1 -outfile test.txt
```
- `-index 0` means the parameter of interests has index 0
- `-poi 1` means only one parameter of interest. 

### Running analysis code

Uses python2.7. Some python packages are required to run analysis code: numpy, matplotlib

