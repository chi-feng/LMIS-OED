# A layered multiple importance sampling scheme for focused optimal Bayesian experimental design

https://arxiv.org/abs/1903.11187

## Build and Install
```
git clone --recursive https://github.com/chi-feng/LMIS-OED.git
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

### Running analysis code

Not compatible with python2.7. Some python packages are required to run analysis code: 
```
$ pip3 install numpy matplotlib ipython
```
