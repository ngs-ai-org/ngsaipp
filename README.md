# NGS-AI C++ library

---------------------------------------------

ngsaipp is the NGS-AI C++ library.

## Table of content

1. [Content](#content)
2. [Dependencies](#dependencies)
3. [Compilation](#compilation)
4. [Running unittests](#running-unittests)
5. [PacBio kinetics and epigenetics](#pacbio-kinetics-and-epigenetics)
6. [Developping with ngsaipp](#developping-with-ngsaipp)
7. [Authors](#authors)


## Content

The library contains the following directories:

1. src/ contains the C++ source code.
2. include/ contains the library header files.
3. data/ contains unittests data.
4. images/ contains doc figures.


## Dependencies

ngsaipp relies on the following third party libraries:

- pthread, normally a system library on linux system.
- [boost](https://www.boost.org/) v1.78 or higher
- [pbcopper](https://github.com/PacificBiosciences/pbcopper) v2.0.0 or higher
- [pbbam](https://github.com/PacificBiosciences/pbbam) v2.0.0 (pbbam relies on pbcopper)
- [gtest](https://github.com/google/googletest) v1.11.0 or higher


We strongly recommend to install these libraries in `/usr/local/lib` and their corresponding header files in `/usr/local/include`.  If done this way, cmake compilation configuration should work out of the box.


## Compilation

The building process uses cmake. It will compile a static and a shared ngsaipp library. To compile and install them, simply type:
```
./build.sh
```

The `CMakeLists.txt` file contains 4 variable defined at its top:

```
# user defined paths
## list of directories containing libraries headers
set(INCLUDE_DIRECTORIES "/usr/local/include/")
## list of directories in which the required libraries are installed
set(LINK_DIRECTORY      "/usr/local/lib")
## path in which the libraries will be installed
set(INSTALL_LIB_DIRECTORY   "/usr/local/lib")
## path in which the library header will be installed
set(INSTALL_HEADER_DIRECTORY   "/usr/local/include")
```

- `INCLUDE_DIRECTORIES` a space separated list of directories in which the required third party library header files are located.

- `LINK_DIRECTORY` a space separated list of directories in which all the required libraries are installed.

- `INSTALL_LIB_DIRECTORY` the directory in which ngsaipp static and shared libraries will be installed. 

- `INSTALL_HEADER_DIRECTORY` the directory in which ngsaipp header files will be installed.


## Running unittests

A set of unit tests will be compiled together with the libraries. It will be located in `<ngsaipp dir>/bin`.

To run it, type:
```
cd <ngsaipp dir>/bin
./ngsaipp_unittests --gtest_color=yes
``` 

Note the `cd` statement. It is important to execute the unittests from the executable directory. Otherwise, it won't find the necessary data in `data/`. 


## PacBio kinetics and epigenetics

An important part of ngsaipp is centered on handling PacBio BAM file epignetics signal data, modelling the signal and computing predictions to determine the presence of epigentics modifications. For an in-depth dive on this topic, please read [the dedicated documention](./PacBio_kinetics.md)


## Developping with ngsaipp

Check the API documentation in [`docs/html/index.html`](docs/html/index.html)

## Authors

* **Romain Groux**
