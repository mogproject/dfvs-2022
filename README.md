# dfvs-2022

This repository stores the source code of my solver for the Exact track of the [PACE 2022](https://pacechallenge.org/2022/) challenge.

### DOI of Version 1

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6604875.svg)](https://doi.org/10.5281/zenodo.6604875)

### Solver Description

- PDF: TBD

### Dependencies

- [WeGotYouCovered](https://github.com/KarlsruheMIS/pace-2019)
  - The winning solver of PACE Challenge 2019 Track A.
  - The modified source code is included in this repository.
- [GNU Make](https://www.gnu.org/software/make/)
- [CMake](https://cmake.org/) version 3.5 or higher
- C++ compiler that supports the C++14 standard.

### How to build

In this directory, run the following command.

```
make build
```

And the executable file `dfvs` will be generated under the `dist` directory.

Notes: A submission file (`dist/dfvs.tgz`) can be created by `make publish`.
