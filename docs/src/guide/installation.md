# Installation

EKRem requires Julia 1.10 or later. The package is not registered; it is installed from its
repository with the Julia package manager:

```julia
pkg> add https://github.com/matteomasi/EKRem
```

## PhreeqcRM library

Geochemical reactions are computed by the
[PhreeqcRM](https://phreeqc-dev.github.io/phreeqcrm/) shared library (version 3.8.6 or
later), which is not distributed with the package. It is built from the
[source distribution](https://github.com/phreeqc-dev/phreeqcrm) with CMake, as described in
its documentation. Models without geochemical reactions do not require the library.

The path of the library is set once, with the project activated:

```julia
using EKRem
set_phreeqcrm_library("C:/path/to/PhreeqcRM.dll")    # libphreeqcrm.so on Linux
```

The path is saved as a preference in `LocalPreferences.toml`, next to the active
`Project.toml`, and it is used automatically in later sessions and by the tests. When the
function is called in the global environment, the setting applies to all projects.
`set_phreeqcrm_library(path; save=false)` loads a library for the current session only, and
[`phreeqcrm_library`](@ref) returns the path in use.

The library is selected in the following order of precedence:

1. the library loaded in the current session;
2. the environment variable `EKREM_PHREEQCRM`;
3. the saved preference.

The PHREEQC databases (for example `phreeqc.dat`) are part of the PhreeqcRM distribution,
usually in the folder `database` next to the folder of the library.
