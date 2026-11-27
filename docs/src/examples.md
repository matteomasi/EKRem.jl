# Examples

Complete simulations are provided in the
[`examples`](https://github.com/matteomasi/EKRem/tree/main/examples) folder of the repository. 
They are organised into two groups: benchmarks, which validate EKRem by comparing its numerical results with analytical solutions or those obtained using other codes, and examples, which demonstrate the application of EKRem for typical electrokinetic systems.

## Benchmarks

The benchmarks (`examples/benchmarks/1D/`) compare the results of EKRem with analytical
solutions or with PHREEQC.

- [Electromigration and electroosmosis of a trace cation](examples/benchmarks/1D/1D_trace_electromigration.md)
  (`1D_trace_electromigration.jl`): transport of a trace cation by electromigration and
  electroosmosis under a uniform electric field. The concentration profiles are compared with
  the analytical solution of Ogata and Banks (1961).
- [Binary electrolyte, fixed voltage](examples/benchmarks/1D/1D_binary_fixed_potential.md)
  (`1D_binary_fixed_potential.jl`): steady-state transport of NaCl in a porous column under a
  potential difference of 0.5 V. The concentration and potential profiles are compared with
  the analytical solution.
- [Binary electrolyte, fixed current](examples/benchmarks/1D/1D_binary_fixed_current.md)
  (`1D_binary_fixed_current.jl`): the same system under a fixed current density.
- [Krabbenhøft and Krabbenhøft benchmark](examples/benchmarks/1D/1D_benchmark_krabbenhoft.md)
  (`1D_benchmark_krabbenhoft.jl`): migration of Na⁺, Cl⁻, K⁺ and OH⁻ described by the coupled
  Nernst–Planck–Poisson equations. The transient and steady-state concentration profiles are
  compared with the solution of Krabbenhøft and Krabbenhøft (2008).
- [PHREEQC Ampere benchmark](examples/benchmarks/1D/1D_benchmark_phreeqc.md)
  (`1D_benchmark_phreeqc.jl`): verification of the coupling with PhreeqcRM. An acid front
  generated at the anode promotes the desorption of nickel from the porous medium. The results
  are compared with the Nernst–Planck module of PHREEQC (Appelo, 2017).

## Other examples

- [1D setup with electrode compartments and controlled voltage](examples/misc/1D_compartments_fixed_voltage.md)
  (`examples/misc/1D_compartments_fixed_voltage.jl`): electrokinetic transport under a constant voltage, with aqueous speciation and water equilibrium computed by PhreeqcRM.

## Running the examples

The examples have their own environment, which also contains the plotting package. They are
run from the root of the repository, for example

```console
julia --project=examples examples/benchmarks/1D/1D_trace_electromigration.jl
```

or all together with `make examples`. The examples that use PhreeqcRM require the library
(see [Installation](guide/installation.md#PhreeqcRM-library)).

## References

- Appelo, C.A.J., 2017. Solute transport solved with the Nernst–Planck equation for concrete
  pores with ‘free’ water and a double layer. Cem. Concr. Res. 101, 102–113.
  [https://doi.org/10.1016/j.cemconres.2017.08.030](https://doi.org/10.1016/j.cemconres.2017.08.030)
- Krabbenhøft, K., Krabbenhøft, J., 2008. Application of the Poisson–Nernst–Planck equations
  to the migration test. Cem. Concr. Res. 38, 77–88.
  [https://doi.org/10.1016/j.cemconres.2007.08.006](https://doi.org/10.1016/j.cemconres.2007.08.006)
- Ogata, A., Banks, R.B., 1961. A solution of the differential equation of longitudinal
  dispersion in porous media. U.S. Geological Survey Professional Paper 411-A, A1–A7.
  [https://doi.org/10.3133/pp411A](https://doi.org/10.3133/pp411A)
