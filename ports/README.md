# scctools native port

C++17/SymEngine port of graph-based quantified charge-flow analysis (QFA). Current supported execution scope is two-phase Dickson converters.

## Requirements

- GNU Make
- C++17 compiler (`g++` by default)
- SymEngine C++ headers and library
- GMP
- MPFR
- Python 3, only for serving the browser workbench
- MATLAB R2021a with Symbolic Math Toolbox, only for MATLAB parity tests

SymEngine headers and library must come from the same installation/version. Defaults expect both under `/usr/local`:

```text
/usr/local/include/symengine
/usr/local/lib/libsymengine.*
```

Override those paths when needed:

```bash
make SYMENGINE_INCLUDE=/path/to/include SYMENGINE_LIB=/path/to/lib
```

## Build native executable

From repository root:

```bash
make -C ports
./ports/scctools
```

Clean build outputs:

```bash
make -C ports clean
```

## Run native tests

```bash
make -C ports \
  native-test \
  native-phase-test \
  native-charge-test \
  native-graph-primitives-test \
  native-matrix-test \
  native-boundary-test
```

Tests use exact symbolic comparisons for canonical two-phase Dickson cases.

## Run MATLAB parity tests

From repository root, replace MATLAB path if installed elsewhere:

```bash
/opt/Polyspace/R2021a/bin/matlab -batch "addpath(pwd); addpath('tests'); test_dickson_hybrid_boundary; test_dickson_mode0_oracle; test_fun_loop_zero_rows; test_solve_charge_vectors_qo; test_native_scope_characterization; disp('MATLAB_PASS')"
```

Expected final output:

```text
MATLAB_PASS
```

## Review browser workbench

Workbench is static and needs HTTP because it uses a Web Worker. From repository root:

```bash
python3 -m http.server 4173 --bind 127.0.0.1 --directory ports/workbench
```

Open:

```text
http://127.0.0.1:4173/
```

Workbench currently renders checked-in fixtures. It does not run native QFA in browser yet.

Run fixture validation:

```bash
node ports/workbench/test-fixtures.js
```

## Frequently asked questions

### Do the parity tests run the `ports/scctools` executable?

No. The current parity tests do not run the `ports/scctools` executable.

The native test targets compile separate test programs. These programs link directly to the production C++ source files. They call functions such as `dickson_hybrid_topology` and `solve_charge_vectors`.

For example, `native-phase-test` builds and runs `/tmp/test_phase`:

```text
ports/test/test_phase.cpp + ports/src/*.cpp -> /tmp/test_phase
```

The parity check has two parts:

1. MATLAB R2021a generates and checks the reference data in `tests/fixtures/dickson_mode0_matlab_r2021a.json`.
2. The native C++ tests check the same matrix shapes, matrix entries, order, duties, symbols, charge vectors (`A`), and conversion ratios (`m`).

The native tests contain the expected MATLAB values. They do not read the JSON fixture at run time.

The `ports/scctools` executable is currently a smoke-test program. It uses one fixed five-capacitor case. It prints only this message:

```text
Topology calculated for 5 capacitors.
```

Thus, the current tests verify parity at the C++ function boundary. They do not verify parity through the final executable interface.

A future integration must add a machine-readable command-line interface. An end-to-end test can then run `ports/scctools`, read its output, and compare that output with the MATLAB reference data. The browser can use the same interface after the native model is connected.
