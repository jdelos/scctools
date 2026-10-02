# scctools native port

C++17/SymEngine port of graph-based quantified charge-flow analysis (QFA). Current supported execution scope is two-phase Dickson converters.

## Requirements

- GNU Make
- C++17 compiler (`g++` by default)
- SymEngine C++ headers and library
- GMP
- MPFR
- Python 3, only for serving the browser workbench
- Docker, for the pinned Emscripten WebAssembly build
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
/opt/Polyspace/R2021a/bin/matlab -batch "addpath(pwd); addpath('tests'); test_dickson_hybrid_boundary; test_dickson_mode0_oracle; test_fun_loop_zero_rows; test_solve_charge_vectors_qo; test_native_scope_characterization; test_generic_qfa_boundary; test_legacy_seeman_boundary; test_compatibility_substitution_order; disp('MATLAB_PASS')"
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

Build the native WebAssembly model and run its Node integration test:

```bash
make -C ports wasm-test
```

The build uses a pinned Emscripten image and checksum-verified GMP 6.2.1 and SymEngine 0.13.0 sources. Generated files and the dependency cache under `ports/wasm/` are ignored.

Run the static request and fixture checks:

```bash
node ports/workbench/test-fixtures.js
```

## Complete symbolic QFA bundle

Version 1 topology requests generate an exact symbolic model with independent duty
`D` and complementary duty `1-D`. Existing request `duty` selects a workbench
preset; it does not silently evaluate the model. Explicit `substitutions` adds an
`evaluated` bundle; the symbolic bundle remains unchanged.

```json
{"version":1,"architecture":"qfa-graph","stages":2,"phases":2,"capacitors":2,"duty":0.5,"substitutions":{"duty":[0.5,0.5],"capacitances":[1,2],"switch_resistances":[0.1,0.2,0.3,0.4],"capacitor_esr":[0.01,0.02],"frequency":[100000]}}
```

Ordered vectors follow `symbols` and one-based `ordering` metadata. Both phase
duties are required, strictly positive, summing to one. Capacitances (F) and
frequency (Hz) must be positive; switch resistance and capacitor ESR (ohms) may
be zero. All values must be finite. Invalid vectors return `INVALID_SUBSTITUTION`.
`matlab_compatible` gives alphabetical MATLAB `symvar` order, not engine-discovered
order. Matrices are row-major arrays of expression strings; evaluated entries
are numerical strings. Compare values/algebra, never expression spelling.

`A` includes supply then capacitor rows; `B` is pumped capacitor charge;
`G` and `r` are aliases for redistributed charge `A_caps-B`. `Ar` has capacitor
rows followed by active-switch rows. Columns follow output-node order.
Capacitor voltage stress is normalized to input voltage; switch voltage is
signed per phase. Switch current stress is signed charge per output, matching
MATLAB `switch_current_ratio_N`, not peak transient current.

`ZSSL` includes `1/fsw`; `ZFSL` includes switch resistance **and capacitor ESR**.
`ZESR` exposes the capacitor-only contribution. `ZSCC = sqrt(ZSSL.^2 + ZFSL.^2)`
is an elementwise analytical root-sum-square approximation, **not transient
impedance**. Loss matrices use output order on both axes.

Live MATLAB/native algebraic parity (two- and three-capacitor models):

```bash
make -C ports matlab-parity-test
```

## Frequently asked questions

### Do the parity tests run the `ports/scctools` executable?

The original parity tests do not run `ports/scctools`. The complete-bundle test
builds `/tmp/scctools-qfa-json`, invokes the real native JSON boundary from MATLAB,
and compares charge/stress fields by exact algebraic differences. Two-capacitor
losses use exact identities. R2021a crashes expanding three-capacitor loss
identities, so those use four unequal-component exact substitutions and
50-digit RSS comparisons (tolerance `1e-40`), as permitted by the parity policy.
Explicit numerical substitutions are also checked against MATLAB.

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

Complete-bundle parity additionally exercises the public JSON boundary through
`/tmp/scctools-qfa-json`. Browser integration calls that same boundary in WASM.
