#!/bin/sh
set -eu
# Reproducible browser build. Generated files stay under ports/wasm/ and remain ignored.
IMAGE='emscripten/emsdk@sha256:8847dad4171ebc8a53d9ae5cda86a2546ef5b2e68834c14dc1ba2b2962e125cc'
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
exec docker run --rm -v "$ROOT:/src" -w /src "$IMAGE" sh -eu -c '
  test -d ports/libs/symengine || { echo "missing ports/libs/symengine source (SymEngine 0.13 required)" >&2; exit 2; }
  test -d ports/libs/gmp-6.2.1 || { echo "missing GMP 6.2.1 source" >&2; exit 2; }
  test -d ports/libs/mpfr-4.2.1 || { echo "missing MPFR 4.2.1 source" >&2; exit 2; }
  cmake -S ports/libs/symengine -B /tmp/scctools-symengine-wasm -DCMAKE_TOOLCHAIN_FILE="$EMSDK/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake" -DBUILD_TESTS=OFF
  cmake --build /tmp/scctools-symengine-wasm --parallel
  em++ -std=c++17 -O2 -Iports/include -I/tmp/scctools-symengine-wasm -s MODULARIZE=1 -s EXPORT_ES6=1 -s EXPORTED_FUNCTIONS="['_malloc','_free','_scctools_submit_json','_scctools_free']" -s EXPORTED_RUNTIME_METHODS="['UTF8ToString']" ports/src/*.cpp -L/tmp/scctools-symengine-wasm -lsymengine -o ports/wasm/scctools.js
'
