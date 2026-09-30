#!/bin/sh
set -eu
# Reproducible browser build. Generated files and dependency cache stay ignored.
IMAGE='emscripten/emsdk@sha256:8847dad4171ebc8a53d9ae5cda86a2546ef5b2e68834c14dc1ba2b2962e125cc'
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
exec docker run --rm -v "$ROOT:/src" -w /src "$IMAGE" sh -eu -c '
  cache=/src/ports/wasm/.cache
  downloads=$cache/downloads; sources=$cache/sources; builds=$cache/builds
  mkdir -p "$downloads" "$sources" "$builds"
  gmp_url=https://ftp.gnu.org/gnu/gmp/gmp-6.2.1.tar.xz
  gmp_sha=fd4829912cddd12f84181c3451cc752be224643e87fac497b69edddadc49b4f2
  sym_url=https://github.com/symengine/symengine/releases/download/v0.13.0/symengine-0.13.0.tar.gz
  sym_sha=f46bcf037529cd1a422369327bf360ad4c7d2b02d0f607a62a5b09c74a55bb59
  fetch() { test -s "$1" || curl -fsSL "$2" -o "$1"; printf "%s  %s\\n" "$3" "$1" | sha256sum -c -; }
  fetch "$downloads/gmp-6.2.1.tar.xz" "$gmp_url" "$gmp_sha"
  fetch "$downloads/symengine-0.13.0.tar.gz" "$sym_url" "$sym_sha"
  image_digest=sha256:8847dad4171ebc8a53d9ae5cda86a2546ef5b2e68834c14dc1ba2b2962e125cc; flags_gmp="--host=wasm32-unknown-emscripten --disable-shared --enable-static --disable-assembly --disable-doc"
  gmp_stamp="$cache/gmp-$image_digest-$gmp_sha.stamp"; gmp_prefix="$builds/gmp-$gmp_sha/install"
  if test ! -f "$gmp_stamp"; then
    rm -rf "$sources/gmp-6.2.1" "$builds/gmp-$gmp_sha"
    mkdir "$sources/gmp-6.2.1"; tar -xJf "$downloads/gmp-6.2.1.tar.xz" --strip-components=1 -C "$sources/gmp-6.2.1"
    mkdir "$builds/gmp-$gmp_sha"
    cd "$builds/gmp-$gmp_sha"; emconfigure "$sources/gmp-6.2.1/configure" $flags_gmp --prefix="$gmp_prefix"
    emmake make -j$(nproc); emmake make install
    emar t "$gmp_prefix/lib/libgmp.a" >/dev/null
    printf "%s\\n" "$image_digest $gmp_sha $flags_gmp" > "$gmp_stamp"
  fi
  sym_flags="-DINTEGER_CLASS=gmp -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DBUILD_SHARED_LIBS=OFF -DWITH_MPFR=OFF"
  sym_stamp="$cache/symengine-$image_digest-$sym_sha-$gmp_sha.stamp"; sym_build="$builds/symengine-$sym_sha-$gmp_sha"
  if test ! -f "$sym_stamp"; then
    rm -rf "$sources/symengine-0.13.0" "$sym_build"
    mkdir "$sources/symengine-0.13.0"; tar -xzf "$downloads/symengine-0.13.0.tar.gz" --strip-components=1 -C "$sources/symengine-0.13.0"
    cmake -S "$sources/symengine-0.13.0" -B "$sym_build" -DCMAKE_TOOLCHAIN_FILE="$EMSDK/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake" $sym_flags -DGMP_ROOT="$gmp_prefix" -DGMP_INCLUDE_DIR="$gmp_prefix/include" -DGMP_LIBRARY="$gmp_prefix/lib/libgmp.a"
    cmake --build "$sym_build" --parallel
    symengine_lib=$(find "$sym_build" -name libsymengine.a -print -quit); test -n "$symengine_lib"
    emar t "$symengine_lib" >/dev/null
    printf "%s\\n" "$image_digest $sym_sha $gmp_sha $sym_flags" > "$sym_stamp"
  fi
  symengine_lib=$(find "$sym_build" -name libsymengine.a -print -quit); gmp_lib="$gmp_prefix/lib/libgmp.a"
  cd /src
  mkdir -p ports/wasm
  em++ -std=c++17 -O2 -Iports/include -I"$sources/symengine-0.13.0" -I"$sym_build" -I"$gmp_prefix/include" -s DISABLE_EXCEPTION_CATCHING=0 -s MODULARIZE=1 -s EXPORT_ES6=0 -s EXPORT_NAME=ScctoolsModule -s EXPORTED_FUNCTIONS="['_malloc','_free','_scctools_submit_json','_scctools_free']" ports/src/*.cpp "$symengine_lib" "$gmp_lib" -o ports/wasm/scctools.js
  test -s ports/wasm/scctools.js
  test -s ports/wasm/scctools.wasm
'
