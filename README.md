# scctools

## Workbench server

Serve repository root, not `ports/workbench`, so worker relative paths can load generated WASM:

```sh
python3 -m http.server 8000 --directory ports
```

Open <http://127.0.0.1:8000/workbench/>. Required paths: `/workbench/worker.js`, `/workbench/worker-core.js`, `/wasm/scctools.js`, `/wasm/scctools.wasm`.

## Native dependency

Native and WASM builds require SymEngine **0.13.0**. `ports/wasm/build.sh` pins release tarball SHA-256 and build image. Native Makefile rejects headers whose `SYMENGINE_MAJOR_VERSION`, `SYMENGINE_MINOR_VERSION`, and `SYMENGINE_PATCH_VERSION` are not 0.13.0.
