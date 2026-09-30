let wasmModule;
try {
  importScripts('worker-core.js', '../wasm/scctools.js');
  wasmModule = ScctoolsModule({ locateFile: (path) => `../wasm/${path}` });
} catch (error) {
  self.postMessage({ version: 1, type: 'error', error: { code: 'INVALID_INPUT', message: `WASM load failed: ${error.message}` } });
}
self.onmessage = async ({ data }) => {
  if (!wasmModule) return self.postMessage({ version: 1, type: 'error', error: { code: 'INVALID_INPUT', message: 'WASM is not ready.' } });
  try { self.postMessage(submitToWasm(data, await wasmModule)); }
  catch (error) { self.postMessage({ version: 1, type: 'error', error: { code: 'INVALID_INPUT', message: error.message } }); }
};
