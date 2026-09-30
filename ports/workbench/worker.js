let wasmInstance;
(async () => {
  try { const response = await fetch('../scctools.wasm'); ({ instance: wasmInstance } = await WebAssembly.instantiate(await response.arrayBuffer(), {})); }
  catch (error) { self.postMessage({ version: 1, type: 'error', error: { code: 'INVALID_INPUT', message: `WASM load failed: ${error.message}` } }); }
})();
self.onmessage = ({ data }) => {
  if (!wasmInstance) return self.postMessage({ version: 1, type: 'error', error: { code: 'INVALID_INPUT', message: 'WASM is not ready.' } });
  self.postMessage(submitToWasm(data, wasmInstance));
};
