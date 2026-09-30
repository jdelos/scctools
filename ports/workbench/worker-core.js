/* Version 1 request/response contract. Native WASM owns computation. */
function validateRequest(data) {
  if (!data || data.version !== 1 || data.architecture !== 'qfa-graph' || data.stages !== 2 || data.phases !== 2 || !Number.isInteger(data.capacitors) || data.capacitors < 2 || data.capacitors > 3 || ![0.25, 0.5].includes(data.duty)) {
    const error = new Error('Request does not match version 1 two-phase schema.');
    error.code = data && data.phases !== 2 ? 'UNSUPPORTED_PHASE_COUNT' : 'INVALID_INPUT';
    throw error;
  }
}
function submitToWasm(data, wasm) {
  try {
    const encoded = new TextEncoder().encode(JSON.stringify(data) + '\0');
    const api = wasm.exports || wasm;
    const memory = api.memory || { buffer: api.HEAPU8.buffer };
    const malloc = api.malloc || api._malloc;
    const free = api.free || api._free;
    const submit = api.scctools_submit_json || api._scctools_submit_json;
    const release = api.scctools_free || api._scctools_free;
    const ptr = malloc(encoded.length);
    new Uint8Array(memory.buffer, ptr, encoded.length).set(encoded);
    const resultPtr = submit(ptr);
    const bytes = new Uint8Array(memory.buffer);
    let end = resultPtr; while (bytes[end]) end += 1;
    const result = JSON.parse(new TextDecoder().decode(bytes.slice(resultPtr, end)));
    release(resultPtr); free(ptr); return result;
  } catch (error) { return { version: 1, type: 'error', error: { code: error.code || 'INVALID_INPUT', message: error.message } }; }
}
if (typeof module !== 'undefined') module.exports = { validateRequest, submitToWasm };
if (typeof self !== 'undefined') self.submitToWasm = submitToWasm;
