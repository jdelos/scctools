const assert = require('node:assert/strict');
const createModule = require('../wasm/scctools.js');
const { submitToWasm } = require('../workbench/worker-core.js');

(async () => {
  const wasm = await createModule();
  const api = wasm.exports || wasm;
  const malloc = api.malloc || api._malloc;
  const free = api.free || api._free;
  const submit = api.scctools_submit_json || api._scctools_submit_json;
  const release = api.scctools_free || api._scctools_free;
  const request = JSON.stringify({ version: 1, architecture: 'qfa-graph', stages: 2, phases: 2, capacitors: 2, duty: 0.5 });
  const input = new TextEncoder().encode(request + '\0');
  const ptr = malloc(input.length);
  new Uint8Array((api.memory || wasm.HEAPU8).buffer, ptr, input.length).set(input);
  const resultPtr = submit(ptr);
  const bytes = new Uint8Array((api.memory || wasm.HEAPU8).buffer);
  let end = resultPtr; while (bytes[end]) end++;
  const result = JSON.parse(new TextDecoder().decode(bytes.slice(resultPtr, end)));
  release(resultPtr); free(ptr);
  assert.equal(result.type, 'result');
  assert.deepEqual(result.ordering, { phases: [0, 1], capacitors: 2 });
  assert.deepEqual(result.A, [
    [['0.75', '0.5', '0.25'], ['0.25', '0.5', '0.25'], ['0.25', '0', '-0.25']],
    [['0', '0', '0'], ['-0.25', '-0.5', '-0.25'], ['-0.25', '0', '0.25']]
  ]);
  assert.deepEqual(result.m, [['0.75'], ['0.5'], ['0.25']]);
  assert.deepEqual(result.metadata, { native: true, provenance: 'native-qfa' });
  assert.equal(submitToWasm({ phases: 3, duty: 0.5, version: 1, capacitors: 2, architecture: 'qfa-graph', stages: 2 }, wasm).error.code, 'UNSUPPORTED_PHASE_COUNT');
  assert.equal(submitToWasm({ ...JSON.parse(request), capacitors: 1 }, wasm).error.code, 'INVALID_INPUT');
  const singular = submitToWasm({
    cutsets: [[[0, 0, 0]], [[0, 0, 0], [0, 0, 0]]],
    duty: [0.5, 0.5], capacitors: 1,
    operation: 'solve-charge-vectors', version: 1
  }, wasm);
  assert.equal(singular.error.code, 'SINGULAR_SYSTEM');
  console.log('generated WASM boundary passed');
})().catch((error) => { console.error(error); process.exitCode = 1; });
