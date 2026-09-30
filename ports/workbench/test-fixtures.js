const assert = require('node:assert/strict');
const fixtures = require('./fixtures.js');
const { validateRequest } = require('./worker-core.js');
const fs = require('node:fs');

assert.equal(fixtures.version, 1);
assert.deepEqual(Object.keys(fixtures.states), ['symbolic', 'numeric']);
assert.deepEqual(fixtures.states.symbolic.matrices, {
  transition: [['a', 'b'], ['c', 'd']], state: ['q0', 'q1']
});
assert.deepEqual(fixtures.states.numeric.matrices, {
  transition: [[1, 0], [0, 1]], state: [1, 0]
});
for (const state of Object.values(fixtures.states)) {
  const { transition, state: vector } = state.matrices;
  assert.equal(transition.length, 2);
  assert.ok(transition.every((row) => row.length === 2));
  assert.equal(vector.length, 2);
  assert.equal(state.parameters.stages, 2);
  assert.equal(state.parameters.phases, 2);
}
const request = { version: 1, architecture: 'qfa-graph', stages: 2, phases: 2, capacitors: 2, duty: 0.5 };
assert.doesNotThrow(() => validateRequest(request));
for (const invalid of [
  { ...request, stages: 3 },
  { ...request, stages: 'nope' },
  { ...request, capacitors: 1 },
  { ...request, phases: 1 }
]) {
  assert.throws(() => validateRequest(invalid), (error) => {
    assert.equal(error.code, invalid.phases === 1 ? 'UNSUPPORTED_PHASE_COUNT' : 'INVALID_INPUT');
    return true;
  });
}
const workerSource = fs.readFileSync(require.resolve('./worker.js'), 'utf8');
assert.match(workerSource, /locateFile:\s*\(path\)\s*=>\s*`\.\.\/wasm\/\$\{path\}`/);
console.log('fixture and worker checks passed');
