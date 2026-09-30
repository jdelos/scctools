const assert = require('node:assert/strict');
const fixtures = require('./fixtures.js');
const handleRequest = require('./worker-core.js');

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
assert.deepEqual(handleRequest({ fixture: 'numeric', stages: '2', phases: '2' }, fixtures).state, fixtures.states.numeric);
for (const request of [
  { fixture: 'missing', stages: 2, phases: 2 },
  { fixture: 'numeric', stages: 3, phases: 2 },
  { fixture: 'numeric', stages: 'nope', phases: 2 },
  { fixture: 'numeric', stages: 2, phases: 1 }
]) {
  const response = handleRequest(request, fixtures);
  assert.equal(response.type, 'error');
  assert.equal(response.error.code, 'INVALID_INPUT');
  assert.ok(response.error.message);
}
console.log('fixture and worker checks passed');
