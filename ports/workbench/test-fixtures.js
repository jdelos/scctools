const fixtures = require('./fixtures.js');
const states = Object.values(fixtures.states);
console.assert(fixtures.version === 1, 'fixture version');
console.assert(states.length === 2, 'representative states');
for (const state of states) {
  console.assert(state.parameters.phases === 2, 'two-phase scope');
  console.assert(state.matrices.transition.length === state.matrices.transition[0].length, 'square transition');
}
console.log('fixture checks passed');
