// Versioned, published states. Workbench never computes these values.
const WORKBENCH_FIXTURES = Object.freeze({
  version: 1,
  states: Object.freeze({
    symbolic: Object.freeze({
      id: 'symbolic', label: 'Two-phase symbolic QFA', kind: 'symbolic', architecture: 'qfa-graph',
      parameters: { stages: 2, phases: 2 },
      matrices: { transition: [['a', 'b'], ['c', 'd']], state: ['q0', 'q1'] },
      notes: ['Graph-based QFA state; separate from Seeman flow.']
    }),
    numeric: Object.freeze({
      id: 'numeric', label: 'Two-phase numeric QFA', kind: 'numeric', architecture: 'qfa-graph',
      parameters: { stages: 2, phases: 2 },
      matrices: { transition: [[1, 0], [0, 1]], state: [1, 0] },
      notes: ['Representative fixture output; no native computation.']
    })
  })
});
if (typeof self !== 'undefined') self.WORKBENCH_FIXTURES = WORKBENCH_FIXTURES;
if (typeof module !== 'undefined') module.exports = WORKBENCH_FIXTURES;
