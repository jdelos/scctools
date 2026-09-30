// Versioned, published states. Workbench never computes these values.
// Matrix rows and columns retain source order; state vectors retain node order.
const WORKBENCH_FIXTURES = Object.freeze({
  version: 1,
  states: Object.freeze({
    symbolic: Object.freeze({
      id: 'symbolic', label: 'Two-phase symbolic QFA', kind: 'symbolic', architecture: 'qfa-graph',
      parameters: { stages: 2, phases: 2, capacitors: 2, duty: 0.5 },
      matrices: { transition: [['a', 'b'], ['c', 'd']], state: ['q0', 'q1'] }, // rows q0,q1; columns q0,q1

      notes: ['Graph-based QFA state; separate from Seeman flow.']
    }),
    numeric: Object.freeze({
      id: 'numeric', label: 'Two-phase numeric QFA', kind: 'numeric', architecture: 'qfa-graph',
      parameters: { stages: 2, phases: 2, capacitors: 2, duty: 0.5 },
      matrices: { transition: [[1, 0], [0, 1]], state: [1, 0] }, // rows/columns q0,q1; vector q0,q1

      notes: ['Representative fixture output; no native computation.']
    })
  })
});
if (typeof self !== 'undefined') self.WORKBENCH_FIXTURES = WORKBENCH_FIXTURES;
if (typeof module !== 'undefined') module.exports = WORKBENCH_FIXTURES;
