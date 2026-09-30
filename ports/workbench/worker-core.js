function handleWorkbenchRequest(data, fixtures) {
  try {
    if (data.version !== 1) throw Object.assign(new Error('Request version must be 1.'), { code: 'INVALID_INPUT' });
    if (data.architecture !== 'qfa-graph') throw Object.assign(new Error('Architecture must be qfa-graph.'), { code: 'INVALID_INPUT' });
    if (!Number.isInteger(Number(data.capacitors)) || Number(data.capacitors) < 1) throw Object.assign(new Error('Capacitors must be positive integer.'), { code: 'INVALID_INPUT' });
    if (![0.25, 0.5].includes(Number(data.duty))) throw Object.assign(new Error('Duty must be 0.25 or 0.5.'), { code: 'INVALID_INPUT' });
    const state = fixtures.states[data.fixture];
    if (!state) throw new Error(`Unknown fixture "${data.fixture}". Choose symbolic or numeric.`);
    const stages = Number(data.stages);
    const phases = Number(data.phases);
    if (!Number.isInteger(stages) || stages < 1) throw new Error('Stages must be positive integer.');
    if (stages !== state.parameters.stages) throw new Error(`Fixture "${data.fixture}" requires ${state.parameters.stages} stages.`);
    if (Number(data.capacitors) !== state.parameters.capacitors || Number(data.duty) !== state.parameters.duty) throw Object.assign(new Error('Fixture parameters do not match request.'), { code: 'INVALID_INPUT' });
    if (phases !== state.parameters.phases) {
      const error = new Error(`Fixture "${data.fixture}" requires ${state.parameters.phases} phases.`);
      error.code = 'UNSUPPORTED_PHASE_COUNT'; throw error;
    }
    return { type: 'result', version: fixtures.version, architecture: 'qfa-graph', parameters: state.parameters, A: state.matrices.state, m: state.matrices.transition, metadata: { fixture: state.id }, ordering: ['q0', 'q1'], state };
  } catch (error) {
    return { version: fixtures.version, type: 'error', error: { code: error.code || 'INVALID_INPUT', message: error.message, action: 'Choose fixture parameters shown by selected fixture.' } };
  }
}
if (typeof module !== 'undefined') module.exports = handleWorkbenchRequest;
if (typeof self !== 'undefined') self.handleWorkbenchRequest = handleWorkbenchRequest;
