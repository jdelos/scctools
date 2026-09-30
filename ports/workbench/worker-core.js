function handleWorkbenchRequest(data, fixtures) {
  try {
    const state = fixtures.states[data.fixture];
    if (!state) throw new Error(`Unknown fixture "${data.fixture}". Choose symbolic or numeric.`);
    const stages = Number(data.stages);
    const phases = Number(data.phases);
    if (!Number.isInteger(stages) || stages < 1) throw new Error('Stages must be positive integer.');
    if (stages !== state.parameters.stages) throw new Error(`Fixture "${data.fixture}" requires ${state.parameters.stages} stages.`);
    if (phases !== state.parameters.phases) throw new Error(`Fixture "${data.fixture}" requires ${state.parameters.phases} phases.`);
    return { type: 'result', version: fixtures.version, state };
  } catch (error) {
    return { type: 'error', error: { code: 'INVALID_INPUT', message: error.message, action: 'Choose fixture parameters shown by selected fixture.' } };
  }
}
if (typeof module !== 'undefined') module.exports = handleWorkbenchRequest;
if (typeof self !== 'undefined') self.handleWorkbenchRequest = handleWorkbenchRequest;
