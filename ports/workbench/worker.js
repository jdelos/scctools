importScripts('fixtures.js');

self.onmessage = ({ data }) => {
  try {
    const state = WORKBENCH_FIXTURES.states[data.fixture];
    if (!state) throw new Error(`Unknown fixture "${data.fixture}". Choose symbolic or numeric.`);
    const stages = Number(data.stages);
    const phases = Number(data.phases);
    if (!Number.isInteger(stages) || stages < 1) throw new Error('Stages must be positive integer.');
    if (phases !== 2) throw new Error('Initial execution scope supports exactly 2 phases.');
    self.postMessage({ type: 'result', version: WORKBENCH_FIXTURES.version, state });
  } catch (error) {
    self.postMessage({ type: 'error', error: { code: 'INVALID_INPUT', message: error.message, action: 'Choose valid fixture, positive stages, and 2 phases.' } });
  }
};
