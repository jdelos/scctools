importScripts('fixtures.js', 'worker-core.js');

self.onmessage = ({ data }) => self.postMessage(handleWorkbenchRequest(data, WORKBENCH_FIXTURES));
