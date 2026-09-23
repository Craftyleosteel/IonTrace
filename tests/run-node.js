/**
 * Console runner for the physics suite.
 *
 * The same modules the browser loads are imported here unchanged, so a green
 * run in the terminal and a green run in tests/index.html are testing exactly
 * the same code. Exits non-zero on failure so CI can gate on it.
 */

import { run } from './harness.js';
import './physics.test.js';
import './beamline.test.js';

const t0 = performance.now();

const { passed, failed, failures } = await run((e) => {
  if (e.type === 'suite') console.log(`\n${e.name}`);
  else if (e.type === 'pass') console.log(`  PASS  ${e.name}  (${e.ms.toFixed(0)} ms)`);
  else if (e.type === 'fail') console.log(`  FAIL  ${e.name}  (${e.ms.toFixed(0)} ms)`);
});

if (failures.length) {
  console.log('\n--- failures ---');
  for (const f of failures) {
    console.log(`\n${f.suite} > ${f.name}\n${f.error.message}`);
  }
}

const secs = ((performance.now() - t0) / 1000).toFixed(1);
console.log(`\n${passed} passed, ${failed} failed in ${secs}s`);
process.exit(failed ? 1 : 0);
