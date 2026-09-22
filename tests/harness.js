/**
 * A deliberately tiny test harness.
 *
 * There is no Node toolchain in this project by design, so tests run in the
 * browser from tests/index.html against the same ES modules the site loads.
 * That removes any chance of the tested code and the shipped code diverging.
 */

const suites = [];
let current = null;

export function describe(name, fn) {
  current = { name, tests: [] };
  suites.push(current);
  fn();
  current = null;
}

export function it(name, fn) {
  if (!current) throw new Error('it() must be called inside describe()');
  current.tests.push({ name, fn });
}

export function assert(condition, message) {
  if (!condition) throw new Error(message || 'Assertion failed');
}

/** Absolute-tolerance comparison. */
export function assertClose(actual, expected, tol, message = '') {
  const diff = Math.abs(actual - expected);
  assert(
    diff <= tol,
    `${message}\n  expected ${expected}\n  actual   ${actual}\n  |diff|   ${diff.toExponential(3)} > tol ${tol.toExponential(3)}`
  );
}

/** Relative-tolerance comparison, falling back to absolute near zero. */
export function assertRelClose(actual, expected, relTol, message = '') {
  const scale = Math.abs(expected);
  const tol = scale > 0 ? relTol * scale : relTol;
  assertClose(actual, expected, tol, message);
}

/**
 * Run every registered suite.
 *
 * @param {(event: object) => void} [report] Called as results arrive, so a UI
 *        can stream them rather than waiting for the whole run.
 */
export async function run(report = () => {}) {
  let passed = 0;
  let failed = 0;
  const failures = [];

  for (const suite of suites) {
    report({ type: 'suite', name: suite.name });
    for (const test of suite.tests) {
      const started = performance.now();
      try {
        await test.fn();
        const ms = performance.now() - started;
        passed++;
        report({ type: 'pass', suite: suite.name, name: test.name, ms });
      } catch (err) {
        const ms = performance.now() - started;
        failed++;
        failures.push({ suite: suite.name, name: test.name, error: err });
        report({
          type: 'fail',
          suite: suite.name,
          name: test.name,
          ms,
          message: err.message,
        });
      }
      // Yield so the page can paint between tests.
      await new Promise((r) => setTimeout(r, 0));
    }
  }

  report({ type: 'done', passed, failed });
  return { passed, failed, failures };
}

/** Attach an informational note to the running report. */
export function note(report, text) {
  report({ type: 'note', text });
}
