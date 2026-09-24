/**
 * The biased tube and the bare conductor.
 *
 * Both are simple enough that the interesting tests are about what they are
 * NOT: a tube that quietly changed the beam energy would be wrong about the
 * hardware, and a solid disc that let ions through would be wrong about the
 * one thing it exists to do.
 */

import { describe, it, assert, assertRelClose } from './harness.js';

import { Beamline } from '../src/beamline.js';
import { createElement } from '../src/elements/index.js';
import { createTube } from '../src/elements/tube.js';
import { createElectrode } from '../src/elements/electrode.js';
import { discBeam, makeIon } from '../src/ion.js';
import { flyBeam, flyIon, kineticEnergy } from '../src/integrator.js';
import { mmToM, mToMm, joulesToEV } from '../src/constants.js';

const ION = { mass: 100, charge: 1, energy: 50 };

/** One element between two drifts. */
function column(element) {
  return new Beamline([
    createElement('drift', { length: 15, bore: 6 }),
    element,
    createElement('drift', { length: 40, bore: 6 }),
  ]);
}

function transmitted(bl, radius = 1.0) {
  const { tracks } = flyBeam(bl, discBeam({ ...ION, count: 9, radius }), {
    cfl: 0.05,
    maxSteps: 400000,
  });
  return tracks.filter((t) => t.stop === 'exited').length;
}

/* ------------------------------------------------------------------ */

describe('Biased tube', () => {
  it('does no net work on an ion that gets through', () => {
    /*
      The claim the docstring makes, and the one a user is most likely to
      assume otherwise: a tube at -500 V does not make a 550 eV beam. Both ends
      of the domain are grounded, so whatever the entrance fringe takes the
      exit fringe gives back.
    */
    const tube = createTube({ voltage: -500 });
    const bl = column(tube);
    const { stop, points } = flyIon(bl, makeIon(ION), { cfl: 0.05 });
    assert(stop === 'exited', `an axial ion should get through, got ${stop}`);

    const last = points[points.length - 1];
    const out = joulesToEV(kineticEnergy(last, ION.mass));
    assertRelClose(out, ION.energy, 2e-3, 'exit energy should equal entry energy');
  });

  it('is at a different energy inside, which is the point of biasing it', () => {
    // If it were not, the element would be doing nothing at all.
    const tube = createTube({ voltage: -500 });
    const mid = tube.potentialAt(0, 0, tube.length / 2);
    assert(
      mid < -100,
      `the axis inside the tube should sit near its bias, got ${mid.toFixed(1)} V`
    );
  });

  it('focuses, and harder at higher voltage', () => {
    // A lens that transmits everything at every setting might simply be inert.
    const tube = createTube({ voltage: 0 });
    const bl = column(tube);
    const spread = (v) => {
      tube.setVoltage(v);
      const { tracks } = flyBeam(bl, discBeam({ ...ION, count: 9, radius: 2.0 }), {
        cfl: 0.05,
        maxSteps: 400000,
      });
      const out = tracks.filter((t) => t.stop === 'exited');
      const r = out.map((t) => {
        const p = t.points[t.points.length - 1];
        return Math.hypot(p.x, p.y ?? 0);
      });
      return Math.max(...r);
    };
    const off = spread(0);
    const on = spread(-400);
    assert(on < off, `biasing should narrow the beam: ${mToMm(on).toFixed(2)} mm vs ${mToMm(off).toFixed(2)} mm off`);
  });

  it('warns when its own end faces are clipping its fringe', () => {
    // No guard cylinders, so this element needs its margin more than any other.
    const tight = createTube({ margin: 4, housingRadius: 16 });
    assert(
      tight.warnings.some((w) => w.includes('decay lengths')),
      'a 4 mm margin on a 16 mm housing should warn'
    );
    const roomy = createTube({ margin: 24, housingRadius: 16 });
    assert(
      !roomy.warnings.some((w) => w.includes('decay lengths')),
      `a 24 mm margin should not warn, got: ${roomy.warnings.join(' | ')}`
    );
  });

  it('transmits from the toolbar defaults', () => {
    assert(transmitted(column(createElement('tube', {}))) === 9, 'defaults should transmit');
  });
});

/* ------------------------------------------------------------------ */

describe('Conductor', () => {
  it('stops the beam when it is a solid disc', () => {
    /*
      r0 = 0 means metal across the whole aperture. The bug this guards is
      specific: declaring `clearBore` would short-circuit the strike test on
      the axis, and an ion aimed at the middle of a beam stop would sail
      through the one element whose entire purpose is to stop it.
    */
    const stop = createElectrode({ r0: 0, r1: 14, voltage: 0 });
    assert(stop.clearBore === undefined, 'a solid disc must not promise a clear bore');
    assert(transmitted(column(stop)) === 0, 'a solid disc should stop every ion');
  });

  it('passes the beam when it is a plate with a hole', () => {
    const plate = createElectrode({ r0: 5, r1: 14, voltage: -75 });
    assert(plate.clearBore === mmToM(5), 'a plate with a hole has a clear bore');
    assert(transmitted(column(plate)) === 9, 'a 5 mm hole should pass a 1 mm beam');
  });

  it('builds a tube, a ring and a plate from the same element', () => {
    // The claim that makes this one element rather than three.
    const shapes = {
      plate: { z0: 18, z1: 20, r0: 5, r1: 14 },
      tube: { z0: 8, z1: 32, r0: 6, r1: 8 },
      ring: { z0: 18, z1: 22, r0: 10, r1: 13 },
    };
    for (const [name, geom] of Object.entries(shapes)) {
      const e = createElectrode({ ...geom, voltage: -60 });
      assert(e.rects.length === 1, `${name} should draw as one piece of metal`);
      assert(transmitted(column(e)) === 9, `${name} should pass a 1 mm beam on axis`);
    }
  });

  it('refuses a geometry that is not a conductor', () => {
    const refuses = (params, why) => {
      let threw = false;
      try {
        createElectrode(params);
      } catch {
        threw = true;
      }
      assert(threw, why);
    };
    refuses({ z0: 10, z1: 10 }, 'zero thickness is not a conductor');
    refuses({ r0: 8, r1: 4 }, 'an inner radius outside the outer one is not a shape');
    refuses({ z0: -2 }, 'metal outside the element is not placeable');
    refuses({ r1: 40, housingRadius: 16 }, 'metal outside the housing is not placeable');
  });

  it('drops the housing from a column solve when the metal reaches it', () => {
    /*
      A conductor taken out to the housing IS the wall there. Leaving a
      grounded housing part in the shared solve would put two potentials on one
      surface, and whichever was painted second would win along the whole run
      rather than just inside this element.
    */
    const across = createElectrode({ r0: 5, r1: 16, housingRadius: 16, voltage: -100 });
    assert(
      !across.parts.some((p) => p.name === 'housing'),
      'a conductor spanning to the wall should not also claim the wall'
    );
    assert(
      across.warnings.some((w) => w.includes('interrupts the grounded wall')),
      'and it should say so'
    );

    const inside = createElectrode({ r0: 5, r1: 13, housingRadius: 16, voltage: -100 });
    assert(
      inside.parts.some((p) => p.name === 'housing'),
      'a conductor clear of the wall keeps the housing'
    );
  });

  it('does not warn about its own margin at the shipped defaults', () => {
    // An element that complains the moment it is placed teaches the user to
    // ignore its warnings, and the next one will matter.
    const e = createElement('electrode', {});
    assert(
      !e.warnings.some((w) => w.includes('decay lengths')),
      `defaults should be roomy enough, got: ${e.warnings.join(' | ')}`
    );
  });
});
