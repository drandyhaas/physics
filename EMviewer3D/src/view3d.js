/*
 * view3d.js — orbit camera and perspective projection.
 *
 * No matrices: the camera is an orbit (azimuth, elevation, distance) about a
 * target, and a frame is an orthonormal basis plus a focal length in pixels.
 * Everything the renderer needs is project() one way and unproject() the other,
 * the second being what lets a handle be dragged in the plane of the screen at
 * its own depth.
 *
 * World convention matches the planar bench: the default circuit lies in z = 0,
 * so z is "out of the old picture" and elevation is measured from that plane.
 *
 * Usable both as a browser global (window.V3) and as a CommonJS module.
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.V3 = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const RAD = Math.PI / 180;
  const clamp = (v, a, b) => (v < a ? a : v > b ? b : v);

  // Elevation is clamped short of the poles: looking straight down the world up
  // axis leaves the right vector undefined.
  const EL_MAX = 89.2;

  /**
   * cam: { az, el, dist, fov, target:{x,y,z} }   az/el/fov in degrees
   * Returns a frame with the camera basis and the pixel focal length.
   */
  function frame(cam, W, H) {
    const az = cam.az * RAD, el = clamp(cam.el, -EL_MAX, EL_MAX) * RAD;
    const ce = Math.cos(el), se = Math.sin(el);
    const t = cam.target;
    const dir = { x: ce * Math.cos(az), y: ce * Math.sin(az), z: se };   // target -> eye
    const eye = { x: t.x + dir.x * cam.dist, y: t.y + dir.y * cam.dist, z: t.z + dir.z * cam.dist };
    const fwd = { x: -dir.x, y: -dir.y, z: -dir.z };                     // eye -> target

    // right = fwd x worldUp, normalised; worldUp is +z
    let rx = fwd.y * 1 - fwd.z * 0, ry = fwd.z * 0 - fwd.x * 1, rz = 0;
    const rn = Math.hypot(rx, ry, rz) || 1e-12;
    rx /= rn; ry /= rn; rz /= rn;
    const right = { x: rx, y: ry, z: rz };
    const up = {                                                          // right x fwd
      x: ry * fwd.z - rz * fwd.y,
      y: rz * fwd.x - rx * fwd.z,
      z: rx * fwd.y - ry * fwd.x
    };
    // Focal length off the SHORTER side, so the field of view is the one that
    // fits: taking it from the height alone makes a tall narrow canvas zoom in
    // horizontally and push the circuit out through the sides.
    const focal = (Math.min(W, H) / 2) / Math.tan(cam.fov * RAD / 2);
    return { eye, fwd, right, up, focal, cx: W / 2, cy: H / 2, near: cam.dist * 0.02 };
  }

  /**
   * World point -> [screenX, screenY, depth]. Depth is distance along the view
   * axis, so it doubles as the painter's-algorithm key. Returns null behind the
   * camera, which every caller must handle rather than drawing a mirrored ghost.
   */
  function project(f, x, y, z) {
    const dx = x - f.eye.x, dy = y - f.eye.y, dz = z - f.eye.z;
    const zv = dx * f.fwd.x + dy * f.fwd.y + dz * f.fwd.z;
    if (zv < f.near) return null;
    const xv = dx * f.right.x + dy * f.right.y + dz * f.right.z;
    const yv = dx * f.up.x + dy * f.up.y + dz * f.up.z;
    const k = f.focal / zv;
    return [f.cx + xv * k, f.cy - yv * k, zv];
  }

  // Screen point at a given depth -> world point. Inverse of project().
  function unproject(f, sx, sy, depth) {
    const k = depth / f.focal;
    const xv = (sx - f.cx) * k, yv = -(sy - f.cy) * k;
    return {
      x: f.eye.x + f.fwd.x * depth + f.right.x * xv + f.up.x * yv,
      y: f.eye.y + f.fwd.y * depth + f.right.y * xv + f.up.y * yv,
      z: f.eye.z + f.fwd.z * depth + f.right.z * xv + f.up.z * yv
    };
  }

  // Depth of a world point along the view axis, without projecting it.
  function depthOf(f, x, y, z) {
    return (x - f.eye.x) * f.fwd.x + (y - f.eye.y) * f.fwd.y + (z - f.eye.z) * f.fwd.z;
  }

  // An orthonormal pair spanning the plane with the given normal.
  function planeAxes(n) {
    const ax = Math.abs(n.x), ay = Math.abs(n.y), az = Math.abs(n.z);
    const seed = ax < ay && ax < az ? { x: 1, y: 0, z: 0 }
               : ay < az           ? { x: 0, y: 1, z: 0 }
                                   : { x: 0, y: 0, z: 1 };
    let ux = seed.y * n.z - seed.z * n.y,
        uy = seed.z * n.x - seed.x * n.z,
        uz = seed.x * n.y - seed.y * n.x;
    const un = Math.hypot(ux, uy, uz) || 1e-12;
    ux /= un; uy /= un; uz /= un;
    return [
      { x: ux, y: uy, z: uz },
      { x: n.y * uz - n.z * uy, y: n.z * ux - n.x * uz, z: n.x * uy - n.y * ux }
    ];
  }

  return { RAD, EL_MAX, clamp, frame, project, unproject, depthOf, planeAxes };
});
