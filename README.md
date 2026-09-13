# physics

Standalone physics explainers — one directory per project. Each holds a
finished document or app plus whatever source or code produces it, and each has
its own README. Built PDFs are committed so they can be read without installing
anything.

**Live at [drandyhaas.github.io/physics](https://drandyhaas.github.io/physics/)** —
the posters and broadsides read in the browser, and Field bench runs there.

| project | what it is |
| --- | --- |
| [**e8-heterotic-posters**](e8-heterotic-posters/) | Two single-page posters on E₈×E₈ heterotic string theory: the symmetry-breaking cascade from the Planck scale down to the Standard Model, and the complete ten-dimensional Lagrangian field by field. Authored as standalone HTML, exported to PDF. |
| [**EMviewer**](EMviewer/) | **Field bench** — an interactive map of **E**, **B** and **S** = E×B/μ₀ everywhere in the plane around a battery–wire–LED circuit, where the wire can be grabbed and reshaped and the fields re-solve as you drag. Shows E meeting an ideal wire at a right angle, and energy flowing to the load through the space beside the wire rather than inside it. Plain HTML, CSS and JavaScript: no build step, no dependencies. |
| [**string-copper-correspondence**](string-copper-correspondence/) | A seven-page broadside on why an open-string scattering amplitude is, integrand for integrand, the thermal-noise statistics of a two-dimensional resistive copper plane. Every number it quotes is produced by a verification script in the same directory, not cited. |

## Conventions

- **One directory per project**, self-contained, with its own README.
- **The built artefact is committed** — a PDF, or a page that runs as-is —
  alongside its source.
- **Numbers are measured, not quoted.** Where a document asserts a numerical
  claim, the directory contains the script that produces it, and the document
  reads its values from that script's output so the two cannot drift apart.
