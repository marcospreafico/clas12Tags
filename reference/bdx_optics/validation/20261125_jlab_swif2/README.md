# JLab SWIF2 tooling validation — 25 November 2026

User request: native ifarm build, many small SWIF2 jobs dedicated to individual
BGO crystal types/distances, and one analysis across many HIPO files. Optical
properties and geometry are unchanged: BGO absorption length remains 10000 mm.

Local verification used the pinned GEMC EL9 container as a native build/test
environment. `farm/build_bdx.sh` compiled a fresh isolated source tree with
Geant4 10.7.4 and packaged a 6.9 MiB immutable runtime. The site config used
ccdb/1.99.6 and hipo/4.5.0 available in that container. The provided ifarm example
uses site-selectable modules; its exact availability on the user's JLab account
has not been verified.

13 unit tests passed, including planning/remainders, deterministic distinct seeds,
mocked SWIF2 invocation, input drift and duplicate-submission rejection, partial
campaigns, ratio-of-totals statistics, zero-arrival events, repeated files/seeds,
and incompatible optics/geometry/beam. The documented 3200-job example has 3200
unique derived seeds.

Four real workers were exercised with staged runtime/job files, including a
symlinked input to mimic SWIF2's input cache. They produced 6 events at one
distance for BGO types 1 and 8, with job sizes 2+1 per type. Every job's HIPO,
requested optical configuration and receipt were validated before publication.
The results were copied as SWIF2 output staging would do; no SWIF2 service was
contacted and no remote jobs were submitted.

The same four HIPO files were pooled on the Mac using the Docker reader and host
Matplotlib. Results: 27253 generated scintillation photons, 22739 whole-base first
arrivals, 105 arrivals in the central 3x3 mm² selection. Generated-photon and
energy denominators were identical under both selections and raw HIPO hashes
were unchanged. These are low-statistics software checks, not experimental
validation or final response predictions.

The native HIPO reducer also read a 16 MB event with 100000 optical first arrivals:
whole-base count 100000, square count 50000, and a second zero-arrival event was
retained. This verifies that photon-level CSV exports are not needed for pooled
analysis. The legacy HIPO check mode remains available.

The minimal simulation container has no Matplotlib: worker execution and HIPO
reduction succeeded without it; graph generation was verified in the host
analysis environment. Farm analysis requires Matplotlib as documented.

Official interfaces and policies checked while implementing:
- https://scicomp.jlab.org/cli/add-job.html
- https://scicomp.jlab.org/cli/swif.html
- https://scicomp.jlab.org/docs/swif2 (worker initializes its own environment)
- https://scicomp.jlab.org/scicomp/ (stdout/stderr must be under /farm_out)

Remote credentials, SWIF2 account, shared storage access, quotas and the first
real job remain to be verified on ifarm. Local mocked submission is not proof
of remote submission capability.
