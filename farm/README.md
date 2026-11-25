# BDX optical production on JLab ifarm / SWIF2

Each job simulates **one crystal type, one distance, and a few independent events**.
No Docker is used on the farm. Build once on an EL9 ifarm environment compatible
with Geant4 10.7.4, then share the frozen release among many jobs. Native JLab
modules/accounts must be checked on your own account; local container tests do
not establish remote access or SWIF2 credentials.

## Clone the BDX branch on ifarm

The BDX implementation is versioned in branch `optical_photons_simulation` of
`marcospreafico/clas12Tags`, starting from the original `frederic` baseline
`1d2313f8cc701e43bf461bee44ea12f43c7d402d`.

```bash
cd /SHARED/PATH
git clone --branch optical_photons_simulation --single-branch \
  https://github.com/marcospreafico/clas12Tags.git clas12Tags
cd clas12Tags
git rev-parse HEAD
```

To update an existing checkout of this branch, use `git pull --ff-only`.
The original `frederic` branch does not include this BDX implementation.
Compile a new native release after changing the source revision; queued jobs
continue to use their frozen runtime archive. Optical `--set` overrides are
independent job inputs. The build records Git commit/status and source hashes.

## Configure and build on ifarm

```bash
cd /SHARED/PATH/clas12Tags
cp farm/ifarm.env.example farm/ifarm.env
# Edit farm/ifarm.env for modules available on this account.
bash farm/build_bdx.sh farm/ifarm.env 2
```

Use only the compilation cores allocated to you. The example module versions
reproduce the local baseline but availability is site-dependent. Every worker
sources the saved config afresh: SWIF2 does not propagate the submitting shell's
environment.

A fresh `farm/build/TIMESTAMP/` contains build.log, source.json, release.json and
runtime.tar.gz. Compilation uses a separate source copy, preventing reuse of Mac
or Docker object files and preserving any existing local executable. The archive
contains the private GEMC, HIPO checker/reducer, worker, environment and source.
`release.json` records source, binary, archive and environment hashes.

## Plan many small jobs

Use a new directory on shared storage visible to SWIF2 and compute nodes, not the
ifarm machine's local `/tmp`. Replace ACCOUNT with your assigned JLab account.

```bash
python3 farm/plan_bdx.py \
  --release farm/build/TIMESTAMP/release.json \
  --account ACCOUNT --workflow bdx_bgo_scan \
  --output /SHARED/PATH/bdx_bgo_scan \
  --types 1 2 3 4 5 6 7 8 --distances-mm 20 80 140 200 \
  --events-per-job 5 --jobs-per-point 100 \
  --set light_yield_per_MeV=800
```

This prepares **3200 jobs / 16000 events**: 8 types × 4 distances × 100 replicas,
5 events each. Each job contains exactly one type/distance. Choose `--types 3`
for a workflow dedicated to type 3; change `--seed` between separate campaigns
whose results will be combined, so their random streams do not repeat.

For an exact event target use `--events-per-point 503 --events-per-job 5` instead
of `--jobs-per-point`: 101 jobs per point, with 3 events in the last job. Distinct
seeds are derived reproducibly from base seed, type, distance and replica; a seed
collision causes an error. Default bulk absorption is **10000 mm**, saved in each
job and HIPO Config bank.

Resource defaults: 1 core, 2 GB RAM, 2 GB local disk, 30 minutes, partition
`production`, constraint `el9`. Override `--ram`, `--disk`, `--time`, `--partition`
and `--constraint` after measuring representative jobs. High light yields and
large energy deposits increase optical-bank memory and output size.

Planning submits nothing. Inspect `plan.json` and the readable `commands.sh`.
The campaign snapshots its archive, bootstrap and job inputs. SWIF2 stages these
with `-input`, launches workers on local scratch, and copies validated results
with `-output`. `/farm_out/USER/WORKFLOW/` holds per-attempt stdout/stderr.

## Submit and monitor

```bash
python3 farm/plan_bdx.py --submit-plan /SHARED/PATH/bdx_bgo_scan
swif2 status bdx_bgo_scan
```

Submission verifies all saved inputs and creates an exclusive submission.json
marker. A second invocation is rejected. After a partial submission failure,
inspect submission.json and SWIF2; do not delete the marker and blindly resubmit,
as that could create repeated jobs. Automatic retry of individual SWIF attempts
uses the same saved job inputs/seed, which is intentional.

Each successful job creates `results/JOB.hipo`, `results/JOB.json` (receipt with
hashes/configuration identifiers and event count), and `results/JOB.logs.tar.gz`.
The worker checks all optical budgets, source links, phase-space coordinates,
event count, type, seed, distance, and saved optical parameters before publishing.
A failed validation returns nonzero; its log tail is included in SWIF stdout/stderr.
No job modifies campaign metadata shared with other jobs.

## Analyze many HIPO files together

Load the environment used to build the reader. For example, after a new login:

```bash
export BDX_RUNTIME_ROOT="$PWD/farm/build/TIMESTAMP/runtime"
source farm/environment.sh "$BDX_RUNTIME_ROOT/farm/ifarm.env"

python3 scripts/analyze_bdx_hipos.py \
  --native --reader "$BDX_RUNTIME_ROOT/build/bdx_optics/bdx_hipo" \
  --campaign /SHARED/PATH/bdx_bgo_scan \
  --output /SHARED/PATH/analysis_full_base

python3 scripts/analyze_bdx_hipos.py \
  --native --reader "$BDX_RUNTIME_ROOT/build/bdx_optics/bdx_hipo" \
  --campaign /SHARED/PATH/bdx_bgo_scan --square-mm 3 \
  --output /SHARED/PATH/analysis_square3mm
```

The analysis needs matplotlib (`requirements-bdx-optics.txt`) in its own Python
environment. Simulation/worker execution does not need matplotlib.
Multiple `--campaign` arguments combine compatible workflows. Every planned job
must have both HIPO and receipt; use `--allow-partial` explicitly to analyze only
completed results. Missing jobs are then listed and the analysis is marked partial.
Invalid files are never silently skipped. Run analysis after staging completes.

You can also pass quoted patterns or newline file lists, without farm metadata:

```bash
python3 scripts/analyze_bdx_hipos.py --native --reader /PATH/TO/bdx_hipo \
  '/SHARED/PATH/results/*.hipo' --output /SHARED/PATH/analysis_all

python3 scripts/analyze_bdx_hipos.py --native --reader /PATH/TO/bdx_hipo \
  --input-list inputs.list --output /SHARED/PATH/analysis_list
```

List paths are relative to the list file unless absolute. On the Mac, omit
`--native` to use the Docker HIPO reader; compile the updated reader first with
`bash scripts/build_bdx_optics.sh`. This mode starts a container per input file;
native analysis on ifarm is preferable for thousands of files.

For each (type, distance), the analysis sums selected scintillation arrivals,
generated scintillation photons and deposited MeV over **all events in all jobs**:

- efficiency = total selected / total generated;
- photons/MeV = total selected / total Edep;
- yield-normalized response = photons/MeV / configured light yield.

It never averages per-job efficiencies. Zero-arrival events remain in denominators;
SEM is calculated using particle events as independent clusters. No error is
invented for a single event. Wilson68 is counting precision only, not a model error.
The C++ reader validates full photon banks and applies the optional square mask,
then emits only event summaries to Python. There are no huge intermediate Flux
CSVs; memory retains event summaries, not individual photon rows.

The analysis rejects duplicated paths/content, repeated seed+configuration at a
scan point, incompatible optics/beam, different geometries for one type, and
mixed build releases when receipts are available. Raw HIPO files alone do not
identify the native archive, so retain their receipts for full provenance.
Outputs in a new directory: collection.csv, analysis.json (inputs/hashes, event
counts and selection), inputs.list, collection.png and collection.svg.

The whole base absorbs each photon at its **first arrival**. A later 3×3 selection
measures first-arrival acceptance; it does not reproduce reflection from paint
outside a physical small SiPM.

## Verification and official references

Local tests cover grid splitting/remainders, seed uniqueness, mocked SWIF2
submission, drift/duplicate submission rejection, partial campaigns, pooled ratios,
zero-arrival events, incompatible inputs and repeated seeds. Local container
checks exercise native release compilation, worker staging and actual HIPO reading;
remote module access, credentials, quotas and the first real job require ifarm.

- https://scicomp.jlab.org/cli/add-job.html
- https://scicomp.jlab.org/cli/swif.html
- https://scicomp.jlab.org/docs/swif2
- https://scicomp.jlab.org/scicomp/ (October 2026 /farm_out stdout/stderr policy)
