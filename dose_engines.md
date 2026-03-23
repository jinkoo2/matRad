# matRad Dose Engines

Three photon/particle dose engines available in `matRad/doseCalc/+DoseEngines/`.

## Quick Comparison

| Feature | SVPB | ompMC | TOPAS |
|---------|------|-------|-------|
| Type | Analytical pencil beam | Monte Carlo (OpenMP) | Monte Carlo (Geant4) |
| Radiation | Photons only | Photons only | Photons, protons, helium, carbon, VHEE |
| Speed | Seconds | Minutes | Minutes–hours |
| Accuracy | Good | Better | Best (full Geant4 physics) |
| External exec | None | MEX compiled | TOPAS binary |
| Biological dose (RBE/LET) | No | No | Yes |
| Remote execution | No | No | Yes (batch export) |

---

## 1. SVPB — SVD Photon Pencil Beam

**File:** `matRad_PhotonPencilBeamSVDEngine.m`
**Inherits:** `matRad_PencilBeamEngineAbstract` → `matRad_DoseEngineBase`
**Algorithm:** Bortfeld et al. SVD pencil beam convolution (PMID: 8497215)

### How It Works
1. Loads three kernel components from machine data: `kernel1/2/3`
2. Each component parameterized as: `D_i = β_i/(β_i - m) × (exp(-m·d) - exp(-β_i·d))`
3. Convolves kernels with Gaussian penumbra (FWHM from machine data) via FFT
4. Optionally convolves with custom primary photon fluence
5. Stores as `griddedInterpolant` for fast per-bixel lookup

### Key Properties

| Property | Default | Description |
|----------|---------|-------------|
| `kernelCutOff` | config | Lateral kernel cutoff [mm] |
| `intConvResolution` | 0.5 | Kernel convolution grid resolution [mm] |
| `enableDijSampling` | true | Spatial sampling of dij |
| `dijSampling.relDoseThreshold` | 0.01 | Relative dose threshold for core region |
| `dijSampling.latCutOff` | 20 | Lateral sampling cutoff [mm] |
| `dijSampling.type` | `'radius'` | `'radius'` or `'dose'` |
| `dijSampling.deltaRadDepth` | 5 | Radiological depth clustering [mm] |
| `useCustomPrimaryPhotonFluence` | config | Use custom fluence vs. Gaussian |
| `randomSeed` | 0 | Mersenne Twister seed for bixel sampling |

### Required Machine Data
- `machine.data.betas` — beta parameters for 3 kernel components
- `machine.data.m` — absorption parameter
- `machine.data.energy` — beam energy
- `machine.data.kernel` — 3D depth-dose kernel matrices at multiple SSDs
- `machine.data.kernelPos` — kernel radial positions [mm]
- `machine.data.penumbraFWHMatIso` — penumbra FWHM at isocenter [mm]
- `machine.meta.SAD` — source-to-axis distance [mm]
- `machine.meta.SCD` — source-to-collimator distance [mm]

### Usage
```matlab
pln.propDoseCalc.engine = 'SVPB';  % or set via factory
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
```

---

## 2. ompMC — Photon OpenMP Monte Carlo

**File:** `matRad_PhotonOmpMCEngine.m`
**Inherits:** `matRad_MonteCarloEngineAbstract` → `matRad_DoseEngineBase`
**Algorithm:** C-based OpenMP MC photon transport via `omc_matrad` MEX function

### How It Works
1. Converts CT HU to one of 4 material indices (AIR, LUNG, TISSUE, BONE) with density
2. Defines beamlet source geometry (corners at SCD plane, Gaussian source width)
3. Calls compiled MEX: `omc_matrad(geom, source, options)` → returns 3D dose array
4. Accumulates into sparse `dij` matrix per beamlet
5. Applies absolute calibration factor

### Key Properties

| Property | Default | Description |
|----------|---------|-------------|
| `numHistoriesPerBeamlet` | config | MC histories per beamlet |
| `numHistoriesDirect` | config | Histories for direct dose calc |
| `absCalibrationFactor` | 3.49056e12 | Calibration to 1 Gy at 5 cm depth, 5×5 cm², SSD=900 mm |
| `useCornersSCD` | true | Use bixel corners at SCD plane |
| `outputMCvariance` | config | Store statistical variance |
| `relativeDosimetricCutOff` | config | Relative dose threshold for storage |
| `scale` | 10 | mm → cm unit conversion for MC geometry |

### Material HU Ranges

| Material | HU Range | Density Range [g/cm³] |
|----------|----------|-----------------------|
| AIR700ICRU | −1024 to −974 | 0.001 – 0.044 |
| LUNG700ICRU | −974 to −724 | 0.044 – 0.302 |
| ICRUTISSUE700ICRU | −724 to 101 | 0.302 – 1.101 |
| ICRPBONE700ICRU | 101 to 1976 | 1.101 – 2.088 |

### External Dependency
- **MEX function:** `thirdParty/ompMC/omc_matrad.[mexext]`
- Auto-compiled from C source with OpenMP flags if missing
- Physics data: `thirdParty/ompMC/data/`, PEGS4 materials: `thirdParty/ompMC/pegs4/`
- Photon spectrum: `thirdParty/ompMC/spectra/mohan6.spectrum` (Varian 6 MV)

### MC Simulation Defaults
- `nSplit` = 20, `nBatches` = 10
- `randomSeeds` = [97, 33]
- `global_ecut` = 0.7 MeV (electron cutoff), `global_pcut` = 0.010 MeV

### Usage
```matlab
pln.propDoseCalc.engine = 'ompMC';
pln.propDoseCalc.numHistoriesPerBeamlet = 5e4;
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
```

---

## 3. TOPAS — Geant4 Monte Carlo Interface

**File:** `matRad_TopasMCEngine.m` (~2700 lines)
**Inherits:** `matRad_MonteCarloEngineAbstract` → `matRad_DoseEngineBase`
**Algorithm:** Full Geant4-based MC transport; writes TOPAS input files, executes binary, reads results

### How It Works
1. Exports CT as binary RSP (Relative Stopping Power) cube + geometry file
2. Generates TOPAS text parameter files from templates (`matRad/doseCalc/topas/`)
3. Optionally executes TOPAS locally or exports files for remote cluster execution
4. Reads back binary/CSV dose (and optional LET/RBE) output
5. Assembles into `dij` sparse matrix or direct dose cubes

### Key Properties

| Property | Default | Description |
|----------|---------|-------------|
| `topasExecCommand` | system path | Command to run TOPAS binary |
| `workingDir` | `{userFolder}/TOPAS/` | Working directory for simulation files |
| `externalCalculation` | `'off'` | `'off'`=local, `'write'`=export files only, folder path=read results |
| `numHistoriesPerBeamlet` | config | Histories per beamlet |
| `numParticlesPerWeight` | 1e6 | Particles per MC weight |
| `minRelWeight` | 1e-5 | Discard beamlets below threshold |
| `numOfRuns` | 1 | Number of simulation batches |
| `numThreads` | 0 | 0=auto (all CPU cores) |
| `parallelRuns` | false | Run beams in parallel |

### Beam Profile Options

| Mode | Profile | Description |
|------|---------|-------------|
| Photons | `'uniform'` | Rectangular flat field |
| Particles | `'biGaussian'` | Bi-Gaussian emittance model |
| Any | `'virtualGaussian'` | Virtual Gaussian |
| Any | `'phasespace'` | Phase-space file source |
| Any | `'mlc'` | MLC-shaped field |

### Scorer Options

| Option | Default | Description |
|--------|---------|-------------|
| `scorer.doseToMedium` | true | Primary dose output |
| `scorer.doseToWater` | false | Dose-to-water scoring |
| `scorer.calcDij` | false | Full dose influence matrix |
| `scorer.LET` | false | Linear energy transfer |
| `scorer.RBE` | false | RBE-weighted dose |
| `scorer.outputType` | `'binary'` | `'binary'` or `'csv'` |

### RBE Models (ions)
- **Protons:** MCN (McNamara), WED (Wedenberg)
- **Carbon/Helium:** LEM1, libamtrack

### Material Conversion

| Option | Default | Values |
|--------|---------|--------|
| `materialConverter.mode` | `'HUToWaterSchneider'` | `'RSP'`, `'HUToWaterSchneider'` |
| `materialConverter.densityCorrection` | `'Schneider_TOPAS'` | `'Schneider_TOPAS'`, `'Schneider_matRad'` |
| `materialConverter.HUSection` | `'advanced'` | `'default'`, `'advanced'` |

### Physics Lists by Modality

| Modality | Physics Modules |
|----------|----------------|
| Photons | `g4em-standard_opt4`, `g4h-phy_QGSP_BIC_HP`, `g4decay` |
| Protons | + `g4h-elastic_HP`, `g4stopping`, `g4ion-QMD`, `g4radioactivedecay` |
| VHEE | + `g4ion-binarycascade`, `g4h-elastic_HP`, `g4stopping` |

### Generated File Structure
```
workingDir/
├── matRad_cube.txt / .dat       # CT geometry + binary voxel data
├── matRad_plan.txt              # Main simulation parameter file
├── world/                       # Geometry definition
├── beamSetup/                   # Beam profile templates
├── materialConverter/           # Schneider HU-to-material tables
├── scorer/                      # Dose/LET/RBE scorer definitions
└── output/
    ├── dose_[beam]_[run].bin    # Binary float32 dose
    └── variance_...             # Optional statistical variance
```

### Usage
```matlab
pln.propDoseCalc.engine = 'TOPAS';
pln.propDoseCalc.topasExecCommand = '/path/to/topas';
pln.propDoseCalc.scorer.LET = true;
pln.propDoseCalc.scorer.RBE = true;

% Local execution
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);

% Export files for cluster, then read results
pln.propDoseCalc.externalCalculation = 'write';
matRad_calcDoseInfluence(ct, cst, stf, pln);
% ... run TOPAS on cluster ...
pln.propDoseCalc.externalCalculation = '/path/to/results/';
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
```

---

## Base Class Hierarchy

```
matRad_DoseEngineBase
├── matRad_PencilBeamEngineAbstract
│   └── matRad_PhotonPencilBeamSVDEngine    (SVPB)
└── matRad_MonteCarloEngineAbstract
    ├── matRad_PhotonOmpMCEngine             (ompMC)
    └── matRad_TopasMCEngine                 (TOPAS)
```

## Engine Selection via Factory

```matlab
% Automatic selection from pln
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);

% Manual via pln
pln.propDoseCalc.engine = 'SVPB';   % or 'ompMC', 'TOPAS'
```

---

## pyMatRad Python Port — Implementation Notes

Python port lives in `pyMatRad/matRad/doseCalc/DoseEngines/`.  The three MATLAB
engines are re-implemented as:

| MATLAB class | Python class | File |
|---|---|---|
| `matRad_PhotonPencilBeamSVDEngine` | `PhotonPencilBeamSVDEngine` | `photon_svd_engine.py` |
| `matRad_PhotonOmpMCEngine` | `PhotonOmpMCEngine` | `photon_ompc_engine.py` |
| `matRad_TopasMCEngine` | `TopasMCEngine` | `topas_mc_engine.py` |

### Python ompMC — Key Differences from MATLAB

The Python ompMC is **not** a Monte Carlo engine.  It is an analytical approximation
using the same TERMA (Total Energy Released per unit MAss) framework:

- **Primary dose:** inverse-square law × exponential attenuation × NIST `μ_en/ρ`
- **Lateral profile:** erf-product (rect bixel convolved with Gaussian penumbra)
- **Scatter correction:** depth-dependent Compton term `f_s × (1 − exp(−d/d₀))`

Physical constants (6 MV, effective ~2 MeV beam, NIST XCOM):

| Constant | Value | Units |
|---|---|---|
| `μ_total/ρ` | 0.0497 | cm²/g |
| `μ_en/ρ` | 0.0270 | cm²/g |
| Scatter fraction | 0.28 | — |
| Scatter buildup depth | 80 | mm |

### Calibration Factor — Bug and Fix (2026-03-23)

**Bug:** `ABS_CALIBRATION_FACTOR = 3.49056e12` was copied from the MATLAB ompMC engine,
where it converts MC *histories* to Gy.  The Python analytical model computes
dimensionless physics ratios, not photon histories — so this constant was off by ~1.5 × 10⁸.

Result: ompMC max dose 706 million Gy/fx vs SVPB 4.7 Gy/fx.

**Fix:** Empirically recalibrate to match SVPB at reference conditions
(water phantom, uniform bixel weights):

```python
# photon_ompc_engine.py
ABS_CALIBRATION_FACTOR = 23220.0   # analytical model calibration

# effective calibration per bixel (bixelWidth = 5 mm):
calib = ABS_CALIBRATION_FACTOR * (bixelWidth / 50.0) ** 2
      = 23220 * 0.01 = 232
```

Derivation: `effective_old / ratio = 3.49e10 / 1.503e8 = 232`;
without the `(bixelWidth/50)²` factor: `ABS_CALIBRATION_FACTOR = 232 / 0.01 = 23220`.

**Validated results after fix:**

| Phantom | SVPB max (Gy/fx) | ompMC max (Gy/fx) | Ratio | ompMC speed |
|---------|-----------------|------------------|-------|------------|
| Water phantom | 4.709 | 4.699 | 0.998 | 3.3× faster |
| TG119 | 5.496 | 5.444 | 0.991 | 2.7× faster |

Per-voxel median ratio: 0.987 (water), 0.922 (TG119).  TG119 difference is expected —
the exponential model does not capture scatter from heterogeneities as accurately as
SVPB's SVD kernels.

### Parallelism

Python SVPB uses `ProcessPoolExecutor` (one worker per beam); ompMC does the same.
The TOPAS engine writes input files then calls the TOPAS binary via `subprocess`.

### CST Loading Pitfall (scipy.io)

When loading MATLAB `cst` (an `N×6` cell array) with `scipy.io.loadmat(squeeze_me=True)`,
the result is a NumPy object array of shape `(N, 6)`.

**Wrong:** `for row in cst_m.flat` — iterates all `N×6` individual cells in C order.
**Correct:** `for i in range(cst_m.shape[0]): row = cst_m[i]` — iterates N rows.
