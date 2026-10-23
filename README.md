# ProtonDNA

A three-stage proton transport and DNA-damage workflow using **Geant4 11.1.3**:

1. **TrackingGrid** transports particles through a water-voxel grid and writes voxel-entry phase space.
2. **DNA-simulation** replays each phase-space entry into a 300 nm DNA target and writes physical/chemical interaction data to ROOT.
3. **Clustering** analyses the DNA output for strand-break clusters, selecting one transport voxel at a time.

The geometry is deliberately simple: the user chooses source positions, dimensions and any absorber. TrackingGrid's optional absorber is retained as commented code. See [source provenance and inherited notices](SOURCE_NOTES.md).

## Contents

| Folder | Purpose |
|---|---|
| `TrackingGrid/` | Transport executable and proton/muon example macros |
| `DNA-simulation/` | DNA transport/chemistry executable `rbe` |
| `Clustering/` | Python analysis and pybind11 clustering extension |
| `geometryFiles/` | Matching sugar/histone geometry, templates and generators |

## Requirements and builds

Use Geant4 **11.1.3**, including its datasets and Qt/OpenGL for the GUI. Set `Geant4_DIR` to the directory containing `Geant4Config.cmake`. Activate that environment before running the executables so the dataset variables are set.

From the repository root:

```sh
cmake -S TrackingGrid -B TrackingGrid/build -DGeant4_DIR="$Geant4_DIR" -DWITH_GEANT4_UIVIS=ON
cmake --build TrackingGrid/build -j 4
cmake -S DNA-simulation -B DNA-simulation/build -DGeant4_DIR="$Geant4_DIR"
cmake --build DNA-simulation/build -j 4
```

For clustering, use Python 3.12 with the packages in `Clustering/requirements.txt`. A separate Python environment is fine:

```sh
python -m pip install -r Clustering/requirements.txt
cmake -S Clustering -B Clustering/build -DPython3_EXECUTABLE="$(command -v python)"
cmake --build Clustering/build -j 4
```

## 1. Generate phase space

```sh
mkdir -p outputs
./TrackingGrid/build/TrackingGrid -mac TrackingGrid/macros/proton.mac -out outputs/proton -seed 1234
```

The default grid is **8 × 8 × 50** water voxels, each 300 nm: x/y = −1.2 to +1.2 μm, z = 0 to 15 μm. The 150 MeV proton starts at `(0.15, 0.15, −10) μm` and travels along +z through vacuum into the grid. No absorber is active. The macro runs one proton; edit `/run/beamOn` for more.

GUI, from the build directory:

```sh
cd TrackingGrid/build
./TrackingGrid -mac ../macros/proton.mac -gui
```

The camera is fitted to the grid and the axes originate at `(0,0,0)` on its entrance face. [Transport options and binary field layout](TrackingGrid/README.md).

`outputs/proton.ps` contains headerless **18-double records with PDG particle codes**, in mm/MeV/ns. There are entry records only, deduplicated per event/track/voxel. ROOT step and entry data are also written.

## 2. Replay into DNA

From the repository root, in the Geant4 environment:

```sh
./DNA-simulation/build/rbe -mac DNA-simulation/rbe_PS.in \
  -ps outputs/proton.ps \
  -sugar geometryFiles/sugarPos_n2_300nm_H75nm_R11nm_33hist_50deg.bin \
  -histone geometryFiles/histonePos_n2_300nm_H75nm_R11nm_33hist_50deg.bin \
  -out outputs/dna -seed 1234
```

The executable automatically replays one event per PS record; do not add `/run/beamOn` to the replay macro. `-chemOFF` disables chemistry for a physical-stage check. Target size must match the transport voxel size and supplied DNA geometry. The supplied geometry is 300 nm; changing that requires appropriate DNA inputs too.

Default DNA physics uses standard EM option4 with the stock Geant4-DNA **Opt2 region activator** in `Target`, plus the inherited application chemistry. `-DNAopt7` selects standalone DNA option7 instead and needs separate energy/species validation. Do not assume every cosmic-ray species supported by TrackingGrid has a valid DNA replay model.

Each transport entry is replayed independently in local voxel coordinates. This is a staged workflow, not a single continuous microscopic trajectory through all voxels. Entry event ID and voxel ID are retained in `part1_EventNum` and `part1_CopyNum`. Clustering retains the inherited replay-event grouping; it does not automatically combine all tracks belonging to one original transport primary.

## 3. Get the Damage

Activate the clustering Python environment, then from the repository root:

```sh
python Clustering/run.py --filename outputs/dna.root \
  --output outputs/damage_voxel36.csv --copy 36 \
  --sugar geometryFiles/sugarPos_n2_300nm_H75nm_R11nm_33hist_50deg.bin
```

`--copy` is the transport voxel ID (`ix + Nx * (iy + Ny * iz)`). The central example trajectory starts in voxel 36; scattering may move it into other columns. Use IDs present in the transport output. The supplied DNA geometry uses discontinuous indexing; `--continuous` is for a separately generated continuous geometry. Thresholds default to 5–37.5 eV for direct damage and 0.405 probability for indirect damage; they are inherited modelling parameters, not fitted values for this project.

 Build directories, caches and generated output are excluded from Git.
