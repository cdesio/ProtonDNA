# TrackingGrid

A standalone Geant4 **11.1.3** simulation that fires a particle through a parameterised water voxel grid and records voxel-entry phase space. Uses stock Geant4, with no custom lithium physics.

## Build and run

Activate an environment containing Geant4 11.1.3 and its physics datasets:

```sh
cmake -S . -B build -DGeant4_DIR=/path/to/Geant4/lib/cmake/Geant4 -DWITH_GEANT4_UIVIS=ON
cmake --build build -j 4
cd build
./TrackingGrid -mac ../macros/proton.mac -gui
```

`-gui` selects Qt and opens an OpenGL view fitted to the voxel grid. The macro fires one proton at startup; clicking Run fires additional events. Omit `-gui` for batch runs. For a build without GUI support use `WITH_GEANT4_UIVIS=OFF`.

```sh
./TrackingGrid -mac ../macros/proton.mac -out proton -ps proton.ps -seed 1234
```

Defaults produce `tracking_grid.ps` and `tracking_grid.root` in the working directory. `-out` changes the output stem; `-ps` overrides the phase-space path (`-bin` remains an alias).

## Geometry and source

- World: vacuum (`G4_Galactic`), half-size 10 cm.
- Grid: water, centred on the z axis with its entrance face at z=0, **8 × 8 × 50** cubic voxels of side **300 nm**; overall extent **2.4 × 2.4 × 15 μm**.
- Example proton: 150 MeV, starts at `(0.15, 0.15, −10) μm`, travelling along +z, 10 μm before the grid front.
- No active absorber. A commented water-block construction is retained in `src/DetectorConstruction.cc`. The student can uncomment it and choose its dimensions/position, checking for overlaps.
- No stopping-depth calibration or event-specific grid placement.

Set dimensions and positions in the macro before `/run/initialize`:

```text
/det/setWorldHalfSize 10 cm
/det/setVoxelSize 300 nm
/det/setGrid 8 8 50
/det/setMaxStep 30 nm
/primary/position 0.15 0.15 -10 um
/primary/energy 150 MeV
/primary/direction 0 0 1
```

By default x and y span −1.2 to +1.2 μm and z spans 0 to 15 μm. The axes originate at (0,0,0), where the central x/y voxel boundaries meet the entrance face. `/det/setGridCentre x y z um` optionally overrides this default placement.

The source should remain outside the grid and inside the world. Move it upstream if the grid grows. Grid numbering is `ix + Nx * (iy + Ny * iz)`. `/det/checkOverlaps true` enables Geant4 overlap checks; parameterised checks may take time.

The muon macro supplies a directed 4 GeV muon example, not a realistic cosmic-ray distribution.

## Phase-space output

The `.ps` file is **binary**, with no header: one record of 18 native double-precision numbers for each first entry per event, track and voxel. No exit-state records are written. Re-entry to the same voxel is deduplicated; particles born inside a voxel have no entry record for their birth voxel.

Field order:

```text
localX, localY, localZ, dirX, dirY, dirZ, kineticEnergy,
eventID, PDGCode, voxelID, globalTime, PDGCode,
worldX, worldY, worldZ, seed, trackID, parentID
```

Binary units: **mm, MeV, ns**. Local coordinates are relative to the entered voxel. PDG codes include negative antiparticle codes; private-ID legacy replay readers need adapting. The repeated PDG field preserves the existing 18-double layout.

ROOT output also contains `Ntuples/primary`, `Ntuples/steps` (voxel energy deposits and step information), and `Ntuples/phase_space`. ROOT units are nm, eV and ns, as labelled. The application is serial. Use distinct output stems for independent runs; repeated interactive runs reuse the ROOT filename and append entry records to the open binary stream.

## Physics

Stock EM option4, decay, radioactive decay, hadron elastic physics, FTFP_BERT inelastic physics, stopping physics, ion physics and a step limiter. The production cut is 1 μm; maximum voxel step is independently set to 30 nm. This is condensed-history transport, not Geant4-DNA track structure. Validate physics and cuts for the intended study.

`-physics transport` disables interactions for geometry/entry-record checks.
