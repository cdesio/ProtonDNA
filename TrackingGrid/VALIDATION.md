# Local validation — 23 October 2026

Current simplified configuration, superseding earlier stopping-region experiments:

- GUI-enabled build succeeded against Geant4 11.1.3 in `/opt/anaconda3/envs/geant4.1`.
- Active geometry: 8 × 8 × 50 water voxels, 300 nm side, centred at the origin in vacuum. Absorber construction is commented out.
- Ran one 150 MeV proton with full physics using `macros/proton.mac`: 52 binary phase-space records, including 50 primary-proton voxel entries spanning all 50 longitudinal layers. Entry energies ranged from 150 MeV to approximately 149.98758 MeV.
- Transportation-only check produced exactly 50 primary-proton records, with constant 150 MeV energy.
- First entry z was −0.0075 mm and final-layer entry z approximately +0.0072 mm. Local positions were within the entered voxel bounds.
- Both runs completed without Geant4 exceptions. `.ps` files were parsed as 18 native doubles per record.
- Qt/OpenGL support remains enabled; the default GUI view still draws and fits `voxel_block`. GUI was not relaunched for this simplification.

Placement update: grid now spans x/y = −1.2 to +1.2 μm and z = 0 to 15 μm, with axes at the world origin. Rebuilt and verified 50 transportation-only entries, from z=0 to the last-layer entrance at z=14.7 μm. Earlier centred-z entry coordinates above describe the previous placement.
