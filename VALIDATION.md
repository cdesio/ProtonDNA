# Local validation — 23 October 2026

- TrackingGrid: Geant4 11.1.3 GUI-enabled build; 8 × 8 × 50 grid, 300 nm cubes, entrance z=0. Transportation-only proton produced 50 entries across all 50 layers, starting at z=0 and ending at the final layer entrance z=14.7 μm. Earlier full-physics check also recorded all 50 primary entries.
- DNA-simulation: built against stock Geant4 11.1.3. Replayed all 50 primary entries from that transport check with chemistry disabled. ROOT EventEdep contained 50 rows with original transport event and voxel IDs and PDG source code 2212.
- A one-entry replay with chemistry enabled completed successfully. This small 150 MeV sample produced no recorded direct/indirect damage hits; it does not establish nonzero chemistry yields or validate damage physics.
- Clustering: compiled with Python 3.12 and pybind11 2.13.6; imported the extension and generated strand-break and DSB-size CSV outputs for voxel 36 from the replay ROOT file. This exercised the empty-damage case.
- Initial inherited hybrid EM replay crashed in G4VEnergyLossProcess. The active DNA list now uses stock standard EM option4 plus the DNA Opt2 activator in the Target region. Physics changes require scientific validation before production use.
- Current GUI changes were compiled; previous Qt/OpenGL launches verified the grid view. No new GUI was launched during packaging.
