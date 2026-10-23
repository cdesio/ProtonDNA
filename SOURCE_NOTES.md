# Source provenance

TrackingGrid was adapted from the local parameterised voxel phase-space transport application. The active project contains no neutron-capture source, custom lithium processes or endpoint-calibration geometry.

`DNA-simulation/`, `Clustering/` and `geometryFiles/` were copied from AlphaGlue_Collab (https://github.com/Lujinabu/AlphaGlue_Collab), source checkout commit `73baaadd5b9dc1c796da971e4101e31f510d821f` (local working files, which may contain changes beyond that commit).

Changes for ProtonDNA include explicit 18-double PDG phase-space replay, generic `-ps` naming, strict file-length/read checks, initialisation of both source pointers, Geant4 11.1.3 version requirements, stock standard EM option4 plus Geant4-DNA Opt2 in the Target region (replacing the inherited custom hybrid EM constructor), build portability and clustering output-path fixes.

Original Geant4 copyright, licence notices and scientific attribution are retained in the imported source files. No new blanket licence is applied to third-party source. Consult the retained notices and upstream authors for redistribution terms. Geant4 licence: https://geant4.web.cern.ch/download/license

The inherited chemistry implementation contains application-specific reaction and IRT classes. The transport and DNA stages use different physics resolutions and require scientific validation for each study. Build and smoke-test results are not a validation of damage yields.
