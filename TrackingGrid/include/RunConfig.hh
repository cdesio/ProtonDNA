#ifndef TRACKINGGRID_RUN_CONFIG_HH
#define TRACKINGGRID_RUN_CONFIG_HH

#include "globals.hh"

struct RunConfig {
  G4String outputStem = "tracking_grid";
  G4String binaryOutput;
  G4long seed = 1234;
};

#endif
