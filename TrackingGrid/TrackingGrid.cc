#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "RunConfig.hh"
#include "TrackingOptions.hh"

#include "G4RunManagerFactory.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#ifdef TRACKINGGRID_WITH_GEANT4_UIVIS
#include "G4VisExecutive.hh"
#endif

#include <cstdlib>
#include <iostream>
#include <sstream>
#include "G4SystemOfUnits.hh"

namespace {

void PrintUsage(const char* program)
{
  std::cout << "Usage: " << program
            << " [-mac macro.mac] [-out output_stem] [-seed integer] [-gui]"
            << " [-physics full|transport] [-ps phase_space.ps]\n";
}

}  // namespace

int main(int argc, char** argv)
{
  G4String macro;
  G4bool gui = false;
  G4String physicsMode = "full";
  RunConfig config;

  for (G4int i = 1; i < argc; ++i) {
    const G4String arg = argv[i];
    if ((arg == "-h") || (arg == "--help")) {
      PrintUsage(argv[0]);
      return 0;
    } else if (arg == "-mac" && i + 1 < argc) {
      macro = argv[++i];
    } else if (arg == "-out" && i + 1 < argc) {
      config.outputStem = argv[++i];
    } else if ((arg == "-ps" || arg == "-bin") && i + 1 < argc) {
      config.binaryOutput = argv[++i];
    } else if (arg == "-seed" && i + 1 < argc) {
      config.seed = std::atol(argv[++i]);
    } else if (arg == "-gui") {
      gui = true;
    } else if (arg == "-physics" && i + 1 < argc) {
      physicsMode = argv[++i];
    } else {
      std::cerr << "Unknown or incomplete argument: " << arg << "\n";
      PrintUsage(argv[0]);
      return 1;
    }
  }

  if (physicsMode != "full" && physicsMode != "transport") {
    std::cerr << "Unknown physics mode: " << physicsMode << "\n";
    PrintUsage(argv[0]);
    return 1;
  }

  if (config.binaryOutput.empty()) {
    config.binaryOutput = config.outputStem + ".ps";
  }

#ifndef TRACKINGGRID_WITH_GEANT4_UIVIS
  if (gui) {
    std::cerr << "GUI support is disabled. Reconfigure with -DWITH_GEANT4_UIVIS=ON.\n";
    return 1;
  }
#endif

  CLHEP::HepRandom::setTheSeed(config.seed);

  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  auto* detector = new DetectorConstruction();
  auto* trackingOptions = new TrackingOptions();
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new PhysicsList(physicsMode == "transport"));
  runManager->SetUserInitialization(new ActionInitialization(config, detector, trackingOptions));

#ifdef TRACKINGGRID_WITH_GEANT4_UIVIS
  G4VisExecutive* visManager = nullptr;
  if (gui) {
    visManager = new G4VisExecutive();
    visManager->Initialize();
  }
#endif

  G4UIExecutive* uiExec = nullptr;
  if (gui) {
    uiExec = new G4UIExecutive(argc, argv, "Qt");
  }

  auto* ui = G4UImanager::GetUIpointer();
  G4int commandStatus = 0;
  if (gui) {
    ui->ApplyCommand("/tracking/storeTrajectory 1");
  }
  if (!macro.empty()) {
    commandStatus = ui->ApplyCommand("/control/execute " + macro);
  }

#ifdef TRACKINGGRID_WITH_GEANT4_UIVIS
  // Batch macros need no visualisation commands: -gui supplies a default viewer.
  // Preserve any viewer explicitly created by the user's macro.
  if (gui && commandStatus == 0 && !visManager->GetCurrentViewer()) {
    ui->ApplyCommand("/run/initialize");
    const auto viewerStatus = ui->ApplyCommand("/vis/open OGL");
    if (viewerStatus == 0) {
      ui->ApplyCommand("/vis/drawVolume voxel_block");
      ui->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 90 0");
      const auto& grid = detector->GetGrid();
      // Show the world origin on the grid entrance face.
      std::ostringstream axes;
      axes << "/vis/scene/add/axes "
           << "0 0 0 "
           << grid.GetSizeX() / um << " um";
      ui->ApplyCommand(axes.str());
      ui->ApplyCommand("/vis/viewer/zoom 1.3");
      ui->ApplyCommand("/vis/scene/add/trajectories smooth");
      ui->ApplyCommand("/vis/viewer/flush");
    } else {
      std::cerr << "Could not open the OpenGL viewer (status " << viewerStatus << ").\n";
      commandStatus = viewerStatus;
    }
  }
#endif

  if (uiExec) {
    uiExec->SessionStart();
    delete uiExec;
  } else if (macro.empty()) {
    PrintUsage(argv[0]);
  }

#ifdef TRACKINGGRID_WITH_GEANT4_UIVIS
  delete visManager;
#endif
  delete runManager;
  delete trackingOptions;
  return commandStatus == 0 ? 0 : 1;
}
