#ifndef TRACKINGGRID_DETECTOR_MESSENGER_HH
#define TRACKINGGRID_DETECTOR_MESSENGER_HH

#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIcmdWithABool;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcommand;
class G4UIcmdWith3VectorAndUnit;

class DetectorMessenger : public G4UImessenger {
public:
  explicit DetectorMessenger(DetectorConstruction* detector);
  ~DetectorMessenger() override;

  void SetNewValue(G4UIcommand* command, G4String value) override;

private:
  DetectorConstruction* fDetector;
  G4UIdirectory* fDirectory;
  G4UIcmdWithADoubleAndUnit* fVoxelSizeCommand;
  G4UIcmdWithADoubleAndUnit* fMaxStepCommand;
  G4UIcmdWithABool* fCheckOverlapsCommand;
  G4UIcommand* fGridCommand;
  G4UIcmdWith3VectorAndUnit* fGridCentreCommand;
  G4UIcmdWithADoubleAndUnit* fWorldCommand;
};

#endif
