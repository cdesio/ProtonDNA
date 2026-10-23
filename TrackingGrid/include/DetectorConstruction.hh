#ifndef TRACKINGGRID_DETECTOR_CONSTRUCTION_HH
#define TRACKINGGRID_DETECTOR_CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "VoxelGrid.hh"

class DetectorMessenger;
class G4LogicalVolume;
class G4VPhysicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction();
  ~DetectorConstruction() override;

  G4VPhysicalVolume* Construct() override;

  void SetGridCentre(const G4ThreeVector& value) { fGridCentre = value; fCustomGridCentre = true; }
  G4ThreeVector GetGridCentre() const {
    return fCustomGridCentre ? fGridCentre : G4ThreeVector(0, 0, fGrid.GetSizeZ()/2);
  }
  void SetWorldHalfSize(G4double value);
  void SetVoxelSize(G4double value);
  void SetGrid(G4int nx, G4int ny, G4int nz);
  void SetMaxStep(G4double value);
  void SetCheckOverlaps(G4bool value);

  const VoxelGrid& GetGrid() const { return fGrid; }
  G4double GetMaxStep() const { return fMaxStep; }
  G4bool GetCheckOverlaps() const { return fCheckOverlaps; }

private:
  VoxelGrid fGrid;
  G4ThreeVector fGridCentre;
  G4bool fCustomGridCentre = false;
  G4double fWorldHalfSize;
  G4double fMaxStep;
  G4bool fCheckOverlaps;
  DetectorMessenger* fMessenger;
  G4LogicalVolume* fVoxelLogical;
};

#endif
