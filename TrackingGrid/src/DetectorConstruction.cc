#include "DetectorConstruction.hh"

#include "DetectorMessenger.hh"
#include "VoxelParameterisation.hh"

#include "G4Box.hh"
#include "G4Exception.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include <cmath>

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),
    fGridCentre(0, 0, 0),
    fWorldHalfSize(10 * cm),
    fMaxStep(30.0 * nm),
    fCheckOverlaps(false),
    fMessenger(new DetectorMessenger(this)),
    fVoxelLogical(nullptr)
{}

DetectorConstruction::~DetectorConstruction()
{
  delete fMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  auto* nist = G4NistManager::Instance();
  auto* water = nist->FindOrBuildMaterial("G4_WATER");

  const auto centre = GetGridCentre();
  if (std::abs(centre.x()) + fGrid.GetSizeX()/2 >= fWorldHalfSize ||
      std::abs(centre.y()) + fGrid.GetSizeY()/2 >= fWorldHalfSize ||
      std::abs(centre.z()) + fGrid.GetSizeZ()/2 >= fWorldHalfSize) {
    G4Exception("DetectorConstruction::Construct", "GeometryBounds", FatalException,
                "Grid must fit strictly inside the world.");
  }
  // Vacuum isolates the voxel transport from upstream scattering.
  auto* worldMaterial = nist->FindOrBuildMaterial("G4_Galactic");
  auto* worldSolid = new G4Box("world", fWorldHalfSize, fWorldHalfSize, fWorldHalfSize);
  auto* worldLogical = new G4LogicalVolume(worldSolid, worldMaterial, "world");
  auto* worldPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), worldLogical,
                                          "world", nullptr, false, 0, fCheckOverlaps);
  // Optional water absorber: uncomment and choose its size/position when needed.
  // Dimensions are full lengths. Keep it clear of the voxel grid and world edges.
  /*
  const G4ThreeVector absorberSize(5 * cm, 5 * cm, 1 * cm);
  const G4ThreeVector absorberCentre(0, 0, -2 * cm);
  auto* absorberSolid = new G4Box("absorber", absorberSize.x()/2,
                                 absorberSize.y()/2, absorberSize.z()/2);
  auto* absorberLogical = new G4LogicalVolume(absorberSolid, water, "absorber");
  new G4PVPlacement(nullptr, absorberCentre, absorberLogical, "absorber",
                    worldLogical, false, 0, fCheckOverlaps);
  absorberLogical->SetVisAttributes(new G4VisAttributes(G4Colour(0, 0.7, 0.9, 0.25)));
  */

  auto* blockSolid = new G4Box("voxel_block",
                               0.5 * fGrid.GetSizeX(),
                               0.5 * fGrid.GetSizeY(),
                               0.5 * fGrid.GetSizeZ());
  auto* blockLogical = new G4LogicalVolume(blockSolid, water, "voxel_block");
  new G4PVPlacement(nullptr, centre, blockLogical, "voxel_block",
                    worldLogical, false, 0, fCheckOverlaps);


  auto* voxelSolid = new G4Box("voxel",
                               0.5 * fGrid.GetVoxelSize(),
                               0.5 * fGrid.GetVoxelSize(),
                               0.5 * fGrid.GetVoxelSize());
  fVoxelLogical = new G4LogicalVolume(voxelSolid, water, "voxel");
  if (fMaxStep > 0.0) {
    fVoxelLogical->SetUserLimits(new G4UserLimits(fMaxStep));
  }

  auto* parameterisation = new VoxelParameterisation(&fGrid);
  new G4PVParameterised("voxel", fVoxelLogical, blockLogical, kUndefined,
                        fGrid.GetNumberOfVoxels(), parameterisation, fCheckOverlaps);

  auto* invis = new G4VisAttributes(false);
  worldLogical->SetVisAttributes(invis);
  blockLogical->SetVisAttributes(invis);

  auto* voxelVis = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.18));
  voxelVis->SetForceSolid(true);
  fVoxelLogical->SetVisAttributes(voxelVis);

  G4cout << "Voxel grid: " << fGrid.GetNx() << " x "
         << fGrid.GetNy() << " x " << fGrid.GetNz()
         << ", voxel size " << fGrid.GetVoxelSize() / nm << " nm"
         << ", max step " << fMaxStep / nm << " nm"
         << ", check overlaps " << (fCheckOverlaps ? "true" : "false")
         << ", total depth " << fGrid.GetSizeZ() / um << " um" << G4endl;

  return worldPhysical;
}

void DetectorConstruction::SetVoxelSize(G4double value)
{
  fGrid.SetVoxelSize(value);
  if (G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

void DetectorConstruction::SetGrid(G4int nx, G4int ny, G4int nz)
{
  fGrid.SetGrid(nx, ny, nz);
  if (G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

void DetectorConstruction::SetMaxStep(G4double value)
{
  fMaxStep = value;
  if (G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

void DetectorConstruction::SetCheckOverlaps(G4bool value)
{
  fCheckOverlaps = value;
  if (G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

void DetectorConstruction::SetWorldHalfSize(G4double value) { fWorldHalfSize = value; }
