#ifndef TRACKINGGRID_PHYSICS_LIST_HH
#define TRACKINGGRID_PHYSICS_LIST_HH
#include "G4VModularPhysicsList.hh"
class PhysicsList : public G4VModularPhysicsList {
public:
  explicit PhysicsList(G4bool transportOnly = false);
  void ConstructParticle() override;
  void ConstructProcess() override;
private:
  G4bool fTransportOnly;
};
#endif
