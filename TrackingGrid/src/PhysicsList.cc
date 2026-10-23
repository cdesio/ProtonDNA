#include "PhysicsList.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList(G4bool transportOnly) : fTransportOnly(transportOnly)
{
  SetVerboseLevel(0);
  SetDefaultCutValue(1 * um);
  RegisterPhysics(new G4EmStandardPhysics_option4());
  RegisterPhysics(new G4DecayPhysics());
  RegisterPhysics(new G4RadioactiveDecayPhysics());
  RegisterPhysics(new G4HadronElasticPhysics());
  RegisterPhysics(new G4HadronPhysicsFTFP_BERT());
  RegisterPhysics(new G4StoppingPhysics());
  RegisterPhysics(new G4IonPhysics());
  auto* limiter = new G4StepLimiterPhysics();
  limiter->SetApplyToAll(true);
  RegisterPhysics(limiter);
}
void PhysicsList::ConstructParticle() { G4VModularPhysicsList::ConstructParticle(); }
void PhysicsList::ConstructProcess()
{
  if (fTransportOnly) AddTransportation();
  else G4VModularPhysicsList::ConstructProcess();
}
