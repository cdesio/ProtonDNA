//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"
#include "CommandLineParser.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4AnalysisManager.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

#include <map>
#include <set>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>



using namespace G4DNAPARSER;
using CLHEP::nanometer;

namespace
{
  G4Mutex messangerInit = G4MUTEX_INITIALIZER;

  // Track bad particle IDs globally (used for logging the warnings at end of run)
  std::set<int> g_badParticleIDs;
  G4Mutex g_badIDsMutex = G4MUTEX_INITIALIZER;
  G4int g_nBadPrimaryWarnings = 0; 
  G4Mutex g_warningCountMutex = G4MUTEX_INITIALIZER;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(G4String PS_data)
    : G4VUserPrimaryGeneratorAction(), fpParticleGun(nullptr), fParticleGun(nullptr), fPS_data(PS_data)
{
  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0); // does this do anything?
  if ((parser->GetCommandIfActive("-ps"))||(parser->GetCommandIfActive("-photonPS")))
  {
    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);
    fParticleGun->SetParticleEnergy(0);
    fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));
  }
  else
  {
    fpParticleGun = new G4GeneralParticleSource();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  // Log all bad particle IDs when the run is complete
  if (!g_badParticleIDs.empty()) {
    G4cout << "\n========================================" << G4endl;
    G4cout << "SUMMARY: Bad Particle IDs Encountered" << G4endl;
    G4cout << "========================================" << G4endl;
    G4cout << "Total unique unknown particle IDs: " << g_badParticleIDs.size() << G4endl;
    G4cout << "Unknown IDs: ";
    for (auto id : g_badParticleIDs) {
      G4cout << id << " ";
    }
    G4cout << "\n========================================\n" << G4endl;
  }else {

    G4cout << "Great no Bad Particle IDs " << G4endl;

  }
  
  // handle cleanup of particle guns
  if (fParticleGun) delete fParticleGun;
  if (fpParticleGun) delete fpParticleGun;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0); // again, does this do anything?
  if (parser->GetCommandIfActive("-ps"))
  {
    G4long eventNum = anEvent->GetEventID();

    std::ifstream ps_file (fPS_data, std::ifstream::binary);

    constexpr size_t recordDoubles = 18;
    ps_file.seekg(eventNum * static_cast<std::streamoff>(recordDoubles * sizeof(double)), ps_file.beg);
    std::vector<double> line(recordDoubles);
    if (!ps_file.read(reinterpret_cast<char*>(line.data()), recordDoubles * sizeof(double))) {
      G4Exception("PrimaryGeneratorAction::GeneratePrimaries", "InvalidPhaseSpace", FatalException,
                  "Could not read a complete 18-double phase-space record.");
      return;
    }
    G4double positionX = line[0];
    G4double positionY = line[1];
    G4double positionZ = line[2];
    G4double momentumX = line[3];
    G4double momentumY = line[4];
    G4double momentumZ = line[5];
    G4double particleEnergy = line[6];
    part1_EventNum = line[7];
    primaryParticle = line[8];
    part1_CopyNum = line[9];
    part1_Time = line[10];
    part1_particleSource = line[11];

    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    // G4ParticleDefinition *particle;
    G4ParticleDefinition *particle = nullptr;  // <-- مهم جداً

  
    const int particleID = static_cast<int>(primaryParticle);
    particle = particleTable->FindParticle(particleID);
    if (!particle && particleID > 1000000000) {
      const auto Z = (particleID / 10000) % 1000;
      const auto A = (particleID / 10) % 1000;
      particle = G4IonTable::GetIonTable()->GetIon(Z, A, 0);
    }
    if (particle) primaryName = particle->GetParticleName();
    // if particle is still null (shouldn't happen if tables are initialised) handle fatal error safely
    if (!particle) {
      G4Exception("PrimaryGeneratorAction::GeneratePrimaries",
              "ParticleCreationFailed",
              FatalException,
              ("Failed to create particle: " + primaryName + 
               " (ID: " + std::to_string(particleID) + ")").c_str());
      return;
    }
    
    // G4cout << primaryName << " position = " << G4ThreeVector(positionX, positionY, positionZ) << " momentum = " << G4ThreeVector(momentumX, momentumY, momentumZ) << " energy = " << particleEnergy << "Ra event number = " << part1_EventNum << " copy number = " << part1_CopyNum << " source particle = " << part1_particleSource <<G4endl;
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticlePosition(G4ThreeVector(positionX, positionY, positionZ));
    fParticleGun->SetParticleEnergy(particleEnergy);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(momentumX, momentumY, momentumZ));

    fParticleGun->GeneratePrimaryVertex(anEvent);

  }
  else if (parser->GetCommandIfActive("-photonPS"))
  {
    G4long eventNum = anEvent->GetEventID();


    std::ifstream ps_file (fPS_data, std::ifstream::binary);

    ps_file.seekg(eventNum*8*8, ps_file.beg); // find position of event data. 8 doubles which are each 8 bytes

    double line[8];

    ps_file.read((char *)&line, sizeof line);

    ps_file.close();

    G4double positionX = line[0];
    G4double positionY = line[1];
    G4double positionZ = line[2];
    G4double momentumX = line[3];
    G4double momentumY = line[4];
    G4double momentumZ = line[5];
    G4double particleEnergy = line[6];
    part1_EventNum = line[7];

    // G4cout <<  " position = " << G4ThreeVector(positionX, positionY, positionZ) << " momentum = " << G4ThreeVector(momentumX, momentumY, momentumZ) << " energy = " << particleEnergy << " Co60 event number = " << part1_EventNum <<G4endl;


    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition *particle = particleTable->FindParticle("e-");
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticlePosition(G4ThreeVector(positionX, positionY, positionZ));
    fParticleGun->SetParticleEnergy(particleEnergy);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(momentumX, momentumY, momentumZ));

    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
  else
  {
    fpParticleGun->GeneratePrimaryVertex(anEvent);
    primaryName = fpParticleGun->GetCurrentSource()->GetParticleDefinition()->GetParticleName();
  }
}
