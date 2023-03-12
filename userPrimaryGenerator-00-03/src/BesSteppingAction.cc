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
/// \file BesSteppingAction.cc
/// \brief Implementation of the BesSteppingAction class

#include "BesSteppingAction.hh"
#include "BesEventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include <fstream>
using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BesSteppingAction::BesSteppingAction(BesEventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0)
{
m_px.clear();
m_py.clear();
m_pz.clear();
m_e.clear();
m_x.clear();
m_y.clear();
m_z.clear();
m_t.clear();
m_particleID.clear();
m_AnSec.clear();
trk1flag = false;
trk2flag = false;
trk3flag = false;
trk4flag = false;
m_eventID = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BesSteppingAction::~BesSteppingAction()
{
myfile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BesSteppingAction::UserSteppingAction(const G4Step* step)
{
		G4int parID = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
		G4double posi_x = step->GetPreStepPoint()->GetPosition().x();
		G4double posi_y = step->GetPreStepPoint()->GetPosition().y();
		G4double posi_z = step->GetPreStepPoint()->GetPosition().z();
		G4double posi_t = step->GetPreStepPoint()->GetGlobalTime();
		G4double mom_x = step->GetPreStepPoint()->GetMomentum().x();
		G4double mom_y = step->GetPreStepPoint()->GetMomentum().y();
		G4double mom_z = step->GetPreStepPoint()->GetMomentum().z();
		G4double mom_e = sqrt(step->GetPreStepPoint()->GetMass() * step->GetPreStepPoint()->GetMass() + step->GetPreStepPoint()->GetMomentum().mag() * step->GetPreStepPoint()->GetMomentum().mag());

		if(fEventAction->GetEventID() != m_eventID){
				trk1flag = false;
				trk2flag = false;
				trk3flag = false;
				trk4flag = false;
				trk5flag = false;
				writeflag = false;
				m_particleID.clear();
				m_px.clear();
				m_py.clear();
				m_pz.clear();
				m_e.clear();
				m_x.clear();
				m_y.clear();
				m_z.clear();
				m_t.clear();
				m_AnSec.clear();
				m_eventID = fEventAction->GetEventID();
		}
//		if(step->GetTrack()->GetTrackID() == 1 || step->GetTrack()->GetTrackID() == 2 || step->GetTrack()->GetTrackID() == 3 || step->GetTrack()->GetParentID() == 2   ){
//				cout << step->GetTrack()->GetParticleDefinition()->GetParticleName() << endl;
//				cout << step->GetTrack()->GetPosition() << endl;
//		}

		if(step->GetTrack()->GetTrackID() == 1 ){
				if(!trk1flag){
						trk1flag = true;
						//		cout << parID << " " << mom_x << " " << mom_y << " " << mom_z << " " << mom_e << " " << " " << posi_x << " " << posi_y << " " << posi_z << " " << posi_t << endl;
						m_particleID.push_back(parID);
						m_px.push_back(mom_x);
						m_py.push_back(mom_y);
						m_pz.push_back(mom_z);
						m_e.push_back(mom_e);
						m_x.push_back(posi_x);
						m_y.push_back(posi_y);
						m_z.push_back(posi_z);
						m_t.push_back(posi_t);
				}
		}
		else if(step->GetTrack()->GetTrackID() == 3){
				if(!trk2flag){
						trk2flag = true;
						m_particleID.push_back(parID);
						m_px.push_back(mom_x);
						m_py.push_back(mom_y);
						m_pz.push_back(mom_z);
						m_e.push_back(mom_e);
						m_x.push_back(posi_x);
						m_y.push_back(posi_y);
						m_z.push_back(posi_z);
						m_t.push_back(posi_t);
				}
		}
		else if(step->GetTrack()->GetTrackID() == 2){
				if(!trk5flag){
						trk5flag = true;
						m_particleID.push_back(parID);
						m_px.push_back(mom_x);
						m_py.push_back(mom_y);
						m_pz.push_back(mom_z);
						m_e.push_back(mom_e);
						m_x.push_back(posi_x);
						m_y.push_back(posi_y);
						m_z.push_back(posi_z);
						m_t.push_back(posi_t);
				}
		}


		else if(step->GetTrack()->GetParentID() == 2){
				if(parID == 2212 ){
						if(!trk3flag){
								trk3flag = true;
								m_particleID.push_back(parID);
								m_px.push_back(mom_x);
								m_py.push_back(mom_y);
								m_pz.push_back(mom_z);
								m_e.push_back(mom_e);
								m_x.push_back(posi_x);
								m_y.push_back(posi_y);
								m_z.push_back(posi_z);
								m_t.push_back(posi_t);
								m_AnSec.push_back(parID);
						}
				}
				else if(parID == -2212){
						if(!trk4flag){
								trk4flag = true;
								m_particleID.push_back(parID);
								m_px.push_back(mom_x);
								m_py.push_back(mom_y);
								m_pz.push_back(mom_z);
								m_e.push_back(mom_e);
								m_x.push_back(posi_x);
								m_y.push_back(posi_y);
								m_z.push_back(posi_z);
								m_t.push_back(posi_t);
								m_AnSec.push_back(parID);
						}
				}
		}


		if(m_px.size() == 5 && m_AnSec.size() == 2 && !writeflag && trk4flag && trk3flag && trk2flag  && trk1flag){
				myfile << 5 << "\n";
				for(int i = 0; i < 5; i++){
						myfile << m_particleID[i] << "  "
								<< m_px[i] << "  "
								<< m_py[i] << "  "
								<< m_pz[i] << "  "
								<< m_e[i] << "  "
								<< m_x[i] << "  "
								<< m_y[i] << "  "
								<< m_z[i] << "  "
								<< m_t[i] << "\n";
				}
				writeflag = true;
		}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BesSteppingAction::SetOutFile(G4String str){
		myfile.open(str);
}

