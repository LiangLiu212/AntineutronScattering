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
/// \file BesSteppingAction.hh
/// \brief Definition of the BesSteppingAction class

#ifndef BesSteppingAction_h
#define BesSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

class BesEventAction;

class G4LogicalVolume;

/// Stepping action class
/// 

class BesSteppingAction : public G4UserSteppingAction
{
  public:
    BesSteppingAction(BesEventAction* eventAction);
    virtual ~BesSteppingAction();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);
	virtual void SetOutFile(G4String str);

  private:
    BesEventAction*  fEventAction;
    G4LogicalVolume* fScoringVolume;
	std::vector<int> m_particleID;
	std::vector<int> m_AnSec;
	std::vector<double> m_px;
	std::vector<double> m_py;
	std::vector<double> m_pz;
	std::vector<double> m_e;
	std::vector<double> m_x;
	std::vector<double> m_y;
	std::vector<double> m_z;
	std::vector<double> m_t;
	ofstream myfile;

	bool trk1flag;
	bool trk2flag;
	bool trk3flag;
	bool trk4flag;
	bool trk5flag;
	bool writeflag;
	int m_eventID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
