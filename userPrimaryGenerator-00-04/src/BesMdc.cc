/*****************************************************************************\
 *  BOSS (BES Offline Software System)                                       *
 *       Filename:  BesMdc.cc
 *                                                                           *
 *    Description:                                                           
 *                                                                           *
 *        Version:  1.0                                                      *
 *        Created:  11/08/2022 03:06:35 PM                                   *
 *       Revision:  none                                                     *
 *       Compiler:  gcc                                                      *
 *                                                                           *
 *         Author:  Liang Liu (USTC), <liangzy@mail.ustc.edu.cn>             *
 *   Organization:  BESIII Collaboration                                     *
 *        License:  GNU Lesser General Public License 3.0, see LICENSE.md.   *
\*****************************************************************************/
#include "BesMdc.hh"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4UnionSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include "BesSubdetector.hh"
using namespace std;
using namespace CLHEP;


BesMdc::BesMdc(){
}


BesMdc::~BesMdc(){

}


void BesMdc::DefineMaterial(){
		G4double density, a, z, fractionmass;
		G4int nel, natoms;
		G4String name, symbol;

		Al = G4Material::GetMaterial("Aluminium");
		G4Element *C  =G4Element::GetElement("Carbon");
		G4Element *H  =G4Element::GetElement("Hydrogen");
		G4Element *O  =G4Element::GetElement("Oxygen");

		density=1.57*g/cm3;
		G4int nElement=3;
		CarbonFiber=new G4Material("CarbonFiber",density,nElement);
		CarbonFiber->AddElement(C,0.697);
		CarbonFiber->AddElement(H,0.0061);
		CarbonFiber->AddElement(O,0.2969);


		a = 14.01*g/mole;
		G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);


		density = 1.290*mg/cm3;
		Air = new G4Material(name="Air",density,nel=2);
		Air->AddElement(N, fractionmass = 0.7);
		Air->AddElement(O, fractionmass = 0.3);




}


void BesMdc::Construct(G4LogicalVolume* logicbes){
		DefineMaterial();

		// the logical volume of MDC inner wall
		//
/*
		G4Tubs* solidMdc = new G4Tubs("solidMdc",62.9, 64.25, 250,0,360);
		logicalMdc = new G4LogicalVolume(solidMdc, G4Material::GetMaterial("Air"),"logicalMdc");
		physicalMdc = new G4PVPlacement(0,G4ThreeVector(0,0,0),logicalMdc,"physicalMdc",logicbes,false,0);

		G4Tubs* solidinnerAl = new G4Tubs("solidinnerAl",62.9, 63.0, 250,0,360);
		logicalinnerAl = new G4LogicalVolume(solidinnerAl, G4Material::GetMaterial("Aluminium"), "logicalinnerAl");
		physicalinnerAl = new G4PVPlacement(0, G4ThreeVector(0,0,0), logicalinnerAl, "physicalinnerAl", logicalMdc, false, 0);


		G4Tubs* solidCarbonFiber = new G4Tubs("solidCarbonFiber",63.0, 64.2, 250,0,360);
		logicalCarbonFiber = new G4LogicalVolume(solidCarbonFiber, G4Material::GetMaterial("CarbonFiber"), "logicalCarbonFiber");
		physicalCarbonFiber = new G4PVPlacement(0, G4ThreeVector(0,0,0), logicalCarbonFiber, "physicalCarbonFiber", logicalMdc, false, 0);


		G4Tubs* solidouterAl = new G4Tubs("solidouterAl",64.2, 64.25, 250,0,360);
		logicalouterAl = new G4LogicalVolume(solidouterAl, G4Material::GetMaterial("Aluminium"), "logicalouterAl");
		physicalouterAl = new G4PVPlacement(0, G4ThreeVector(0,0,0), logicalouterAl, "physicalouterAl", logicalMdc, false, 0);


	
		G4VisAttributes *visMdc= new G4VisAttributes(G4Colour(0.,0.,1.));
		logicalMdc->SetVisAttributes(visMdc);

		G4VisAttributes *visAl= new G4VisAttributes(G4Colour(0.,1.,0.));
		logicalinnerAl->SetVisAttributes(visAl);
		logicalouterAl->SetVisAttributes(visAl);

		G4VisAttributes *visCarbonFiber= new G4VisAttributes(G4Colour(1.,0.,0.));
		logicalCarbonFiber->SetVisAttributes(visCarbonFiber);
*/

		G4Tubs* solidMdc = new G4Tubs("solidMdc",63.0, 64.2, 440,0,360);
		logicalMdc = new G4LogicalVolume(solidMdc, G4Material::GetMaterial("Air"),"logicalMdc");
		physicalMdc = new G4PVPlacement(0,G4ThreeVector(0,0,0),logicalMdc,"physicalMdc",logicbes,false,0);

		G4Tubs* solidCarbonFiber = new G4Tubs("solidCarbonFiber",63.0, 64.2, 440,0,360);
		logicalCarbonFiber = new G4LogicalVolume(solidCarbonFiber, G4Material::GetMaterial("CarbonFiber"), "logicalCarbonFiber");
		physicalCarbonFiber = new G4PVPlacement(0, G4ThreeVector(0,0,0), logicalCarbonFiber, "physicalCarbonFiber", logicalMdc, false, 0);

	
		G4VisAttributes *visMdc= new G4VisAttributes(G4Colour(0.,0.,1.));
		logicalMdc->SetVisAttributes(visMdc);

		G4VisAttributes *visCarbonFiber= new G4VisAttributes(G4Colour(1.,0.,0.));
		logicalCarbonFiber->SetVisAttributes(visCarbonFiber);


}
