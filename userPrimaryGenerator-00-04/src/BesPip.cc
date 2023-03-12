//---------------------------------------------------------------------------//
//      BOOST --- BESIII Object_Oreiented Simulation Tool                    //
//---------------------------------------------------------------------------//
//Descpirtion: BES beam pipe 
//Author: Liuhm
//Created: May 21, 2003
//Comment:
//---------------------------------------------------------------------------//
//
#include "BesPip.hh"

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
#include "BesPipParameter.hh"
#include "BesSubdetector.hh"
using namespace std;
using namespace CLHEP;
BesPip::BesPip()
{
  pipPar = new BesPipParameter();
  pipPar->ReadData();
  for(G4int i=0; i<3; i++)
  {
    goldLayer[i] = pipPar->getGoldLayer(i);
    innerBe[i] = pipPar->getInnerBe(i);
    oilLayer[i] = pipPar->getOilLayer(i);
    outerBe[i] = pipPar->getOuterBe(i);
  }

  for(G4int i=0; i<5; i++)
  {
    innerBeSide[i] = pipPar->getInnerBeSide(i);
    innerAl[i] = pipPar->getInnerAl(i);
    AlRing[i] = pipPar->getAlRing(i);
    outerAl[i] = pipPar->getOuterAl(i);
    AgLayer[i] = pipPar->getAgLayer(i);
    CuLayer[i] = pipPar->getCuLayer(i);
    AlCover[i] = pipPar->getAlCover(i);
  }

  logicalPip = 0;
  physicalPip = 0;

  logicalgoldLayer = 0;
  physicalgoldLayer = 0;

  logicalinnerBe = 0;
  physicalinnerBe = 0;

  logicaloilLayer = 0;
  physicaloilLayer = 0;

  logicalouterBe = 0;
  physicalouterBe = 0;

   logicalinnerBeSide = 0;
   physicalinnerBeSide = 0;

   logicalinnerAl = 0;
   physicalinnerAl = 0;

   logicalAlRing = 0;
   physicalAlRing = 0;

   logicalouterAl = 0;
   physicalouterAl = 0;
  
   logicalAgLayer = 0;
   physicalAgLayer = 0;

   logicalCuLayer = 0;
   physicalCuLayer = 0;

   logicalAlCover = 0;
   physicalAlCover = 0;

  Au = 0;
  Ag = 0;
  Oil = 0;
}

void BesPip::DefineMaterial()
{
  G4double density, a, z,fractionmass;
  G4int nel,natoms;
  G4String name, symbol;

  density=19.32*CLHEP::g/CLHEP::cm3;
  a = 196.967*CLHEP::g/CLHEP::mole;
  Au= new G4Material("Gold",79,a,density);

  density=10.5*CLHEP::g/CLHEP::cm3;
  a = 107.9*CLHEP::g/CLHEP::mole;
  Ag= new G4Material("Silver",47,a,density);

  a = 1.01*CLHEP::g/CLHEP::mole;
  G4Element* H  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);
    
  a = 12.01*CLHEP::g/CLHEP::mole;
  G4Element* C  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

  density = 0.810*CLHEP::g/CLHEP::cm3;
  Oil = new G4Material(name="Oil",density,nel=2);
  Oil->AddElement(C, natoms=18);
  Oil->AddElement(H, natoms=38);

}

void BesPip::Construct(G4LogicalVolume* logicalbes)
{
  DefineMaterial();
cout << "Ok 1"<< endl;
    G4double a, z, density,abundance,fractionmass;
    G4double temperature, pressure;
    G4int iz,n,nel,ncomponents,natoms;
    G4String name,symbol;
    G4Isotope* U5 = new G4Isotope(name="U235", iz=92, n=235, a=235.01*g/mole);
    G4Isotope* U8 = new G4Isotope(name="U238", iz=92, n=238, a=238.03*g/mole);

    G4Element* U  = new G4Element(name="enriched Uranium",symbol="U",ncomponents=2);
    U->AddIsotope(U5, abundance= 90.*perCent);
    U->AddIsotope(U8, abundance= 10.*perCent);
    
    a = 1.01*g/mole;
    G4Element* H  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);
    
    a = 12.01*g/mole;
    G4Element* C  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

    a = 14.01*g/mole;
    G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

    a = 16.00*g/mole;
    G4Element* O  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

    a=28.09*g/mole;
    G4Element* Si = new G4Element(name="Silicon",symbol="Si",z=14.,a);
    
    a=126.90447*g/mole;
    G4Element* I = new G4Element(name="Iodine",symbol="I",z=53.,a);

    a=132.90545*g/mole;
    G4Element* Cs = new G4Element(name="Cesium",symbol="Cs",z=55.,a);

    density = 1.848*g/cm3;
    a = 9.012182*g/mole;
    G4Material* Beryllium = new G4Material(name="Beryllium",z=4.0,a,density);

    density = 2.70*g/cm3;
    a = 26.98*g/mole;
    G4Material* Aluminium = new G4Material(name="Aluminium",z=13.,a,density);

    density = 7.84*g/cm3;
    a = 55.845*g/mole;
    G4Material* Iron = new G4Material(name="Iron",z=26.0,a,density);

    density = 8.96*g/cm3;
    a = 63.546*g/mole;
    G4Material* Copper = new G4Material(name="Copper",z=29.0,a,density);

    density = 0.0001664*g/cm3;
    a = 4.0026*g/mole;
    G4Material* Hegas = new G4Material(name="Hegas",z=2.0,a,density);

    density = 0.001879*g/cm3;
    G4Material* Propane = new G4Material(name="Propane",density,nel=2);
    Propane->AddElement(C, natoms=3);
    Propane->AddElement(H, natoms=8);

    density = 4.53*g/cm3;
    G4Material* Cesiumiodide = new G4Material(name="Cesiumiodide",density,nel=2);
    Cesiumiodide->AddElement(I, natoms=1);
    Cesiumiodide->AddElement(Cs, natoms=1);

    density = 0.00085144*g/cm3;
    G4Material* Mdcgas = new G4Material(name="Mdcgas",density,nel=2);
    Mdcgas->AddMaterial(Hegas, fractionmass = 0.1173);
    Mdcgas->AddMaterial(Propane, fractionmass = 0.8827);

    density = 1.390*g/cm3;
    a = 39.95*g/mole;
    G4Material* liquidArgon = new G4Material(name="liquidArgon",z=18.0,a,density);

    density = 11.35*g/cm3; 
    a = 207.19*g/mole; 
    G4Material* Lead = new G4Material(name="Lead",z=82.,a,density);

    density = 1.0*g/cm3;
    G4Material* Water = new G4Material(name="Water", density, nel=2);
    Water->AddElement(H, natoms=2);
    Water->AddElement(O, natoms=1);

    density = 1.032*g/cm3;
    G4Material* Scintillator = new G4Material(name="Scintillator", density, nel=2);
    Scintillator->AddElement(C, natoms=9);
    Scintillator->AddElement(H, natoms=10);  

    density = 2.200*g/cm3;
    G4Material* SiO2 = new G4Material(name="SiO2", density, nel=2);
    SiO2->AddElement(Si, natoms=1);
    SiO2->AddElement(O, natoms=2);

    density = 1.290*mg/cm3;
    G4Material* Air = new G4Material(name="Air",density,nel=2);
    Air->AddElement(N, fractionmass = 0.7);
    Air->AddElement(O, fractionmass = 0.3);
    
    density = 0.200*g/cm3;
    G4Material* Aerogel = new G4Material(name="Aerogel",density,nel=3);
    Aerogel->AddMaterial(SiO2, fractionmass = 0.625);
    Aerogel->AddMaterial(Water, fractionmass = 0.374);
    Aerogel->AddElement(C, fractionmass = 0.001);
    
    density     = 27.0*mg/cm3;
    pressure    = 50.0*atmosphere;
    temperature = 325.0*kelvin;
    G4Material* CarbonicGas = new G4Material(name="CarbonicGas",density,nel=2,kStateGas,temperature,pressure);
    CarbonicGas->AddElement(C,natoms=1);
    CarbonicGas->AddElement(O,natoms=2);

    density     = 0.3*mg/cm3;
    pressure    = 2.0*atmosphere;
    temperature = 500.0*kelvin;
    G4Material* WaterSteam = new G4Material(name="WaterSteam",density,nel=1,kStateGas,temperature,pressure);
    WaterSteam->AddMaterial(Water,fractionmass=1);

    density     = universe_mean_density;
    pressure    = 3.e-18*pascal;
    temperature = 2.73*kelvin;
    G4Material* Galactic = new G4Material(name="Galactic", z=1., a=1.01*g/mole,
                                       density,kStateGas,temperature,pressure);

    density     = 1.0e-5*g/cm3;
    pressure    = 2.e-2*bar;
    G4Material* Beam = new G4Material(name="Beam",density,nel=1,kStateGas,STP_Temperature,pressure);
    Beam->AddMaterial(Air,fractionmass=1.0);



   //G4RotationMatrix* xRot = new G4RotationMatrix;
  //xRot->rotateX(90*deg);

  //the logical volume of beam pipe
  G4Tubs* solidPip1 = new G4Tubs("solidPip1",0.,33.7,134,0,360);
  G4Tubs* solidPip2 = new G4Tubs("solidPip2",0.,48,66,0,360);
  G4UnionSolid* solidPip_tmp = new G4UnionSolid("solidPip_tmp",solidPip1,solidPip2,0,G4ThreeVector(0,0,-167));
  G4UnionSolid* solidPip = new G4UnionSolid("solidPip",solidPip_tmp,solidPip2,0,G4ThreeVector(0,0,167));
  logicalPip = new G4LogicalVolume(solidPip, G4Material::GetMaterial("Beam"),"logicalPip");
  physicalPip = new G4PVPlacement(0,G4ThreeVector(0,0,0),logicalPip,"physicalPip",logicalbes,false,0);

cout << "Ok 2"<< endl;
  //the volume of gold layer
  G4Tubs* solidgoldLayer = new G4Tubs("solidgoldLayer",goldLayer[0],goldLayer[1],goldLayer[2]/2,0,360);
  logicalgoldLayer = new G4LogicalVolume(solidgoldLayer, Au,"logicalgoldLayer");
  physicalgoldLayer = new G4PVPlacement(0,G4ThreeVector(0,0,0),logicalgoldLayer,"physicalgoldLayer",logicalPip,false,0);

cout << "Ok 3"<< endl;
  //the volume of inner Be pipe
  G4Tubs* solidinnerBe = new G4Tubs("solidinnerBe",innerBe[0],innerBe[1],innerBe[2]/2,0,360);
  logicalinnerBe = new G4LogicalVolume(solidinnerBe, G4Material::GetMaterial("Beryllium"),"logicalinnerBe");
  physicalinnerBe = new G4PVPlacement(0,G4ThreeVector(0,0,0),logicalinnerBe,"physicalinnerBe",logicalPip,false,0);

cout << "Ok 4"<< endl;
  //the volume of oil layer
  G4Tubs* solidoilLayer = new G4Tubs("solidoilLayer",oilLayer[0],oilLayer[1],oilLayer[2]/2,0,360);
  logicaloilLayer = new G4LogicalVolume(solidoilLayer, Oil,"logicaloilLayer");
  physicaloilLayer = new G4PVPlacement(0,G4ThreeVector(0,0,0),logicaloilLayer,"physicaloilLayer",logicalPip,false,0);

cout << "Ok 5" << endl;
  //the volume of outer Be pipe
   G4Tubs* solidouterBe = new G4Tubs("solidouterBe",outerBe[0],outerBe[1],outerBe[2]/2,0,360);
  logicalouterBe = new G4LogicalVolume(solidouterBe, G4Material::GetMaterial("Beryllium"),"logicalouterBe");
  physicalouterBe = new G4PVPlacement(0,G4ThreeVector(0,0,0),logicalouterBe,"physicalouterBe",logicalPip,false,0);

cout << "Ok 6"<< endl;
  //the volume of inner side Be layer
   G4Tubs* solidinnerBeSide = new G4Tubs("solidinnerBeSide",innerBeSide[0],innerBeSide[1],innerBeSide[2]/2,0,360);
  logicalinnerBeSide = new G4LogicalVolume(solidinnerBeSide, G4Material::GetMaterial("Beryllium"),"logicalinnerBeSide");
  physicalinnerBeSide = new G4PVPlacement(0,G4ThreeVector(0,0,innerBeSide[3]),logicalinnerBeSide,"physicalinnerBeSide1",logicalPip,false,0);
                        new G4PVPlacement(0,G4ThreeVector(0,0,innerBeSide[4]),logicalinnerBeSide,"physicalinnerBeSide2",logicalPip,false,1);

  //the volume of inner Al layer
   G4Tubs* solidinnerAl = new G4Tubs("solidinnerAl",innerAl[0],innerAl[1],innerAl[2]/2,0,360);
  logicalinnerAl = new G4LogicalVolume(solidinnerAl, G4Material::GetMaterial("Aluminium"),"logicalinnerAl");
  physicalinnerAl = new G4PVPlacement(0,G4ThreeVector(0,0,innerAl[3]),logicalinnerAl,"physicalinnerAl1",logicalPip,false,0);
                        new G4PVPlacement(0,G4ThreeVector(0,0,innerAl[4]),logicalinnerAl,"physicalinnerAl2",logicalPip,false,1);

  //the volume of Al ring
   G4Tubs* solidAlRing = new G4Tubs("solidAlRing",AlRing[0],AlRing[1],AlRing[2]/2,0,360);
  logicalAlRing = new G4LogicalVolume(solidAlRing, G4Material::GetMaterial("Aluminium"),"logicalAlRing");
  physicalAlRing = new G4PVPlacement(0,G4ThreeVector(0,0,AlRing[3]),logicalAlRing,"physicalAlRing1",logicalPip,false,0);
                        new G4PVPlacement(0,G4ThreeVector(0,0,AlRing[4]),logicalAlRing,"physicalAlRing2",logicalPip,false,1);
cout << "Ok 6"<< endl;

  //the volume of outer Al layer
   G4Tubs* solidouterAl = new G4Tubs("solidouterAl",outerAl[0],outerAl[1],outerAl[2]/2,0,360);
  logicalouterAl = new G4LogicalVolume(solidouterAl, G4Material::GetMaterial("Aluminium"),"logicalouterAl");
  physicalouterAl = new G4PVPlacement(0,G4ThreeVector(0,0,outerAl[3]),logicalouterAl,"physicalouterAl1",logicalPip,false,0);
                        new G4PVPlacement(0,G4ThreeVector(0,0,outerAl[4]),logicalouterAl,"physicalouterAl2",logicalPip,false,1);

  //the volume of Ag layer
   G4Tubs* solidAgLayer = new G4Tubs("solidAgLayer",AgLayer[0],AgLayer[1],AgLayer[2]/2,0,360);
  logicalAgLayer = new G4LogicalVolume(solidAgLayer, Ag,"logicalAgLayer");
  physicalAgLayer = new G4PVPlacement(0,G4ThreeVector(0,0,AgLayer[3]),logicalAgLayer,"physicalAgLayer1",logicalPip,false,0);
                        new G4PVPlacement(0,G4ThreeVector(0,0,AgLayer[4]),logicalAgLayer,"physicalAgLayer2",logicalPip,false,1);
cout << "Ok 6"<< endl;

  //the volume of  Cu layer
   G4Tubs* solidCuLayer = new G4Tubs("solidCuLayer",CuLayer[0],CuLayer[1],CuLayer[2]/2,0,360);
  logicalCuLayer = new G4LogicalVolume(solidCuLayer, G4Material::GetMaterial("Copper"),"logicalCuLayer");
  physicalCuLayer = new G4PVPlacement(0,G4ThreeVector(0,0,CuLayer[3]),logicalCuLayer,"physicalCuLayer1",logicalPip,false,0);
                        new G4PVPlacement(0,G4ThreeVector(0,0,CuLayer[4]),logicalCuLayer,"physicalCuLayer2",logicalPip,false,1);

  //the volume of Al cover
   G4Tubs* solidAlCover = new G4Tubs("solidAlCover",AlCover[0],AlCover[1],AlCover[2]/2,0,360);
  logicalAlCover = new G4LogicalVolume(solidAlCover, G4Material::GetMaterial("Aluminium"),"logicalAlCover");
  physicalAlCover = new G4PVPlacement(0,G4ThreeVector(0,0,AlCover[3]),logicalAlCover,"physicalAlCover1",logicalPip,false,0);
                        new G4PVPlacement(0,G4ThreeVector(0,0,AlCover[4]),logicalAlCover,"physicalAlCover2",logicalPip,false,1);

cout << "Ok 6"<< endl;
  G4VisAttributes* visPip = new G4VisAttributes(G4Colour(0.,0.,1.));
  logicalPip->SetVisAttributes(visPip);
  //logicalPip->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* visgold = new G4VisAttributes(G4Colour(1.,1.,0.));
  logicalgoldLayer->SetVisAttributes(visgold);
//  logicalgoldLayer->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* visBe = new G4VisAttributes(G4Colour(0.,1.,0.));
  logicalinnerBe->SetVisAttributes(visBe);
//  logicalinnerBe->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* visOil = new G4VisAttributes(G4Colour(1.,1.,1.));
  logicaloilLayer->SetVisAttributes(visOil);
//  logicaloilLayer->SetVisAttributes(G4VisAttributes::Invisible);
  logicalouterBe->SetVisAttributes(visBe);
//  logicalouterBe->SetVisAttributes(G4VisAttributes::Invisible);
  logicalinnerBeSide->SetVisAttributes(visBe);
//  logicalinnerBeSide->SetVisAttributes(G4VisAttributes::Invisible);
  logicalinnerAl->SetVisAttributes(visOil);
//  logicalinnerAl->SetVisAttributes(G4VisAttributes::Invisible);
  logicalAlRing->SetVisAttributes(visOil);
//  logicalAlRing->SetVisAttributes(G4VisAttributes::Invisible);
  logicalouterAl->SetVisAttributes(visOil);
//  logicalouterAl->SetVisAttributes(G4VisAttributes::Invisible);
  logicalAgLayer->SetVisAttributes(visBe);
//  logicalAgLayer->SetVisAttributes(G4VisAttributes::Invisible);
  logicalCuLayer->SetVisAttributes(visPip);
//  logicalCuLayer->SetVisAttributes(G4VisAttributes::Invisible);
  logicalAlCover->SetVisAttributes(visOil);
//  logicalAlCover->SetVisAttributes(G4VisAttributes::Invisible);
cout << "Ok 7"<< endl;
}
