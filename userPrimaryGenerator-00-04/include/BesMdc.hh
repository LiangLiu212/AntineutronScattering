/*****************************************************************************\
 *  BOSS (BES Offline Software System)                                       *
 *       Filename:  BesMdc.hh
 *                                                                           *
 *    Description:                                                           
 *                                                                           *
 *        Version:  1.0                                                      *
 *        Created:  11/08/2022 03:03:47 PM                                   *
 *       Revision:  none                                                     *
 *       Compiler:  gcc                                                      *
 *                                                                           *
 *         Author:  Liang Liu (USTC), <liangzy@mail.ustc.edu.cn>             *
 *   Organization:  BESIII Collaboration                                     *
 *        License:  GNU Lesser General Public License 3.0, see LICENSE.md.   *
\*****************************************************************************/


#ifndef BES_MDC_H
#define BES_MDC_H

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"
class G4VPhysicalVolume;
#include "BesSubdetector.hh"
class G4LogicalVolume;

class BesMdc : public BesSubdetector {

		public:
				BesMdc();
				~BesMdc();

				void Construct(G4LogicalVolume* logicbes);
				void DefineMaterial();

		private:

				G4Material *Al, *CarbonFiber, *Air;

				G4LogicalVolume *logicalMdc;
				G4VPhysicalVolume *physicalMdc;

				G4LogicalVolume *logicalinnerAl;
				G4VPhysicalVolume *physicalinnerAl;

				G4LogicalVolume *logicalCarbonFiber;
				G4VPhysicalVolume *physicalCarbonFiber;

				G4LogicalVolume *logicalouterAl;
				G4VPhysicalVolume *physicalouterAl;


};

#endif

