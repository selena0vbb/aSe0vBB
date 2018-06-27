// General Simulation Settings

#ifndef SimulationSettings_H
#define SimulationSettings_H

#include "G4SystemofUnits.hh"

#define OUTPUTDIR "C:\\Users\\alexp\\Documents\\UW\\Research\\Selenium\\aSe0vBB\\particle\\selenium_detector\\ouput\\"


// Defining geometries of the problem
const G4double seleniumXY = 10. * mm;
const G4double seleniumDepth = 0.2 * mm;
const G4double goldElectrodeArea = 2 * mm2;
const G4double goldElectrodeDepth = 10 * nm;

// Sensitive Detector parameters
const G4String seleniumHitsCollectionName = "SeleniumHitsCollection";

#endif