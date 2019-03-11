// General Simulation Settings

#ifndef SimulationSettings_H
#define SimulationSettings_H

#include "G4SystemOfUnits.hh"

#define OUTPUTDIR "C:\\Users\\alexp\\Documents\\UW\\Research\\Selenium\\aSe0vBB\\particle\\selenium_detector\\ouput\\"


// Defining geometries of the problem
const G4double seleniumXY = 10. * mm;
const G4double seleniumDepth = 0.2 * mm;
const G4double glassDepth = 1. * mm;
const G4double copperDepth = 20. * mm;
const G4double copperXY = 40. * mm;
const G4double airSpacingZ = 30. * mm;
const G4double electrodeDepth = 1. * um;
const G4double pixelRadius = 1. * mm;
const G4double guardInnerRadius = pixelRadius + 0.088 * mm;
const G4double guardOuterRadius = guardInnerRadius + 0.1 * mm;

// Sensitive Detector parameters
const G4String seleniumHitsCollectionName = "SeleniumHitsCollection";

// Defining default units to use to keep analysis consistent
const double DEFAULT_LENGTH = mm;
const double DEFAULT_ENERGY = keV;

const bool SIMPLE = false;

#endif
