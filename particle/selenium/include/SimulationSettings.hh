// General Simulation Settings

#ifndef SimulationSettings_H
#define SimulationSettings_H

#include "G4SystemOfUnits.hh"

#define OUTPUTDIR "C:\\Users\\alexp\\Documents\\UW\\Research\\Selenium\\aSe0vBB\\particle\\selenium_detector\\ouput\\"


// Defining geometries of the problem
const G4double seleniumXY = 10. * mm;
const G4double seleniumDepth = 0.5 * mm;
const G4double glassDepth = 1. * mm;
const G4double copperZ = 12.5 * mm;
const G4double copperX1 = 17.8 * mm;
const G4double copperX2 = 9.5 * mm;
const G4double copperY = 19.1 * mm;
const G4double copperTraceWidth = 0.3 * mm;
const G4double airSpacingZ = 5.4 * mm;
const G4double electrodeDepth = 1. * um;
const G4double pixelRadius = 1. * mm;
const G4double guardInnerRadius = pixelRadius + 0.088 * mm;
const G4double guardOuterRadius = guardInnerRadius + 0.1 * mm;
const G4double collimatorRadius = 0.5 * mm;

// Sensitive Detector parameters
const G4String seleniumHitsCollectionName = "SeleniumHitsCollection";

// Defining default units to use to keep analysis consistent
const double DEFAULT_LENGTH = mm;
const double DEFAULT_ENERGY = keV;

const bool SIMPLE = false;

#endif
