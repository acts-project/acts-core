// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Fatras/FatrasSimulation.hpp"

#include <stdexcept>

#include "boost/program_options.hpp"

void ActsExamples::FatrasSimulation::addOptions(Options::Description& desc) {
  using boost::program_options::bool_switch;
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("fatras-pmin-gev", value<double>()->default_value(0.5),
      "Minimum momentum for simulated particles in GeV");
  opt("fatras-em-scattering", value<bool>()->default_value(true),
      "Simulate multiple scattering of charged particles");
  opt("fatras-em-ionisation", value<bool>()->default_value(true),
      "Simulate ionisiation/excitation energy loss of charged particles");
  opt("fatras-em-radiation", value<bool>()->default_value(true),
      "Simulate radiative energy loss of charged particles");
  opt("fatras-em-photonconversion", value<bool>()->default_value(true),
      "Simulate electron-positron pair production by photon conversion");
  opt("fatras-hits",
      value<std::string>()
          ->value_name("none|sensitive|material|all")
          ->default_value("sensitive"),
      "Which surfaces should record charged particle hits");
}

ActsExamples::FatrasSimulation::Config
ActsExamples::FatrasSimulation::readConfig(const Options::Variables& vars) {
  using namespace Acts::UnitLiterals;

  ActsExamples::FatrasSimulation::Config cfg;

  cfg.pMin = vars["fatras-pmin-gev"].as<double>() * 1_GeV;

  cfg.emScattering = vars["fatras-em-scattering"].as<bool>();
  cfg.emEnergyLossIonisation = vars["fatras-em-ionisation"].as<bool>();
  cfg.emEnergyLossRadiation = vars["fatras-em-radiation"].as<bool>();
  cfg.emPhotonConversion = vars["fatras-em-photonconversion"].as<bool>();

  // select hit surfaces for charged particles
  const std::string hits = vars["fatras-hits"].as<std::string>();
  if (hits == "sensitive") {
    cfg.generateHitsOnSensitive = true;
    cfg.generateHitsOnMaterial = false;
    cfg.generateHitsOnPassive = false;
  } else if (hits == "material") {
    cfg.generateHitsOnSensitive = false;
    cfg.generateHitsOnMaterial = true;
    cfg.generateHitsOnPassive = false;
  } else if (hits == "all") {
    cfg.generateHitsOnSensitive = true;
    cfg.generateHitsOnMaterial = true;
    cfg.generateHitsOnPassive = true;
  } else {
    throw std::runtime_error("Invalid Fatras hits selection '" + hits + "'");
  }

  return cfg;
}