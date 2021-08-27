// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Options/Pythia8Options.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/PdgParticle.hpp"
#include "ActsExamples/Generators/MultiplicityGenerators.hpp"
#include "ActsExamples/Generators/Pythia8ProcessGenerator.hpp"
#include "ActsExamples/Generators/VertexGenerators.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <vector>

#include <boost/program_options.hpp>

void ActsExamples::Options::addPythia8Options(Description& desc) {
  using boost::program_options::bool_switch;
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("gen-nhard", value<size_t>()->default_value(1u),
      "Number of hard interactions, zero to disable");
  opt("gen-npileup", value<double>()->default_value(200.0),
      "Mean number of pile-up interactions, zero to disable");
  opt("gen-vertex-xy-std-mm", value<double>()->default_value(0.0125),
      "Transverse vertex standard deviation in mm");
  opt("gen-vertex-z-std-mm", value<double>()->default_value(55.5),
      "Longitudinal vertex standard deviation in mm");
  opt("gen-vertex-t-std-ns", value<double>()->default_value(5.0),
      "Temporal vertex standard deviation in ns");
  opt("gen-cms-energy-gev", value<double>()->default_value(14000.0),
      "Center-of-mass energy collision in GeV");
  opt("gen-pdg-beam0",
      value<int32_t>()->default_value(Acts::PdgParticle::eProton),
      "PDG number of the first beam particle");
  opt("gen-pdg-beam1",
      value<int32_t>()->default_value(Acts::PdgParticle::eProton),
      "PDG number of the second beam particle");
  opt("gen-hard-process",
      value<std::vector<std::string>>()->default_value({"HardQCD:all = on"},
                                                       "HardQCD:all = on"),
      "Pythia8 process string for the hard interactions. Can be given multiple "
      "times.");
  opt("gen-pileup-process",
      value<std::vector<std::string>>()->default_value({"SoftQCD:all = on"},
                                                       "SoftQCD:all = on"),
      "Pythi8 process string for the pile-up interactions. Can be given "
      "multiple times.");
}

ActsExamples::EventGenerator::Config ActsExamples::Options::readPythia8Options(
    const Variables& vars, Acts::Logging::Level lvl) {
  using namespace Acts::UnitLiterals;

  const auto nhard = vars["gen-nhard"].as<size_t>();
  const auto npileup = vars["gen-npileup"].as<double>();
  const auto pdgBeam0 = static_cast<Acts::PdgParticle>(
      vars["gen-pdg-beam0"].template as<int32_t>());
  const auto pdgBeam1 = static_cast<Acts::PdgParticle>(
      vars["gen-pdg-beam1"].template as<int32_t>());
  const auto cmsEnergy = vars["gen-cms-energy-gev"].as<double>() * 1_GeV;

  auto vertexGen = std::make_shared<GaussianVertexGenerator>();
  vertexGen->stddev[Acts::ePos0] =
      vars["gen-vertex-xy-std-mm"].as<double>() * 1_mm;
  vertexGen->stddev[Acts::ePos1] =
      vars["gen-vertex-xy-std-mm"].as<double>() * 1_mm;
  vertexGen->stddev[Acts::ePos2] =
      vars["gen-vertex-z-std-mm"].as<double>() * 1_mm;
  vertexGen->stddev[Acts::eTime] =
      vars["gen-vertex-t-std-ns"].as<double>() * 1_ns;

  EventGenerator::Config cfg;
  if (0u < nhard) {
    Pythia8Generator::Config hard;
    hard.pdgBeam0 = pdgBeam0;
    hard.pdgBeam1 = pdgBeam1;
    hard.cmsEnergy = cmsEnergy;
    hard.settings = vars["gen-hard-process"].as<std::vector<std::string>>();

    cfg.generators.push_back(
        {std::make_shared<FixedMultiplicityGenerator>(nhard), vertexGen,
         std::make_shared<Pythia8Generator>(hard, lvl)});
  }
  if (0.0 < npileup) {
    Pythia8Generator::Config pileup;
    pileup.pdgBeam0 = pdgBeam0;
    pileup.pdgBeam1 = pdgBeam1;
    pileup.cmsEnergy = cmsEnergy;
    pileup.settings = vars["gen-pileup-process"].as<std::vector<std::string>>();

    cfg.generators.push_back(
        {std::make_shared<PoissonMultiplicityGenerator>(npileup), vertexGen,
         std::make_shared<Pythia8Generator>(pileup, lvl)});
  }

  return cfg;
}
