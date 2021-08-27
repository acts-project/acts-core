// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Options/ParticleGunOptions.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Generators/MultiplicityGenerators.hpp"
#include "ActsExamples/Generators/ParametricParticleGenerator.hpp"
#include "ActsExamples/Generators/VertexGenerators.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <boost/program_options.hpp>

void ActsExamples::Options::addParticleGunOptions(Description& desc) {
  using boost::program_options::bool_switch;
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("gen-vertex-xy-std-mm", value<double>()->default_value(0.0),
      "Transverse vertex standard deviation in mm");
  opt("gen-vertex-z-std-mm", value<double>()->default_value(0.0),
      "Longitudinal vertex standard deviation in mm");
  opt("gen-vertex-t-std-ns", value<double>()->default_value(0.0),
      "Temporal vertex standard deviation in ns");
  opt("gen-phi-degree",
      value<Interval>()->value_name("MIN:MAX")->default_value({0.0, 360.0}),
      "Transverse direction angle generation range in degree");
  opt("gen-eta",
      value<Interval>()->value_name("MIN:MAX")->default_value({-4.0, 4.0}),
      "Pseudo-rapidity generation range");
  opt("gen-mom-gev",
      value<Interval>()->value_name("MIN:MAX")->default_value({1.0, 10.0}),
      "Absolute (or transverse) momentum generation range in GeV");
  opt("gen-mom-transverse", bool_switch(),
      "Momentum referse to transverse momentum");
  opt("gen-pdg", value<int32_t>()->default_value(Acts::PdgParticle::eMuon),
      "PDG number of the particle, will be adjusted for charge flip.");
  opt("gen-randomize-charge", bool_switch(),
      "Flip the charge and change the PDG number accordingly.");
  opt("gen-nparticles", value<size_t>()->default_value(1u),
      "Number of generated particles per vertex");
  opt("gen-nvertices", value<size_t>()->default_value(1u),
      "Number of generated vertices");
}

ActsExamples::EventGenerator::Config
ActsExamples::Options::readParticleGunOptions(const Variables& vars) {
  using namespace Acts::UnitLiterals;

  // access user config w/ unit conversion
  auto getValue = [&](const char* name, auto unit) {
    return vars[name].as<double>() * unit;
  };
  auto getRange = [&](const char* name, auto unit, auto& lower, auto& upper) {
    auto interval = vars[name].as<Options::Interval>();
    lower = interval.lower.value() * unit;
    upper = interval.upper.value() * unit;
  };

  auto vertexGen = std::make_shared<GaussianVertexGenerator>();
  vertexGen->stddev[Acts::ePos0] = getValue("gen-vertex-xy-std-mm", 1_mm);
  vertexGen->stddev[Acts::ePos1] = getValue("gen-vertex-xy-std-mm", 1_mm);
  vertexGen->stddev[Acts::ePos2] = getValue("gen-vertex-z-std-mm", 1_mm);
  vertexGen->stddev[Acts::eTime] = getValue("gen-vertex-t-std-ns", 1_ns);

  ParametricParticleGenerator::Config pgCfg;
  getRange("gen-phi-degree", 1_degree, pgCfg.phiMin, pgCfg.phiMax);
  // user config sets eta but the generator takes theta
  double etaMin, etaMax;
  getRange("gen-eta", 1.0, etaMin, etaMax);
  pgCfg.thetaMin = 2 * std::atan(std::exp(-etaMin));
  pgCfg.thetaMax = 2 * std::atan(std::exp(-etaMax));
  getRange("gen-mom-gev", 1_GeV, pgCfg.pMin, pgCfg.pMax);
  pgCfg.pTransverse = vars["gen-mom-transverse"].template as<bool>();
  pgCfg.pdg =
      static_cast<Acts::PdgParticle>(vars["gen-pdg"].template as<int32_t>());
  pgCfg.randomizeCharge = vars["gen-randomize-charge"].template as<bool>();
  pgCfg.numParticles = vars["gen-nparticles"].as<size_t>();

  if (pgCfg.numParticles > std::pow(2, ActsFatras::Barcode::bits(2))) {
    throw std::runtime_error{
        "Too many particles per vertex requested for Fatras Barcode"};
  }

  size_t nVertices = vars["gen-nvertices"].as<size_t>();

  if (nVertices > std::pow(2, ActsFatras::Barcode::bits(0))) {
    throw std::runtime_error{"Too many vertices requested for Fatras Barcode"};
  }

  auto mGen = std::make_shared<FixedMultiplicityGenerator>();
  mGen->n = nVertices;

  EventGenerator::Config cfg;
  cfg.generators = {
      {std::move(mGen), std::move(vertexGen),
       std::make_shared<ParametricParticleGenerator>(pgCfg)},
  };

  return cfg;
}
