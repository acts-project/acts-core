// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"

#include <type_traits>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;

namespace Acts::Python {

void addEventData(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  py::class_<Acts::ParticleHypothesis>(m, "ParticleHypothesis")
      .def(py::init([](Acts::PdgParticle absPdg, float mass, float absCharge) {
             return Acts::ParticleHypothesis(absPdg, mass,
                                             AnyCharge{absCharge});
           }),
           py::arg("pdg"), py::arg("mass"), py::arg("absCharge"))
      .def(py::init([](std::underlying_type_t<Acts::PdgParticle> absPdg,
                       float mass, float absCharge) {
             return Acts::ParticleHypothesis(
                 static_cast<Acts::PdgParticle>(absPdg), mass,
                 AnyCharge{absCharge});
           }),
           py::arg("absPdg"), py::arg("mass"), py::arg("absCharge"))
      .def("__str__",
           [](const Acts::ParticleHypothesis& particleHypothesis) {
             std::stringstream os;
             particleHypothesis.toStream(os);
             return os.str();
           })
      .def("absolutePdg",
           [](const Acts::ParticleHypothesis& p) { return p.absolutePdg(); })
      .def("mass", [](const Acts::ParticleHypothesis& p) { return p.mass(); })
      .def("absoluteCharge",
           [](const Acts::ParticleHypothesis& p) { return p.absoluteCharge(); })
      .def_property_readonly_static("muon",
                                    [](py::object /* self */) {
                                      return Acts::ParticleHypothesis::muon();
                                    })
      .def_property_readonly_static("pion",
                                    [](py::object /* self */) {
                                      return Acts::ParticleHypothesis::pion();
                                    })
      .def_property_readonly_static(
          "electron",
          [](py::object /* self */) {
            return Acts::ParticleHypothesis::electron();
          })
      .def_property_readonly_static("kaon",
                                    [](py::object /* self */) {
                                      return Acts::ParticleHypothesis::kaon();
                                    })
      .def_property_readonly_static("proton",
                                    [](py::object /* self */) {
                                      return Acts::ParticleHypothesis::proton();
                                    })
      .def_property_readonly_static(
          "geantino",
          [](py::object /* self */) {
            return Acts::ParticleHypothesis::geantino();
          })
      .def_property_readonly_static(
          "chargedGeantino", [](py::object /* self */) {
            return Acts::ParticleHypothesis::chargedGeantino();
          });
}

}  // namespace Acts::Python
