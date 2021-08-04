// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepDetectorOptions.hpp"
#include "ActsExamples/Geant4/G4DetectorConstructionFactory.hpp"
#include "ActsExamples/Geant4DD4hep/DD4hepDetectorConstruction.hpp"

#include <functional>
#include <memory>

#include <G4VUserDetectorConstruction.hh>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

PYBIND11_MODULE(ActsPythonBindingsGeant4DD4hep, m) {
  py::module_::import("acts.ActsPythonBindingsGeant4");

  py::class_<DD4hepDetectorConstructionFactory, G4DetectorConstructionFactory,
             std::shared_ptr<DD4hepDetectorConstructionFactory>>(
      m, "DD4hepDetectorConstructionFactory")
      .def(py::init([](DD4hep::DD4hepGeometryService& geometrySvc) {
        return std::make_shared<DD4hepDetectorConstructionFactory>(
            *geometrySvc.lcdd());
      }));
}
