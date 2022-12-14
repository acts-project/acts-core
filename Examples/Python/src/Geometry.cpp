// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <memory>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace Acts::Python {
void addGeometry(Context& ctx) {
  auto m = ctx.get("main");
  {
    py::class_<Acts::Surface, std::shared_ptr<Acts::Surface>>(m, "Surfaces")
        .def("volumeId",
             [](Acts::Surface& self) { return self.geometryId().volume(); })
        .def("boundaryId",
             [](Acts::Surface& self) { return self.geometryId().boundary(); })
        .def("layerId",
             [](Acts::Surface& self) { return self.geometryId().layer(); })
        .def("approachId",
             [](Acts::Surface& self) { return self.geometryId().approach(); })
        .def("sensitiveId",
             [](Acts::Surface& self) { return self.geometryId().sensitive(); })
        .def("extraId",
             [](Acts::Surface& self) { return self.geometryId().extra(); });
  }
  {
    py::class_<Acts::TrackingGeometry, std::shared_ptr<Acts::TrackingGeometry>>(
        m, "TrackingGeometry")
        .def("visitSurfaces",
             [](Acts::TrackingGeometry& self, py::function& func) {
               self.visitSurfaces(func);
             });
  }
  {
    py::class_<Acts::GeometryIdentifierHook, std::shared_ptr<Acts::GeometryIdentifierHook>>(
        m, "GeometryIdentifierHook");
  }
}

}  // namespace Acts::Python
