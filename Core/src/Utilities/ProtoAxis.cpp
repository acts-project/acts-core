// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/ProtoAxis.hpp"

namespace {
void checkConsistency(Acts::AxisDirection aDir, Acts::AxisBoundaryType abType) {
  if (abType == Acts::AxisBoundaryType::Closed &&
      aDir != Acts::AxisDirection::AxisPhi &&
      aDir != Acts::AxisDirection::AxisRPhi) {
    std::string msg =
        "ProtoBinning: Invalid axis boundary type 'Closed' for direction '";
    msg += axisDirectionName(aDir) + "'.";
    throw std::invalid_argument(msg);
  }
}
}  // namespace

Acts::ProtoAxis::ProtoAxis(AxisDirection aDir, Acts::AxisBoundaryType abType,
                           const std::vector<double>& edges)
    : m_axisDir(aDir), m_axis(IAxis::create(abType, edges)) {
  checkConsistency(aDir, abType);
}

Acts::ProtoAxis::ProtoAxis(AxisDirection aDir, AxisBoundaryType abType,
                           double minE, double maxE, std::size_t nbins)
    : m_axisDir(aDir), m_axis(IAxis::create(abType, minE, maxE, nbins)) {
  checkConsistency(aDir, abType);
}

Acts::ProtoAxis::ProtoAxis(AxisDirection aDir, AxisBoundaryType abType,
                           std::size_t nbins)
    : m_axisDir(aDir),
      m_axis(IAxis::create(abType, 0., 1., nbins)),
      m_autorange(true) {
  checkConsistency(aDir, abType);
}

Acts::AxisDirection Acts::ProtoAxis::getAxisDirection() const {
  return m_axisDir;
}

const Acts::IAxis& Acts::ProtoAxis::getAxis() const {
  return *m_axis;
}

bool Acts::ProtoAxis::isAutorange() const {
  return m_autorange;
}

void Acts::ProtoAxis::toStream(std::ostream& os) const {
  os << toString();
}

std::string Acts::ProtoAxis::toString() const {
  std::stringstream ss;
  ss << "ProtoAxis: " << getAxis().getNBins() << " bins in "
     << axisDirectionName(m_axisDir);
  ss << (getAxis().getType() == AxisType::Variable ? ", variable "
                                                   : ", equidistant ");
  if (!m_autorange) {
    const auto& edges = getAxis().getBinEdges();
    ss << "within [" << edges.front() << ", " << edges.back() << "]";
  } else {
    ss << "within automatic range";
  }
  return ss.str();
}
