// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvPropagationStepsWriter.hpp"

namespace ActsExamples {

ActsExamples::ProcessCode CsvPropagationStepsWriter::writeT(
    const ActsExamples::AlgorithmContext& context,
    const StepsVector& stepCollection) {
  // open per-event file
  std::string path = ActsExamples::perEventFilepath(
      m_cfg.outputDir, "propagation-steps.csv", context.eventNumber);

  std::ofstream os(path, std::ofstream::out | std::ofstream::trunc);
  if (!os) {
    throw std::ios_base::failure("Could not open '" + path + "' to write");
  }

  os << "ssaccuracy,ssactor,ssaborter,ssuser,x,y,z,px,py,pz,surface,volume\n";
  os << std::setprecision(m_cfg.outputPrecision);

  for (const auto& steps : stepCollection) {
    for (const auto& step : steps) {
      os << step.stepSize.value(Acts::ConstrainedStep::Type::accuracy) << ","
         << step.stepSize.value(Acts::ConstrainedStep::Type::actor) << ","
         << step.stepSize.value(Acts::ConstrainedStep::Type::aborter) << ","
         << step.stepSize.value(Acts::ConstrainedStep::Type::user) << ","
         << step.position.x() << "," << step.position.y() << ","
         << step.position.z() << "," << step.momentum.x() << ","
         << step.momentum.y() << "," << step.momentum.z() << ",";

      if (step.surface)
        os << step.surface->geometryId().value() << ",";
      else
        os << "na,";

      if (step.volume)
        os << step.volume->volumeName() << "\n";
      else
        os << "na\n";
    }
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
