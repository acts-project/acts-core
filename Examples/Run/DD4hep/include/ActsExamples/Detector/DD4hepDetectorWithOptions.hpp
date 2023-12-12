// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

#include <memory>
#include <utility>

namespace Acts {
class IMaterialDecorator;
}  // namespace Acts
namespace boost {
namespace program_options {
class options_description;
}  // namespace program_options
}  // namespace boost

namespace ActsExamples {

class DD4hepDetectorWithOptions : public IBaseDetector {
  DD4hep::DD4hepDetector m_detector;

 public:
  void addOptions(
      boost::program_options::options_description& opt) const override;

  auto finalize(const boost::program_options::variables_map& vm,
                std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
      -> std::pair<TrackingGeometryPtr, ContextDecorators> override;
};

}  // namespace ActsExamples
