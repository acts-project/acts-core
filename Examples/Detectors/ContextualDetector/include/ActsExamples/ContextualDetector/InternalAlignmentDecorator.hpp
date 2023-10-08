// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/ContextualDetector/AlignmentDecorator.hpp"
#include "ActsExamples/ContextualDetector/InternallyAlignedDetectorElement.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include <mutex>
#include <vector>

namespace ActsExamples {

namespace Contextual {

/// @brief A mockup service that rotates the modules in a
/// simple tracking geometry
///
/// It acts on the PayloadDetectorElement, i.e. the
/// geometry context carries the full transform store (payload)
class InternalAlignmentDecorator : public AlignmentDecorator {
 public:
  using LayerStore =
      std::vector<std::shared_ptr<InternallyAlignedDetectorElement>>;
  using DetectorStore = std::vector<LayerStore>;

  /// @brief nested configuration struct
  struct Config : public AlignmentDecorator::Config {
    /// The detector store (filled at creation creation)
    DetectorStore detectorStore;
  };

// Adding definition of applyMisalignment
void applyMisalignment(std::shared_ptr<InternallyAlignedDetectorElement>& element,
                       const AlignmentMisalignment& misalignment);


  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param logger The logging framework
  InternalAlignmentDecorator(const Config& cfg,
                             std::unique_ptr<const Acts::Logger> logger =
                                 Acts::getDefaultLogger("AlignmentDecorator",
                                                        Acts::Logging::INFO));

  /// Virtual destructor
  virtual ~InternalAlignmentDecorator() = default;

  /// @brief decorates (adds, modifies) the AlgorithmContext
  /// with a geometric rotation per event
  ///
  /// @note If decorators depend on each other, they have to be
  /// added in order.
  ///
  /// @param context the bare (or at least non-const) Event context
  ProcessCode decorate(AlgorithmContext& context) final override;

  /// @brief decorator name() for screen output
  const std::string& name() const final override { return m_name; }

 private:
  Config m_cfg;                                  ///< the configuration class
  std::unique_ptr<const Acts::Logger> m_logger;  ///!< the logging instance
  std::string m_name = "AlignmentDecorator";

  ///< Protect multiple alignments to be loaded at once
  std::mutex m_alignmentMutex;
  struct IovStatus {
    size_t lastAccessed;
  };
  std::unordered_map<unsigned int, IovStatus> m_activeIovs;
  size_t m_eventsSeen{0};

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }


   // Struct that takes care of misalignment (random translation + rotation)
  struct AlignmentMisalignment {
    Acts::Translation3 translation; 
    Acts::Rotation3 rotation;       
  };

// Adding function that generates misalignment
  AlignmentMisalignment generateMisalignment();

  // Adding function that applies  misalignment and fills transform deltas
  void applyMisalignment(std::shared_ptr<InternallyAlignedDetectorElement>& element,
                         const AlignmentMisalignment& misalignment);


  AlignmentMisalignment m_alignmentMisalignment;
};

}  // namespace Contextual
}  // namespace ActsExamples
