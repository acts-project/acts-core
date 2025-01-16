// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/FpeMonitoring/FpeMonitor.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/IWriter.hpp"
#include "ActsExamples/Framework/SequenceElement.hpp"
#include "ActsExamples/Utilities/tbbWrap.hpp"
#include <Acts/Utilities/Logger.hpp>

#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <tbb/enumerable_thread_specific.h>

namespace ActsExamples {
class DataHandleBase;
class IAlgorithm;
class IContextDecorator;
class IReader;
class IWriter;
class SequenceElement;

using IterationCallback = void (*)();

/// Custom exception class so FPE failures can be caught
class FpeFailure : public std::runtime_error {
  using std::runtime_error::runtime_error;
};

class SequenceConfigurationException : public std::runtime_error {
 public:
  SequenceConfigurationException()
      : std::runtime_error{"Sequence configuration error"} {}
};

/// A simple algorithm sequencer for event processing.
///
/// This is the backbone of the framework. It reads events from file,
/// runs the configured algorithms for each event, and writes selected data
/// back to a file.
class Sequencer {
 public:
  struct FpeMask {
    std::string file;
    std::pair<std::size_t, std::size_t> lines;
    Acts::FpeType type;
    std::size_t count;
  };

  struct Config {
    /// number of events to skip at the beginning
    std::size_t skip = 0;
    /// number of events to process, std::numeric_limits<std::size_t>::max() to
    /// process all available events
    std::optional<std::size_t> events = std::nullopt;
    /// logging level
    Acts::Logging::Level logLevel = Acts::Logging::INFO;
    /// number of parallel threads to run, negative for automatic
    /// determination
    int numThreads = -1;
    /// output directory for timing information, empty for working directory
    std::string outputDir;
    /// output name of the timing file
    std::string outputTimingFile = "timing.csv";
    /// Callback that is invoked in the event loop.
    /// @warning This function can be called from multiple threads and should therefore be thread-safe
    IterationCallback iterationCallback = []() {};

    bool trackFpes = true;
    std::vector<FpeMask> fpeMasks{};
    bool failOnFirstFpe = false;
    std::size_t fpeStackTraceLength = 8;
  };

  explicit Sequencer(const Config &cfg);

  /// Add a context decorator to the set of context decorators.
  ///
  /// @throws std::invalid_argument if the decorator is NULL.
  void addContextDecorator(std::shared_ptr<IContextDecorator> decorator);

  /// Add a reader to the set of readers.
  ///
  /// @throws std::invalid_argument if the reader is NULL.
  void addReader(std::shared_ptr<IReader> reader);

  /// Append an algorithm to the sequence of algorithms.
  ///
  /// @throws std::invalid_argument if the algorithm is NULL.
  void addAlgorithm(std::shared_ptr<IAlgorithm> algorithm);

  /// Append a sequence element to the sequence
  ///
  /// @throws std::invalid_argument if the element is NULL.
  void addElement(const std::shared_ptr<SequenceElement> &element);

  /// Add a writer to the set of writers.
  ///
  /// @throws std::invalid_argument if the writer is NULL.
  void addWriter(std::shared_ptr<IWriter> writer);

  /// Add an alias to the whiteboard.
  void addWhiteboardAlias(const std::string &aliasName,
                          const std::string &objectName);

  Acts::FpeMonitor::Result fpeResult() const;

  /// Run the event loop.
  ///
  /// @return status code compatible with the `main()` return code
  /// @returns EXIT_SUCCESS when everying worked without problems
  /// @returns EXIT_FAILURE if something went wrong
  ///
  /// @note If the number of events to process is undefined, the sequencer
  /// will process events until the first reader signals the end-of-file. If
  /// given, it sets an upper bound.
  ///
  /// This function is intended to be run as the last thing in the tool
  /// main function and its return value can be used directly as the program
  /// return value, i.e.
  ///
  ///     int main(int argc, char* argv[])
  ///     {
  ///         Sequencer::Config cfg;
  ///         ... // configure the sequencer
  ///         Sequencer seq;
  ///         ... // set up the algorithms
  ///         return seq.run();
  ///     }
  ///
  /// This will run the start-of-run hook for all configured services, run all
  /// configured readers, algorithms, and writers for each event, then invoke
  /// the end-of-run hook for all configured writers.
  int run();

  /// Get const access to the config
  const Config &config() const { return m_cfg; }

 private:
  /// List of all configured algorithm names.
  std::vector<std::string> listAlgorithmNames() const;
  /// Determine range of (requested) events;
  /// [std::numeric_limits<std::size_t>::max(),
  /// std::numeric_limits<std::size_t>::max()) for error.
  std::pair<std::size_t, std::size_t> determineEventsRange() const;

  std::pair<std::string, std::size_t> fpeMaskCount(
      const boost::stacktrace::stacktrace &st, Acts::FpeType type) const;

  void fpeReport() const;

  struct SequenceElementWithFpeResult {
    std::shared_ptr<SequenceElement> sequenceElement;
    tbb::enumerable_thread_specific<Acts::FpeMonitor::Result> fpeResult{};
  };

  Config m_cfg;
  tbbWrap::task_arena m_taskArena;
  std::vector<std::shared_ptr<IContextDecorator>> m_decorators;
  std::vector<std::shared_ptr<IReader>> m_readers;
  std::vector<SequenceElementWithFpeResult> m_sequenceElements;
  std::unique_ptr<const Acts::Logger> m_logger;

  std::unordered_multimap<std::string, std::string> m_whiteboardObjectAliases;

  std::unordered_map<std::string, const DataHandleBase *> m_whiteBoardState;

  std::atomic<std::size_t> m_nUnmaskedFpe = 0;

  const Acts::Logger &logger() const { return *m_logger; }
};

std::ostream &operator<<(std::ostream &os,
                         const ActsExamples::Sequencer::FpeMask &m);

}  // namespace ActsExamples
