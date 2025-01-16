// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Framework/Sequencer.hpp"

#include "Acts/Plugins/FpeMonitoring/FpeMonitor.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/IWriter.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/SequenceElement.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <ostream>
#include <ratio>
#include <regex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <typeinfo>

#include <TROOT.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/core/demangle.hpp>
#include <boost/stacktrace/stacktrace.hpp>

namespace ActsExamples {

namespace {

std::string_view getAlgorithmType(const SequenceElement& element) {
  if (dynamic_cast<const IWriter*>(&element) != nullptr) {
    return "Writer";
  }
  if (dynamic_cast<const IReader*>(&element) != nullptr) {
    return "Reader";
  }
  return "Algorithm";
}

// Saturated addition that does not overflow and exceed
// std::numeric_limits<std::size_t>::max().
//
// From http://locklessinc.com/articles/sat_arithmetic/
std::size_t saturatedAdd(std::size_t a, std::size_t b) {
  std::size_t res = a + b;
  res |= -static_cast<int>(res < a);
  return res;
}

/// Shorten some common but lengthy C++ constructs
std::string demangleAndShorten(std::string name) {
  name = boost::core::demangle(name.c_str());

  // Remove std::allocator from vector
  const static std::regex vector_pattern(
      R"??(std::vector<(.*), std::allocator<(\1\s*)>\s*>)??");
  name = std::regex_replace(name, vector_pattern, "std::vector<$1>");

  // Shorten Acts::BoundVariantMeasurement
  const static std::regex variant_pattern(
      R"??(std::variant<(Acts::Measurement<Acts::BoundIndices, [0-9]ul>(,|)\s+)+>)??");
  name = std::regex_replace(name, variant_pattern,
                            "Acts::BoundVariantMeasurement");

  // strip namespaces
  boost::algorithm::replace_all(name, "std::", "");
  boost::algorithm::replace_all(name, "boost::container::", "");
  boost::algorithm::replace_all(name, "Acts::", "");
  boost::algorithm::replace_all(name, "ActsExamples::", "");
  boost::algorithm::replace_all(name, "ActsFatras::", "");

  return name;
}

}  // namespace

Sequencer::Sequencer(const Sequencer::Config& cfg)
    : m_cfg(cfg),
      m_taskArena((m_cfg.numThreads < 0) ? tbb::task_arena::automatic
                                         : m_cfg.numThreads),
      m_logger(Acts::getDefaultLogger("Sequencer", m_cfg.logLevel)) {
  if (m_cfg.numThreads == 1) {
    ACTS_INFO("Create Sequencer (single-threaded)");
  } else {
    ROOT::EnableThreadSafety();
    ACTS_INFO("Create Sequencer with " << m_cfg.numThreads << " threads");
  }

  const char* envvar = std::getenv("ACTS_SEQUENCER_DISABLE_FPEMON");
  if (envvar != nullptr) {
    ACTS_INFO(
        "Overriding FPE tracking Sequencer based on environment variable "
        "ACTS_SEQUENCER_DISABLE_FPEMON");
    m_cfg.trackFpes = false;
  }

  if (m_cfg.trackFpes && !m_cfg.fpeMasks.empty() &&
      !Acts::FpeMonitor::canSymbolize()) {
    ACTS_ERROR("FPE monitoring is enabled but symbolization is not available");
    throw std::runtime_error(
        "FPE monitoring is enabled but symbolization is not available");
  }
}

void Sequencer::addContextDecorator(
    std::shared_ptr<IContextDecorator> decorator) {
  if (!decorator) {
    throw std::invalid_argument("Can not add empty/NULL context decorator");
  }
  m_decorators.push_back(std::move(decorator));
  ACTS_INFO("Added context decorator '" << m_decorators.back()->name() << "'");
}

void Sequencer::addReader(std::shared_ptr<IReader> reader) {
  if (!reader) {
    throw std::invalid_argument("Can not add empty/NULL reader");
  }
  m_readers.push_back(reader);
  addElement(std::move(reader));
}

void Sequencer::addAlgorithm(std::shared_ptr<IAlgorithm> algorithm) {
  if (!algorithm) {
    throw std::invalid_argument("Can not add empty/NULL algorithm");
  }

  addElement(std::move(algorithm));
}

void Sequencer::addWriter(std::shared_ptr<IWriter> writer) {
  if (!writer) {
    throw std::invalid_argument("Can not add empty/NULL writer");
  }
  addElement(std::move(writer));
}

void Sequencer::addElement(const std::shared_ptr<SequenceElement>& element) {
  if (!element) {
    throw std::invalid_argument("Can not add empty/NULL element");
  }

  m_sequenceElements.push_back({element});

  std::string elementType{getAlgorithmType(*element)};
  std::string elementTypeCapitalized = elementType;
  elementTypeCapitalized[0] = std::toupper(elementTypeCapitalized[0]);
  ACTS_INFO("Add " << elementType << " '" << element->name() << "'");

  auto symbol = [&](const char* in) {
    std::string s = demangleAndShorten(in);
    std::size_t pos = 0;
    while (pos + 80 < s.size()) {
      ACTS_INFO("   " + s.substr(pos, pos + 80));
      pos += 80;
    }
    ACTS_INFO("   " + s.substr(pos));
  };

  bool valid = true;

  for (const auto* handle : element->readHandles()) {
    if (!handle->isInitialized()) {
      continue;
    }

    ACTS_INFO("<- " << handle->name() << " '" << handle->key() << "':");
    symbol(handle->typeInfo().name());

    if (auto it = m_whiteBoardState.find(handle->key());
        it != m_whiteBoardState.end()) {
      const auto& source = *it->second;
      if (!source.isCompatible(*handle)) {
        ACTS_ERROR("Adding "
                   << elementType << " " << element->name() << ":"
                   << "\n-> white board will contain key '" << handle->key()
                   << "'"
                   << "\nat this point in the sequence (source: "
                   << source.fullName() << "),"
                   << "\nbut the type will be\n"
                   << "'" << demangleAndShorten(source.typeInfo().name()) << "'"
                   << "\nand not\n"
                   << "'" << demangleAndShorten(handle->typeInfo().name())
                   << "'");
        valid = false;
      }
    } else {
      ACTS_ERROR("Adding " << elementType << " " << element->name() << ":"
                           << "\n-> white board will not contain key"
                           << " '" << handle->key()
                           << "' at this point in the sequence."
                           << "\n   Needed for read data handle '"
                           << handle->name() << "'");
      valid = false;
    }
  }

  if (valid) {  // only record outputs this if we're valid until here
    for (const auto* handle : element->writeHandles()) {
      if (!handle->isInitialized()) {
        continue;
      }

      ACTS_INFO("-> " << handle->name() << " '" << handle->key() << "':");
      symbol(handle->typeInfo().name());

      if (auto it = m_whiteBoardState.find(handle->key());
          it != m_whiteBoardState.end()) {
        const auto& source = *it->second;
        ACTS_ERROR("White board will already contain key '"
                   << handle->key() << "'. Source: '" << source.fullName()
                   << "' (cannot overwrite)");
        valid = false;
        break;
      }

      m_whiteBoardState.emplace(std::pair{handle->key(), handle});

      if (auto it = m_whiteboardObjectAliases.find(handle->key());
          it != m_whiteboardObjectAliases.end()) {
        ACTS_DEBUG("Key '" << handle->key() << "' aliased to '" << it->second
                           << "'");
        m_whiteBoardState[it->second] = handle;
      }
    }
  }

  if (!valid) {
    throw SequenceConfigurationException{};
  }
}

void Sequencer::addWhiteboardAlias(const std::string& aliasName,
                                   const std::string& objectName) {
  const auto range = m_whiteboardObjectAliases.equal_range(objectName);
  for (auto it = range.first; it != range.second; ++it) {
    const auto& [key, value] = *it;
    if (value == aliasName) {
      ACTS_INFO("Key '" << objectName << "' aliased to '" << aliasName
                        << "' already set");
      return;
    }
  }

  m_whiteboardObjectAliases.insert({objectName, aliasName});

  auto oit = m_whiteBoardState.find(objectName);
  if (oit == m_whiteBoardState.end()) {
    ACTS_ERROR("Key '" << objectName << "' does not exist");
    return;
  }

  ACTS_INFO("Key '" << objectName << "' aliased to '" << aliasName << "'");
  m_whiteBoardState[aliasName] = oit->second;
}

std::vector<std::string> Sequencer::listAlgorithmNames() const {
  std::vector<std::string> names;

  // WARNING this must be done in the same order as in the processing
  for (const auto& decorator : m_decorators) {
    names.push_back("Decorator:" + decorator->name());
  }
  for (const auto& [algorithm, fpe] : m_sequenceElements) {
    names.push_back(std::string(getAlgorithmType(*algorithm)) + ":" +
                    algorithm->name());
  }

  return names;
}

std::pair<std::size_t, std::size_t> Sequencer::determineEventsRange() const {
  constexpr auto kInvalidEventsRange =
      std::make_pair(std::numeric_limits<std::size_t>::max(),
                     std::numeric_limits<std::size_t>::max());

  // Note on skipping events:
  //
  // Previously, skipping events was only allowed when readers where
  // available, since only readers had a `.skip()` functionality. The
  // `.skip()` interface has been removed in favour of telling the readers the
  // event they are requested to read via the algorithm context. Skipping can
  // now also be used when no readers are configured, e.g. for generating only
  // a few specific events in a simulation setup.

  // determine intersection of event ranges available from readers
  std::size_t beg = 0u;
  std::size_t end = std::numeric_limits<std::size_t>::max();
  for (const auto& reader : m_readers) {
    auto available = reader->availableEvents();
    beg = std::max(beg, available.first);
    end = std::min(end, available.second);
  }

  // since we use event ranges (and not just num events) they might not
  // overlap
  if (end < beg) {
    ACTS_ERROR("Available events ranges from readers do not overlap (beg="
               << beg << ", end=" << end << ")");
    return kInvalidEventsRange;
  }
  // configured readers without available events makes no sense
  // TODO could there be a use-case for zero events? run only setup functions?
  if (beg == end) {
    ACTS_ERROR("No events available");
    return kInvalidEventsRange;
  }
  // trying to skip too many events must be an error
  if (end <= saturatedAdd(beg, m_cfg.skip)) {
    ACTS_ERROR("Less events available than requested to skip");
    return kInvalidEventsRange;
  }
  // events range was not defined by either the readers or user command line.
  if ((beg == 0u) && (end == std::numeric_limits<std::size_t>::max()) &&
      (!m_cfg.events.has_value())) {
    ACTS_ERROR("Could not determine number of events");
    return kInvalidEventsRange;
  }

  // take user selection into account
  auto begSelected = saturatedAdd(beg, m_cfg.skip);
  auto endSelected = end;
  if (m_cfg.events.has_value()) {
    auto endRequested = saturatedAdd(begSelected, m_cfg.events.value());
    endSelected = std::min(end, endRequested);
    if (end < endRequested) {
      ACTS_INFO("Restrict requested number of events to available ones");
    }
  }

  return {begSelected, endSelected};
}

// helpers for per-algorithm timing information
namespace {
using Clock = std::chrono::high_resolution_clock;
using Duration = Clock::duration;
using Timepoint = Clock::time_point;
using Seconds = std::chrono::duration<double>;
using NanoSeconds = std::chrono::duration<double, std::nano>;

// RAII-based stopwatch to time execution within a block
struct StopWatch {
  Timepoint start;
  Duration& store;

  explicit StopWatch(Duration& s) : start(Clock::now()), store(s) {}
  ~StopWatch() { store += Clock::now() - start; }
};

// Convert duration to a printable string w/ reasonable unit.
template <typename D>
inline std::string asString(D duration) {
  double ns = std::chrono::duration_cast<NanoSeconds>(duration).count();
  if (1e9 < std::abs(ns)) {
    return std::to_string(ns / 1e9) + " s";
  } else if (1e6 < std::abs(ns)) {
    return std::to_string(ns / 1e6) + " ms";
  } else if (1e3 < std::abs(ns)) {
    return std::to_string(ns / 1e3) + " us";
  } else {
    return std::to_string(ns) + " ns";
  }
}

// Convert duration scaled to one event to a printable string.
template <typename D>
inline std::string perEvent(D duration, std::size_t numEvents) {
  return asString(duration / numEvents) + "/event";
}

void storeTiming(const std::vector<std::string>& identifiers,
                 const std::vector<Duration>& durations, std::size_t numEvents,
                 const std::string& path) {
  std::ofstream file(path);

  file << "identifier,time_total_s,time_perevent_s\n";

  for (std::size_t i = 0; i < identifiers.size(); ++i) {
    const auto time_total_s =
        std::chrono::duration_cast<Seconds>(durations[i]).count();
    file << identifiers[i] << "," << time_total_s << ","
         << time_total_s / numEvents << "\n";
  }
  file << "\n";
}
}  // namespace

int Sequencer::run() {
  // measure overall wall clock
  Timepoint clockWallStart = Clock::now();
  // per-algorithm time measures
  std::vector<std::string> names = listAlgorithmNames();
  std::vector<Duration> clocksAlgorithms(names.size(), Duration::zero());
  tbbWrap::queuing_mutex clocksAlgorithmsMutex;

  // processing only works w/ a well-known number of events
  // error message is already handled by the helper function
  std::pair<std::size_t, std::size_t> eventsRange = determineEventsRange();
  if ((eventsRange.first == std::numeric_limits<std::size_t>::max()) &&
      (eventsRange.second == std::numeric_limits<std::size_t>::max())) {
    return EXIT_FAILURE;
  }

  ACTS_INFO("Processing events [" << eventsRange.first << ", "
                                  << eventsRange.second << ")");
  ACTS_INFO("Starting event loop with " << m_cfg.numThreads << " threads");
  ACTS_INFO("  " << m_decorators.size() << " context decorators");
  ACTS_INFO("  " << m_sequenceElements.size() << " sequence elements");

  std::size_t nWriters = 0;
  std::size_t nReaders = 0;
  std::size_t nAlgorithms = 0;
  for (const auto& [alg, fpe] : m_sequenceElements) {
    if (dynamic_cast<const IWriter*>(alg.get()) != nullptr) {
      nWriters++;
    } else if (dynamic_cast<const IReader*>(alg.get()) != nullptr) {
      nReaders++;
    } else if (dynamic_cast<const IAlgorithm*>(alg.get()) != nullptr) {
      nAlgorithms++;
    } else {
      throw std::runtime_error{"Unknown sequence element type"};
    }
  }

  ACTS_INFO("  " << nReaders << " readers");
  ACTS_INFO("  " << nAlgorithms << " algorithms");
  ACTS_INFO("  " << nWriters << " writers");

  ACTS_VERBOSE("Initialize sequence elements");
  for (auto& [alg, fpe] : m_sequenceElements) {
    ACTS_VERBOSE("Initialize " << getAlgorithmType(*alg) << ": "
                               << alg->name());
    if (alg->initialize() != ProcessCode::SUCCESS) {
      ACTS_FATAL("Failed to initialize " << getAlgorithmType(*alg) << ": "
                                         << alg->name());
      throw std::runtime_error("Failed to process event data");
    }
  }

  // execute the parallel event loop
  std::atomic<std::size_t> nProcessedEvents = 0;
  std::size_t nTotalEvents = eventsRange.second - eventsRange.first;
  m_taskArena.execute([&] {
    tbbWrap::parallel_for(
        tbb::blocked_range<std::size_t>(eventsRange.first, eventsRange.second),
        [&](const tbb::blocked_range<std::size_t>& r) {
          std::vector<Duration> localClocksAlgorithms(names.size(),
                                                      Duration::zero());

          for (std::size_t event = r.begin(); event != r.end(); ++event) {
            ACTS_DEBUG("start processing event " << event);
            m_cfg.iterationCallback();
            // Use per-event store
            WhiteBoard eventStore(
                Acts::getDefaultLogger("EventStore#" + std::to_string(event),
                                       m_cfg.logLevel),
                m_whiteboardObjectAliases);
            // If we ever wanted to run algorithms in parallel, this needs to
            // be changed to Algorithm context copies
            AlgorithmContext context(0, event, eventStore);
            std::size_t ialgo = 0;

            /// Decorate the context
            for (auto& cdr : m_decorators) {
              StopWatch sw(localClocksAlgorithms[ialgo++]);
              ACTS_VERBOSE("Execute context decorator: " << cdr->name());
              if (cdr->decorate(++context) != ProcessCode::SUCCESS) {
                throw std::runtime_error("Failed to decorate event context");
              }
            }

            ACTS_VERBOSE("Execute sequence elements");

            for (auto& [alg, fpe] : m_sequenceElements) {
              std::optional<Acts::FpeMonitor> mon;
              if (m_cfg.trackFpes) {
                mon.emplace();
                context.fpeMonitor = &mon.value();
              }
              StopWatch sw(localClocksAlgorithms[ialgo++]);
              ACTS_VERBOSE("Execute " << getAlgorithmType(*alg) << ": "
                                      << alg->name());
              if (alg->internalExecute(++context) != ProcessCode::SUCCESS) {
                ACTS_FATAL("Failed to execute " << getAlgorithmType(*alg)
                                                << ": " << alg->name());
                throw std::runtime_error("Failed to process event data");
              }

              if (mon) {
                auto& local = fpe.local();

                for (const auto& [count, type, st] :
                     mon->result().stackTraces()) {
                  auto [maskLoc, nMasked] = fpeMaskCount(*st, type);
                  if (nMasked < count) {
                    std::stringstream ss;
                    ss << "FPE of type " << type
                       << " exceeded configured per-event threshold of "
                       << nMasked << " (mask: " << maskLoc
                       << ") (seen: " << count << " FPEs)\n"
                       << Acts::FpeMonitor::stackTraceToString(
                              *st, m_cfg.fpeStackTraceLength);

                    m_nUnmaskedFpe += (count - nMasked);

                    if (m_cfg.failOnFirstFpe) {
                      ACTS_ERROR(ss.str());
                      local.merge(mon->result());  // merge so we get correct
                                                   // results after throwing
                      throw FpeFailure{ss.str()};
                    } else if (!local.contains(type, *st)) {
                      ACTS_INFO(ss.str());
                    }
                  }
                }

                local.merge(mon->result());
              }
              context.fpeMonitor = nullptr;
            }

            nProcessedEvents++;
            if (logger().level() <= Acts::Logging::DEBUG) {
              ACTS_DEBUG("finished event " << event);
            } else if (nTotalEvents <= 100) {
              ACTS_INFO("finished event " << event);
            } else if (nProcessedEvents % 100 == 0) {
              ACTS_INFO(nProcessedEvents << " / " << nTotalEvents
                                         << " events processed");
            }
          }

          // add timing info to global information
          {
            tbbWrap::queuing_mutex::scoped_lock lock(clocksAlgorithmsMutex);
            for (std::size_t i = 0; i < clocksAlgorithms.size(); ++i) {
              clocksAlgorithms[i] += localClocksAlgorithms[i];
            }
          }
        });
  });

  ACTS_VERBOSE("Finalize sequence elements");
  for (auto& [alg, fpe] : m_sequenceElements) {
    ACTS_VERBOSE("Finalize " << getAlgorithmType(*alg) << ": " << alg->name());
    if (alg->finalize() != ProcessCode::SUCCESS) {
      ACTS_FATAL("Failed to finalize " << getAlgorithmType(*alg) << ": "
                                       << alg->name());
      throw std::runtime_error("Failed to process event data");
    }
  }

  fpeReport();

  // summarize timing
  Duration totalWall = Clock::now() - clockWallStart;
  Duration totalReal = std::accumulate(
      clocksAlgorithms.begin(), clocksAlgorithms.end(), Duration::zero());
  std::size_t numEvents = eventsRange.second - eventsRange.first;
  ACTS_INFO("Processed " << numEvents << " events in " << asString(totalWall)
                         << " (wall clock)");
  ACTS_INFO("Average time per event: " << perEvent(totalReal, numEvents));
  ACTS_DEBUG("Average time per algorithm:");
  for (std::size_t i = 0; i < names.size(); ++i) {
    ACTS_DEBUG("  " << names[i] << ": "
                    << perEvent(clocksAlgorithms[i], numEvents));
  }

  if (!m_cfg.outputDir.empty()) {
    storeTiming(names, clocksAlgorithms, numEvents,
                joinPaths(m_cfg.outputDir, m_cfg.outputTimingFile));
  }

  if (m_nUnmaskedFpe > 0) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

void Sequencer::fpeReport() const {
  if (!m_cfg.trackFpes) {
    return;
  }

  for (auto& [alg, fpe] : m_sequenceElements) {
    auto merged = std::accumulate(
        fpe.begin(), fpe.end(), Acts::FpeMonitor::Result{},
        [](const auto& lhs, const auto& rhs) { return lhs.merged(rhs); });
    if (!merged.hasStackTraces()) {
      // no FPEs to report
      continue;
    }
    ACTS_INFO("-----------------------------------");
    ACTS_INFO("FPE summary for " << getAlgorithmType(*alg) << ": "
                                 << alg->name());
    ACTS_INFO("-----------------------------------");

    std::vector<std::reference_wrapper<const Acts::FpeMonitor::Result::FpeInfo>>
        sorted;
    std::transform(merged.stackTraces().begin(), merged.stackTraces().end(),
                   std::back_inserter(sorted),
                   [](const auto& f) -> const auto& { return f; });
    std::ranges::sort(sorted, std::greater{},
                      [](const auto& s) { return s.get().count; });

    std::vector<std::reference_wrapper<const Acts::FpeMonitor::Result::FpeInfo>>
        remaining;

    for (const auto& el : sorted) {
      const auto& [count, type, st] = el.get();
      auto [maskLoc, nMasked] = fpeMaskCount(*st, type);
      ACTS_INFO("- " << type << ": (" << count << " times) "
                     << (nMasked > 0 ? "[MASKED: " + std::to_string(nMasked) +
                                           " per event by " + maskLoc + "]"
                                     : "")
                     << "\n"
                     << Acts::FpeMonitor::stackTraceToString(
                            *st, m_cfg.fpeStackTraceLength));
    }
  }

  if (m_nUnmaskedFpe > 0) {
    ACTS_ERROR("Encountered " << m_nUnmaskedFpe << " unmasked FPEs");
  } else {
    ACTS_INFO("No unmasked FPEs encountered");
  }
}

std::pair<std::string, std::size_t> Sequencer::fpeMaskCount(
    const boost::stacktrace::stacktrace& st, Acts::FpeType type) const {
  for (const auto& frame : st) {
    std::string loc = Acts::FpeMonitor::getSourceLocation(frame);
    auto it = loc.find_last_of(':');
    std::string locFile = loc.substr(0, it);
    unsigned int locLine = std::stoi(loc.substr(it + 1));
    for (const auto& [file, lines, fType, count] : m_cfg.fpeMasks) {
      const auto [start, end] = lines;
      if (boost::algorithm::ends_with(locFile, file) &&
          (start <= locLine && locLine < end) && fType == type) {
        std::string ls = start + 1 == end ? std::to_string(start)
                                          : "(" + std::to_string(start) + ", " +
                                                std::to_string(end) + "]";
        return {file + ":" + ls, count};
      }
    }
  }
  return {"NONE", 0};
}

Acts::FpeMonitor::Result Sequencer::fpeResult() const {
  Acts::FpeMonitor::Result merged;
  for (auto& [alg, fpe] : m_sequenceElements) {
    merged.merge(std::accumulate(
        fpe.begin(), fpe.end(), Acts::FpeMonitor::Result{},
        [](const auto& lhs, const auto& rhs) { return lhs.merged(rhs); }));
  }
  return merged;
}

std::ostream& operator<<(std::ostream& os,
                         const ActsExamples::Sequencer::FpeMask& m) {
  os << "FpeMask(" << m.file << ":";

  if (m.lines.first + 1 == m.lines.second) {
    os << m.lines.first;
  } else {
    os << "(" << m.lines.first << ", " << m.lines.second << "]";
  }
  os << ", " << m.type << " <= " << m.count << ")";
  return os;
}

}  // namespace ActsExamples
