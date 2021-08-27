// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Options/CommonOptions.hpp"

#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <exception>
#include <fstream>
#include <regex>
#include <system_error>

using namespace boost::program_options;

boost::program_options::options_description
ActsExamples::Options::makeDefaultOptions(std::string caption) {
  options_description opt(caption);

  opt.add_options()("help,h", "Produce help message");
  opt.add_options()(
      "loglevel,l", value<size_t>()->default_value(2),
      "The output log level. Please set the wished number (0 = VERBOSE, 1 = "
      "DEBUG, 2 = INFO, 3 = WARNING, 4 = ERROR, 5 = FATAL).");
  opt.add_options()(
      "response-file", value<std::string>()->default_value(""),
      "Configuration file (response file) replacing command line options.");

  return opt;
}

void ActsExamples::Options::addSequencerOptions(
    boost::program_options::options_description& opt) {
  // sequencer options
  opt.add_options()("events,n", value<size_t>(),
                    "The number of events to process. If not given, all "
                    "available events will be processed.")(
      "skip", value<size_t>()->default_value(0),
      "The number of events to skip")("jobs,j", value<int>()->default_value(-1),
                                      "Number of parallel jobs, negative for "
                                      "automatic.");
}

void ActsExamples::Options::addRandomNumbersOptions(
    boost::program_options::options_description& opt) {
  opt.add_options()("rnd-seed", value<uint64_t>()->default_value(1234567890u),
                    "Random numbers seed.");
}

void ActsExamples::Options::addGeometryOptions(
    boost::program_options::options_description& opt) {
  opt.add_options()("geo-surface-loglevel", value<size_t>()->default_value(3),
                    "The outoput log level for the surface building.")(
      "geo-layer-loglevel", value<size_t>()->default_value(3),
      "The output log level for the layer building.")("geo-volume-loglevel",
                                                      value<size_t>()
                                                          ->default_value(3),
                                                      "The output log level "
                                                      "for the volume "
                                                      "building.");
}

void ActsExamples::Options::addMaterialOptions(
    boost::program_options::options_description& opt) {
  opt.add_options()(
      "mat-input-type", value<std::string>()->default_value("build"),
      "The way material is loaded: 'none', 'build', 'proto', 'file'.")(
      "mat-input-file", value<std::string>()->default_value(""),
      "Name of the material map input file, supported: '.json', '.cbor' or "
      "'.root'.")("mat-output-file", value<std::string>()->default_value(""),
                  "Name of the material map output file (without extension).")(
      "mat-output-sensitives", value<bool>()->default_value(true),
      "Write material information of sensitive surfaces.")(
      "mat-output-approaches", value<bool>()->default_value(true),
      "Write material information of approach surfaces.")(
      "mat-output-representing", value<bool>()->default_value(true),
      "Write material information of representing surfaces.")(
      "mat-output-boundaries", value<bool>()->default_value(true),
      "Write material information of boundary surfaces.")(
      "mat-output-volumes", value<bool>()->default_value(true),
      "Write material information of volumes.")(
      "mat-output-dense-volumes", value<bool>()->default_value(false),
      "Write material information of dense volumes.")(
      "mat-output-allmaterial", value<bool>()->default_value(false),
      "Add protoMaterial to all surfaces and volume for the mapping.");
}

void ActsExamples::Options::addOutputOptions(
    boost::program_options::options_description& opt,
    OutputFormat formatFlags) {
  // Add specific options for this example
  opt.add_options()("output-dir", value<std::string>()->default_value(""),
                    "Output directory location.");

  if (ACTS_CHECK_BIT(formatFlags, OutputFormat::Root))
    opt.add_options()("output-root", bool_switch(),
                      "Switch on to write '.root' output file(s).");

  if (ACTS_CHECK_BIT(formatFlags, OutputFormat::Csv))
    opt.add_options()("output-csv", bool_switch(),
                      "Switch on to write '.csv' output file(s).");

  if (ACTS_CHECK_BIT(formatFlags, OutputFormat::Obj))
    opt.add_options()("output-obj", bool_switch(),
                      "Switch on to write '.obj' ouput file(s).");

  if (ACTS_CHECK_BIT(formatFlags, OutputFormat::Json))
    opt.add_options()("output-json", bool_switch(),
                      "Switch on to write '.json' ouput file(s).");

  if (ACTS_CHECK_BIT(formatFlags, OutputFormat::Cbor))
    opt.add_options()("output-cbor", bool_switch(),
                      "Switch on to write '.cbor' ouput file(s).");

  if (ACTS_CHECK_BIT(formatFlags, OutputFormat::Txt))
    opt.add_options()("output-txt", bool_switch(),
                      "Switch on to write '.txt' ouput file(s).");
}

void ActsExamples::Options::addInputOptions(
    boost::program_options::options_description& opt) {
  // Add specific options for this example
  opt.add_options()("input-dir", value<std::string>()->default_value(""),
                    "Input directory location.")(
      "input-files", value<std::vector<std::string>>(),
      "Input files, can occur multiple times.")(
      "input-root", value<bool>()->default_value(false),
      "Switch on to read '.root' file(s).")(
      "input-csv", value<bool>()->default_value(false),
      "Switch on to read '.csv' file(s).")("input-obj",
                                           value<bool>()->default_value(false),
                                           "Switch on to read '.obj' file(s).")(
      "input-json", value<bool>()->default_value(false),
      "Switch on to read '.json' file(s).")(
      "input-cbor", value<bool>()->default_value(false),
      "Switch on to read '.cbor' file(s).");
}

boost::program_options::variables_map ActsExamples::Options::parse(
    const boost::program_options::options_description& opt, int argc,
    char* argv[]) noexcept(false) {
  variables_map vm;
  store(command_line_parser(argc, argv).options(opt).run(), vm);
  notify(vm);

  if (vm.count("response-file") and
      not vm["response-file"].template as<std::string>().empty()) {
    // Load the file and tokenize it
    std::ifstream ifs(vm["response-file"].as<std::string>().c_str());
    if (!ifs) {
      throw(std::system_error(std::error_code(),
                              "Could not open response file."));
    }
    // Read the whole file into a string
    std::stringstream ss;
    ss << ifs.rdbuf();
    std::string rString = ss.str();
    std::vector<std::string> args;
    const std::regex rgx("[ \t\r\n\f]");
    std::sregex_token_iterator iter(rString.begin(), rString.end(), rgx, -1);
    std::sregex_token_iterator end;
    for (; iter != end; ++iter) {
      if (std::string(*iter).empty()) {
        continue;
      }
      args.push_back(*iter);
    }
    // Parse the file and store the options
    store(command_line_parser(args).options(opt).run(), vm);
  }

  // Automatically handle help
  if (vm.count("help")) {
    std::cout << opt << std::endl;
    vm.clear();
  }
  return vm;
}

Acts::Logging::Level ActsExamples::Options::readLogLevel(
    const boost::program_options::variables_map& vm) {
  return Acts::Logging::Level(vm["loglevel"].as<size_t>());
}

ActsExamples::Sequencer::Config ActsExamples::Options::readSequencerConfig(
    const boost::program_options::variables_map& vm) {
  Sequencer::Config cfg;
  cfg.skip = vm["skip"].as<size_t>();
  if (not vm["events"].empty()) {
    cfg.events = vm["events"].as<size_t>();
  }
  cfg.logLevel = readLogLevel(vm);
  cfg.numThreads = vm["jobs"].as<int>();
  if (not vm["output-dir"].empty()) {
    cfg.outputDir = vm["output-dir"].as<std::string>();
  }
  return cfg;
}

// Read the random numbers config.
ActsExamples::RandomNumbers::Config
ActsExamples::Options::readRandomNumbersConfig(
    const boost::program_options::variables_map& vm) {
  ActsExamples::RandomNumbers::Config cfg;
  cfg.seed = vm["rnd-seed"].as<uint64_t>();
  return cfg;
}
