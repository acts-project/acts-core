# Building with Spack

[Spack](https://spack.io/) is a package manager designed to handle software
used in scientific computing, including computing in high-energy physics. The
Spack package manager can be used to more easily install ACTS' dependencies.
Builds using Spack are not officially supported. In order to start using Spack,
follow the instructions on the [official Spack
documentation](https://spack.readthedocs.io/en/latest/getting_started.html).

## Installing dependencies

Once Spack has been installed, we are ready to create a development
environment; for example:

```console
$ spack env create acts
```

This new environment can then be activated and deactivated as follows:

```console
$ spack env activate acts
$ despacktivate
```

:::{tip}
The creation of a Spack environment is persistent, i.e. you should only need to
create your environment once. The activation of the environment, on the other
hand, affects only your current shell. As a result, you will need to activate
the environment every time you open a new shell.
:::

With the Spack environment activated, you are ready to install ACTS'
dependencies. You can either do this manually, or you can rely on Spack's
definition of how to do this. To use the Spack definition, execute the
following command (tuned for ACTS version 36.1.0; make sure you update the
version number to match the version of ACTS you want to build):

```console
$ spack install --add --only dependencies "acts@36.1.0"
```

This specification will suffice to build the core ACTS library. If you want to
build ACTS plugins, you will need to add additional _variants_; the following
is a spec which can be used to build ACTS and all its plugins:

```console
$ spack install --add --only dependencies "acts@36.1.0+tgeo+geomodel+podio+edm4hep+examples+odd+fatras+json+dd4hep+geant4+fatras_geant4+hepmc3+pythia8+python+svg+traccc"
```

The string starting with `acts` in the commands above is the so-called _spec_,
and it consists of three parts. The first part is the name of the package,
which is `acts`. The second part is the version of the package, in this case
`36.1.0`; Spack will pull in a set of dependencies suitable for that given
version, so it is recommended to tune the version number to match the version
of ACTS you want to build. The final part is a series of variants, marked with
the `+` symbol. The command above allows you to build virtually all components
of ACTS, but it also pulls in a large number of dependencies. You can remove
variants you do not need in order to speed up the Spack build.

Once Spack has concretized and installed the ACTS dependencies, you can build
ACTS as normal using CMake.

### Running DD4hep

In order to run code depending on DD4hep, you may need to run an additional
shell script when you activate your Spack environment. This shell script is
called `thisdd4hep_only.sh` and resides in the binary directory of your Spack
view. An easy way to activate this shell script is to run the following
command:

```console
$ source $(dirname $(which ddsim))/thisdd4hep_only.sh
```

## Build caches

Work on build caches for ACTS is in progress.

## Support

Spack support is experimental and not officially supported. Because ACTS'
dependencies are constantly changing, it is possible that Spack's database may
at times be out of date. If you encounter any problems using Spack to build
ACTS, unofficial support is available in the
[#Spack](https://mattermost.web.cern.ch/acts/channels/spack) channel on the
ACTS Mattermost.

## Troubleshooting

The following are some common problems you may encounter using Spack to build
ACTS.

### DD4hep not finding factories

If you run into an error like the following, DD4hep is failing to find the
object files required to parse some its required factories. Remember to
initialize DD4hep as described in ["Running DD4hep"](#running-dd4hep):

```text
dd4hep: Failed to locate plugin to interprete files of type "lccdd" - no factory:lccdd_XML_reader.
        No factory with name Create(lccdd_XML_reader) for type lccdd_XML_reader found.
        Please check library load path and/or plugin factory name.
```
