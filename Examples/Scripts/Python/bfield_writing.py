#!/usr/bin/env python3

import acts
import acts.examples

u = acts.UnitConstants


solenoid = acts.SolenoidBField(
    radius=1200 * u.mm, length=6000 * u.mm, bMagCenter=2 * u.T, nCoils=1194
)
field = acts.solenoidFieldMap(
    rlim=(0, 1200 * u.mm),
    zlim=(-5000 * u.mm, 5000 * u.mm),
    nbins=(10, 10),
    field=solenoid,
)

print("Solenoid ready")

cfg = acts.examples.RootBFieldWriter.Config()
cfg.bField = field
cfg.gridType = acts.examples.RootBFieldWriter.GridType.rz
cfg.fileName = "solenoid.root"
cfg.treeName = "solenoid"

acts.examples.RootBFieldWriter.run(cfg, acts.logging.VERBOSE)

print("Done")
