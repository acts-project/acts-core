#!/usr/bin/env python3
import os

import acts
import acts.examples

u = acts.UnitConstants


def runPropagation(trackingGeometry, field, outputDir, s=None):
    s = s or acts.examples.Sequencer(events=10)

    rnd = acts.examples.RandomNumbers(seed=42)

    nav = acts.Navigator(trackingGeometry=trackingGeometry)

    stepper = acts.EigenStepper(field)
    # stepper = acts.AtlasStepper(field)
    # stepper = acts.StraightLineStepper()

    print("We're running with:", type(stepper).__name__)
    prop = acts.examples.ConcretePropagator(acts.Propagator(stepper, nav))

    alg = acts.examples.PropagationAlgorithm(
        propagatorImpl=prop,
        level=acts.logging.INFO,
        randomNumberSvc=rnd,
        ntests=1000,
        sterileLogger=True,
        propagationStepCollection="propagation-steps",
    )

    s.addAlgorithm(alg)

    # Output
    s.addWriter(
        acts.examples.ObjPropagationStepsWriter(
            level=acts.logging.INFO,
            collection="propagation-steps",
            outputDir=outputDir + "/obj",
        )
    )

    s.addWriter(
        acts.examples.RootPropagationStepsWriter(
            level=acts.logging.INFO,
            collection="propagation-steps",
            filePath=outputDir + "/propagation_steps.root",
        )
    )

    return s


if "__main__" == __name__:
    matDeco = None
    # matDeco = acts.IMaterialDecorator.fromFile("material.json")
    # matDeco = acts.IMaterialDecorator.fromFile("material.root")

    (
        detector,
        trackingGeometry,
        contextDecorators,
    ) = acts.examples.GenericDetector.create(mdecorator=matDeco)

    # dd4hepCfg = acts.examples.DD4hepDetector.Config()
    # dd4hepCfg.xmlFileNames = ["thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
    # detector = acts.examples.DD4hepDetector()
    # trackingGeometry, contextDecorators = detector.finalize(dd4hepCfg, None)

    # field = acts.NullBField()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))

    # solenoid = acts.SolenoidBField(
    #     radius = 1200*u.mm,
    #     length = 6000*u.mm,
    #     bMagCenter = 2*u.T,
    #     nCoils = 1194
    # )
    # field = acts.solenoidFieldMap(
    #     rlim=(0, 1200*u.mm),
    #     zlim=(-5000*u.mm, 5000*u.mm),
    #     nbins=(50, 50),
    #     field=solenoid
    # )

    runPropagation(trackingGeometry, field, os.getcwd()).run()
