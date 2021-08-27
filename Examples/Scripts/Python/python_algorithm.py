#!/usr/bin/env python3

import acts
import acts.examples


s = acts.examples.Sequencer(events=1000, numThreads=10, logLevel=acts.logging.VERBOSE)


class PyAlg(acts.examples.BareAlgorithm):
    def __init__(self, name, level):
        acts.examples.BareAlgorithm.__init__(self, name, level)

    # def name(self):
    #     return "PyAlg"

    def execute(self, context):
        return acts.examples.ProcessCode.SUCCESS


class PyAlg1(acts.examples.IAlgorithm):
    def name(self):
        return "PyAlg1"

    def execute(self, context):
        print(context.eventStore.exists("blubb"))
        return acts.examples.ProcessCode.SUCCESS


# def add():
#     alg = PyAlg(name="BA", level=acts.logging.INFO)
#     # s.addAlgorithm(PyAlg(name="BA", level=acts.logging.INFO))
#     s.addAlgorithm(alg)


# add()

# a1 = PyAlg1()
# print(a1.name())
s.addAlgorithm(PyAlg(name="blubb", level=acts.logging.INFO))

print("alg go")
s.run()
