Run the alignment examples
===============================

Prerequisites
-------------

Acts must be build with activated examples, Pythia8 and alignment support
(``ACTS_BUILD_EXAMPLES_PYTHIA8=on`` and ``ACTS_BUILD_ALIGNMENT=on``) to enable the fast simulation. ``<build>``
is used to identify the path to the build directory.

We assume that you have generated a simulation dataset based on the TrackML detector as described in :ref:`simulate-TrackML`. An example dataset for running the alignment would be e.g. a single muon sample in a 2T magnetic field. Suppose the generated single muon sample is available at ``data/sim_trackML/single_muon``. 

Run track fitting with misaligned detector 
----------------------

The effects of misalignment on the track parameters estimation can be mimicked by decorating some misalignment to the detector. The AlignedDetector and PayloaDetector are two 'misaligned' versions of the TrackML detector. The difference between them is that they have different handling of the alignment parameters for the detector modules. In addition to the options as used for truth fitting in in :ref:`truth-fit-TrackML`, additional options for misalignment decoration must be specified to run the track fiting with the AlignedDetector: 

.. code-block:: console

   $ <build>/bin/ActsExampleDetectorAlignContextual \
       --input-dir=data/sim_trackML/single_muon \
       --digi-config-file <source>/Examples/Algorithms/Digitization/share/default-smearing-config-generic.json \
       --bf-constant-tesla=0:0:2 \
       --align-sigma-iplane 100 \ 
       --align-sigma-oplane 0 \
       --align-sigma-irot 10 \
       --align-sigma-orot 0 \
       --output-dir=data/reco_misalignedTrackML_nonaligned/single_muon


The ``--digi-config-file`` specifies the path for the digitization configuration file. The magnetic field setup should be consistent between simulation and the track fitting. 
The ``--align-sigma-iplane`` and ``--align-sigma-oplane`` specifies the standard deviation of the shifts of the detector element in um within the plane and out of plane, respectively. 
The ``--align-sigma-oplane`` and ``--align-sigma-orot`` specifies the standard deviation of the rotation of the detector element in mrad around the direction normal and parallel to the plane, respectively.


Run detector alignment and refit with misalignment corrected 
----------------------

Alignment can be run to estimate the misalignment, which can then be further used to correct the alignment parameters of the detector elements. To run track fitting with some misalignment decorated as above and also corrected by runing alignment prior to it:

.. code-block:: console

   $ <build>/bin/ActsExampleDetectorAlignContextual \
       --input-dir=data/sim_trackML/single_muon \
       --digi-config-file <source>/Examples/Algorithms/Digitization/share/default-smearing-config-generic.json \
       --bf-constant-tesla=0:0:2 \
       --align-sigma-iplane 100 \
       --align-sigma-oplane 0 \
       --align-sigma-irot 10 \
       --align-sigma-orot 0 \
       --reco-with-misalignment-correction true \
       --alignment-geo-config-file <source>/Examples/Algorithms/Alignment/share/alignment-geo-contextualDetector.json \
       --output-dir=data/reco_misalignedTrackML_aligned/single_muon


The ``--reco-with-misalignment-correction`` must be true to turn on the alignment, and the ``--alignment-geo-config-file`` is a jason file to specify which detector objects are to be aligned. Currently, only module level alignment is possible. 
