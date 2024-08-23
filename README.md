# SpSpSel
A toy project that aims at an introductory tutorial to spontaneous spectral spatial selectivity in MRI


This tutorial will develop itself with the following logic:

  0. Building infrastructure for the main context: Bloch equation solver and extended Bloch-McConnell equation. In the future we may add Hilbert space LvN equation solver and Liouvillian space LvN Eqn solver for more complicated simulations.
  1. Single active absoption with simple Bloch equation

    * Demonstrating basic spectral selectivity from Dante train with various modification of overall r.f. envelope.

    * Demonstrating basic spatial selectivity.

    * Demonstrating the phase refocusing in both spectral and spatial domain from slice rephasing. 

    * Demonstrating the spontaneous dual selectivity in both spectral and spatial domain

  2. Two active species with Bloch-McConnel Eqn

    - Repeating what we have done in Step 1

  3. Wrap the stepwise demo into module, then apply naive BlackBox optimization for fitting kinetic parameters.  
