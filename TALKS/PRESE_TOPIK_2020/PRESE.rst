:title: Introduction to Research Topics (2020)
:data-transition-duration: 1500
:css: PRESE.css

----

Introduction to Research Topics (2020)
======================================

Solving Kohn-Sham Problem
-------------------------

Fadjar Fathurrahman

----

:data-x: -800
:data-y: -800
:data-scale: 0.1

Overview
========

- What is it

- Various approach to solve Kohn-Sham equation

- What are going to be done

----

:data-x: r200
:data-y: r0
:data-scale: 0.1

Density functional theory
=========================

KS equation arises in the density functional theory (DFT) of electronic structure of materials.
It is now commonly used as a foundation to explore various
properties of materials.

Many research topics in CMD Lab are using Kohn-Sham density functional
theory (KSDFT).

- Fuel cells and batteries

- Drugs

- :math:`\mathrm{CO}_{2}` capture

- surface science and heterogenous catalysis

Many concrete examples by Dr. Haris and Dr. Ganes.


----

Kohn-Sham energy functional
===========================

Minimize with respect to single-particle wave functions 
:math:`\{\psi_{i}(\mathbf{r})\}`

.. math::

    E\left[\{\psi_{i}(\mathbf{r})\}\right] = E_{\mathrm{kin}} + E_{\mathrm{ext}} + E_{\mathrm{Ha}} + E_{\mathrm{xc}}

where:

.. math::

    \begin{aligned}
    E_{\mathrm{kin}} & = -\frac{1}{2} \int \psi_{i}(\mathbf{r}) \nabla^{2} \psi_{i}(\mathbf{r})\,\mathrm{d}\mathbf{r} \\
    E_{\mathrm{ext}} & = \int \rho(\mathbf{r}) V_{\mathrm{ext}}(\mathbf{r})\,\mathrm{d}\mathbf{r} \\
    E_{\mathrm{Ha}} & = \frac{1}{2} \int \int \frac{\rho(\mathbf{r}) \rho(\mathbf{r}')}{\left|\mathbf{r}-\mathbf{r}'\right|}\, \mathrm{d}\mathbf{r}\,\mathrm{d}\mathbf{r}' \\
    \rho(\mathbf{r}) & = \sum_{i} f_{i}\, \psi_{i}^{*}(\mathbf{r}) \psi_{i}(\mathbf{r})
    \end{aligned}

----

Kohn-Sham equation
==================

Minimum of the Kohn-Sham energy functional also can be found by solving the Kohn-Sham
equation:

.. math::

    \hat{H}_{\mathrm{KS}}\,\psi_{i}(\mathbf{r}) & = \epsilon_{i}\,\psi_{i}(\mathbf{r})


Kohn-Sham Hamiltonian:

.. math::

    \begin{aligned}
    \hat{H}_{\mathrm{KS}} & = -\frac{1}{2}\nabla^{2} + V_{\mathrm{ext}}(\mathbf{r}) + V_{\mathrm{Ha}}(\mathbf{r}) + V_{\mathrm{xc}}(\mathbf{r}) \\
    V_{\mathrm{Ha}}(\mathbf{r}) & = \int \frac{\rho(\mathbf{r})}{\mathbf{r} - \mathbf{r}'}\,\mathrm{d}\mathbf{r}' \\
    V_{\mathrm{ext}}(\mathbf{r}) & = \sum_{I} \frac{Z_{I}}{\left|\mathbf{r} - \mathbf{R}_{I}\right|} \\
    V_{\mathrm{xc}}(\mathbf{r}) & = \frac{\delta E[\rho(\mathbf{r})]}{\delta \rho(\mathbf{r})}
    \end{aligned}


----

Using DFT in research
=====================

Nowadays we can solve the Kohn-Sham by using many software packages.

We need to learn the softwares:

- learn the tutorials and do the tutorials

- read the manual

- practice, practice, and practice


----

Common research steps in CMD
============================

**STEP 1** Choose the system

Preferably having significance or related to problems in industries, experiments, 
or of global importance. (energy, health, environment, etc)

**STEP 2** Model and calculate

Choose a software and learn how to use the software to calculate the
properties we want to investigate. Create models of the system and calculate
the properties we want to learn.

**STEP 3** Analyze the results.

Find the parameters (structure, types of atoms, etc)
of the system that strongly affect the properties of the system

**STEP 4** Optimize

Vary, optimize, or tune the parameters until we get the properties that
we want.

----

Problems
========

- Most interesting systems are quite difficult to model: too big for
  DFT, requires a long time for the calculation to finish.

- Competition from other researchers: The system we are studying is popular.
  There is a chance that there are already similar calculations being done with
  bigger size and more sophisticated methods.

- There is "disconnects" between the equations that we read in the books
  and the softwares. We are not really sure what are actually calculated.
  We rely on the softwares to deliver the results that we want.

----

My researches
=============

My current researches are centered around implementing and exploring
various approaches to solve KS equation, i.e. I choose to write my own KSDFT solver.

Motivation:

- Not satisfied with current available softwares

- Want to learn the inside of the black box

- Avoid big atomistic systems (reducing the load of computing facilities in ITB)

- Education: any details of solving KS problem are not well described in literatures, especially
  in terms of the actual code.

----

Some consequences
=================

- Need to learn how to read and write codes. A lot of codes.

- Need to pay more attention to the equations presented in the literature.

- Need to learn and implement many numerical algorithms. Some of them are quite
  advanced, such as nonlinear optimization and iterative diagonalization.

- Need more time investment. Deal with bugs in the codes

- Limited to small systems.

- More difficult to publish papers.

----

Several ongoing and planned
===========================

- Plane wave basis: PWDFT.jl

- Real space methods: ffr-LFDFT

- FLAPW: (planned)

- Spectral finite element: (planned)

----

About the inhouse softwares
===========================

About PWDFT.jl (https://github.com/f-fathurrahman/PWDFT.jl):

- Based on plane wave basis set and pseudopotentials.
- Similar to Quantum ESPRESSO, ABINIT, VASP, etc.
- Written in Julia programming language.

About ffr-LFDFT (https://github.com/f-fathurrahman/ffr-LFDFT):

- Based on Lagrange functions (LF)
- Not many similar programs yet ...
- Written in Fortran language


----

Some topics in PWDFT.jl
=======================

- Development of geometry optimization and molecular dynamics methods

- Improvement of Kohn-Sham solvers:

  - Advanced mixing methods: Anderson accelaration, preconditioning, etc

  - Direct minimization for metallic systems

- Advanced XC functionals: meta-GGA, SCAN, vdw-DF, etc.

- Exact-exchange calculations

- Parallelization: threads, MPI, and GPU (with CUDA)

----

Some topics in real-space methods
=================================

- Porting from Fortran to Julia (a book is being written, a draft is available for those who are
  interested).

- Similar to PWDFT.jl

- TDDFT

----

More about PWDFT.jl
===================

Not a program

like toolbox (MATLAB) or library/package (Python)

No need for input file, learn the API (application programming interface)

----

Example
=======

.. code:: julia
    
    using PWDFT
    # crystalline structure
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))
    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
    # Solve the SCF
    KS_solve_SCF!( Ham, mix_method="rpulay", betamix=0.5 )


----


Julia programming language
==========================

Hello I am efefer again.

.. code:: julia
    :class: hidden

    using PyPlot
    α = β + 2

----

Testing Raw HTML
================

.. raw:: html

    <strong>This is efefer</strong>

Canvas example

.. raw:: html

    <p style="text-align: center;">
    <canvas id="myCanvas" width="200" height="200" style="border:1px solid #000000;">
    </canvas>
    </p>
    <script>
    var c = document.getElementById("myCanvas");
    var ctx = c.getContext("2d");
    ctx.beginPath();
    ctx.arc(95, 50, 60, 0, 2*Math.PI);
    ctx.stroke(); 
    </script>
