:title: Introduction to Research Topics (2020)
:data-transition-duration: 1500
:css: PRESE.css

----

CMD Lab Meeting, 12 July 2020

Introduction to Research Topics (2020)
======================================

Solving Kohn-Sham Problem
-------------------------

Fadjar Fathurrahman

.. raw:: html

   <br>
   <br>

Overview:

- What it is

- What have been done

- What are going to be done

----

:data-x: -800
:data-y: -800
:data-scale: 0.1

Density functional theory
=========================

Kohn-Sham problem or more specifically, Kohn-Sham equation arises
in the density functional theory (DFT) of electronic structure of materials.
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

:data-x: r200
:data-y: r0
:data-scale: 0.1


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

Extensive list of software packages can be found in
Wikipedia page:

.. raw:: html
    
    <iframe
    src="https://en.wikipedia.org/wiki/List_of_quantum_chemistry_and_solid-state_physics_software"
    height="250" width="1200">
    </iframe>


We need to learn the softwares:

- learn the tutorials and do the tutorials

- read the manual

- practice, practice, and practice


----

Common research steps
=====================

**STEP 1** Choose the system

Preferably having significance or related to problems in industries, experiments, 
or of global importance. (energy, health, environment, etc)

**STEP 2** Model and calculate

Create models of the system and calculate
the properties we want to learn. Choose a software and learn how to use the
software to calculate the properties. 

**STEP 3** Analyze the results.

Find the parameters (structure, types of atoms, etc)
of the system that strongly affect the properties of the system

**STEP 4** Optimize

Vary, optimize, or tune the parameters until we get the properties that
we want.

----

Advantages
==========

- No need to pay attention on the nitty-gritty details of solving the Kohn-Sham equations.
  No need to invent the wheel.

- Directly calculate the properties we are interested in.

- Focus on the physics or chemistry of materials

- Faster to get the results (and get the publications done)

----

Problems
========

- Most interesting systems are quite difficult to model: too big for
  DFT, requires a long time for the calculation to finish.

- Competition from other researchers: The system we are studying is popular.
  There is a chance that there are already similar calculations being done with
  bigger size and more sophisticated methods.

- Need to wait or rely on the "software developers" for new developments.

- There is "disconnects" between the equations that we read in the books
  and the softwares. We are not really sure what are actually calculated.
  We rely on the softwares to deliver the results that we want.

- Need to implement specialized/custom post-processing (if not available)

----

My research
===========

My current research is centered around implementing and exploring
various approaches to solve KS equation, i.e. I choose to write
my own KSDFT solver from scratch.

Motivation:

- Not satisfied with current available softwares

- Want to learn the inside of the black box

- Avoid big atomistic systems (reducing the load of computing facilities in ITB)

- Education: many details of solving KS problem are not well described in literatures, especially
  in terms of the actual code.

----

Some consequences
=================

- Need to learn how to read and write codes. A lot of codes.

- Need to pay more attention to the equations presented in the literature.

- Need to learn and implement many numerical algorithms. Some of them are quite
  advanced, such as nonlinear optimization and iterative diagonalization.

- Need more time investment. Deal with bugs in the codes ...

- Limited features. Lack of man power.

- More difficult to publish papers.


However I choose to do it anyway to catch up with other researchers.


----

Several ongoing and planned
===========================

- Plane wave basis: PWDFT.jl

- Real space methods: ffr-LFDFT

- FLAPW: (planned)

- Spectral finite element: (planned)

- Discontinuous Galerkin finite element: (planned)



----

The inhouse softwares
=====================

About PWDFT.jl (https://github.com/f-fathurrahman/PWDFT.jl):

- Based on plane wave basis set and pseudopotentials.
- Similar to Quantum ESPRESSO, ABINIT, VASP, etc.
- Written in Julia programming language.
- The method is standard in solid solid state physics. Useful for comparison
  with other methods.

About ffr-LFDFT (https://github.com/f-fathurrahman/ffr-LFDFT):

- Based on Lagrange functions (LF). Shares many similarities with finite difference method.
- Historically, I developed this before writing PWDFT.jl
- Not many similar programs yet ...
- Written in Fortran language


----

Some topics in PWDFT.jl
=======================

There are a lot of things that can be done:

- Implementation of geometry optimization and molecular dynamics methods

- Improvement of Kohn-Sham solvers:

  - Advanced mixing methods: Anderson accelaration, preconditioning, etc

  - Direct minimization for metallic systems

- Advanced XC functionals: exact-exchange (EXX), meta-GGA, SCAN, vdw-DF,
  double hybrids, etc.

- Parallelization: threads, MPI, and GPU (with CUDA)

- Extension to USPP and PAW


I recommend development of molecular dynamics and geometry optimization.
Forces are already there, you only need to implement the integration of
equation of motion (Born-Oppenheimer and/or Car-Parrinello).




----

Some topics in real-space methods
=================================

- Porting from Fortran to Julia (a book is being written, a draft is available for those who are
  interested).

- The method (finite difference) is more intuitive than plane wave basis.

- For better parallel scaling, real-space methods are very promising compared to more traditional
  approach such as PW or APW.

- Futher development is mostly similar to PWDFT.jl:

  - Parallelization: MPI via PETSc and SLEPc.

- A simple implementation for time-dependent density functional theory
  is being worked on.



----

More about PWDFT.jl
===================

- Implemented not as a program but more like a toolbox (MATLAB) or library package (Python).

- No need for input file, learn the API (application programming interface)

- Present features:
  
  - Plane wave with norm-conserving pseudopotentials (Goedecker-Teter-Huter)
  - KS solvers: SCF with density or potential mixing and direct minimization
    (for non-metallic systems)
  - LDA VWN and GGA PBE functionals
  - k-points sampling
  - spin-polarized systems
  - force calculation

- All quantities are directly accessible: Kohn-Sham orbitals, electron density, potentials.
  No need for special post-processing step codes to obtain them.

- Notable missing features: stress tensor calculation and parallelization


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

- A rather new programming language, first announced in 2012.

- Developed by: Jeff Bezanson, Alan Edelman, Stefan Karpinski and Viral B. Shah.

- The libraries are not as extensive as Python yet, but it is improving.

.. image:: images/julia-language-developers-mit-00_0.png
    :height: 426px
    :width: 644px

http://news.mit.edu/2018/julia-language-co-creators-win-james-wilkinson-prize-numerical-software-1226

----

Additional Topics
=================

Machine learning related (very hot topics right now, many research groups are developing):

- Development of neural network potential for atomistic systems (software: AMP).
  It will be used for molecular dynamics.

- Gaussian process regression for surface reactions (software: Catlearn). Vieri is
  doing it now.

I encourage you to (re)write the software rather than using it directly:

- They are mostly being developed recently.

- Not yet robust. Many bugs.

- Codes, algorithms, procedures (computational protocols) are very important,
  not just the result. Reproducibility is also very important.

----

General strategy that you can adopt
===================================

Learn the package.

Read original literatures, study the equations and their
physical meaning.

Find the parts of the code that are most interested.

Refactor the code.

Rewrite the software using you own style. 

Find bugs, simplify, and improve


----

:data-x: 2000
:data-y: 0
:data-scale: 1.0

.. raw:: html

    <h1 style="text-align: center;">Thank you for your attention</h1>

.. raw:: html

    <img src="images/questions_cat.jpg" alt="Cat Questions" style="height: 50%">

    
