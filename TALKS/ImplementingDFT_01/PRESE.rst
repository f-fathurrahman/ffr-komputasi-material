:title: Density Functional Theory Calculations using Finite Difference
:data-transition-duration: 1500
:css: PRESE.css

----

Density Functional Theory Calculations using Finite Difference
==============================================================

Fadjar Fathurrahman
-------------------

----

What is density functional theory (DFT)
=======================================

Density functional theory (DFT) may refer to a theory or a method to investigate
many-body systems. The systems can be nucleons in an atomic nuclei, or electrons
in atoms, molecules, and solids. We will exclusively talk about electrons.


----

Many-body Schroedinger equation
===============================

Time-independent many-body Schroedinger equation:

.. math::

    \hat{H}\,\Psi\left(\{ \mathbf{r}_{i} \}\right) =
    E\,\Psi\left(\{ \mathbf{r}_{i} \}\right)

The many-body wave function :math:`\Psi\left(\{ \mathbf{r}_{i} \}\right)` is
very complicated quantity. In DFT we replace the role of many-body wave function
with (single-particle) electron density.

.. math::

    \rho(r) = \sum_{i} \psi_{i}(\mathbf{r})


----

Kohn-Sham equations
===================

Nowadays, most DFT calculations are carry out via solving the Kohn-Sham equations:

.. math::

    \hat{H}_{\mathrm{KS}} \psi_{i}(\mathbf{r}) = \epsilon_{i} \psi_{i}(\mathbf{r})

where the Kohn-Sham Hamiltonian :math:`\hat{H}_{\mathrm{KS}}` is defined as:

.. math::

    \hat{H}_{\mathrm{KS}} = -\frac{1}{2}\nabla^2 + V_{\mathrm{KS}}(\mathbf{r})

and the Kohn-Sham potential :math:`V_{\mathrm{KS}}(\mathbf{r})` consists of three terms:

.. math::

    V_{\mathrm{KS}}(\mathbf{r}) = V_{\mathrm{ion}}(\mathbf{r}) +
    V_{\mathrm{Ha}}(\mathbf{r}) + V_{\mathrm{xc}}(\mathbf{r})

----

The Kohn-Sham potential
=======================

About V-ion, V-Ha, and V-xc


----

Total energy
============

Kohn-Sham total energy.


----

Solutions to the Kohn-Sham equations
====================================

Self-consistent

----

Objective and approach
=======================

- 1d Schroedinger equation
- 2d Schroedinger equation
- 3d Schroedinger equation
- Poisson equation
- XC energy and potential calculation

----

Test math
=========

.. math::

    x & = 2\alpha \\
    y & = 3 + 2\beta

Test braket:

.. math::

    \left\langle \nabla^2 \Psi \right\rangle

    \left\langle \frac{1}{2} \middle| 1 \right\rangle

----

Test HTML
=========

.. raw:: html

    <h1>This is header 1</h1>

    <p style="color: blue;">This is a paragraph</p>

    <strong>I am stronk</strong>

----

Test HTML5 Canvas
=================

.. raw:: html

    <canvas id="myCanvas" width="200" height="100" style="border:1px solid #000000;">
    </canvas>
