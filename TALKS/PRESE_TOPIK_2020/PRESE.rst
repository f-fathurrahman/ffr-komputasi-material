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

Overview
========

- What is it

- Various approach to solve Kohn-Sham equation

- What are going to be done

----

About KS equation
=================

KS equation arises in the density functional theory of electronic structure of materials.

It is now commonly used as a foundation to explore various
properties of materials.

- structure

- chemical reactivity

----

Kohn-Sham problems
==================

Minimize:

.. math::

    E\left[\{\psi_{i}(\mathbf{r})\}\right] = -\frac{1}{2} \int \psi_{i}(\mathbf{r}) \nabla^{2} \psi_{i}(\mathbf{r})\,\mathrm{d}\mathbf{r}

----

Test code
=========

Hello I am efefer again.

.. code:: julia
    :class: hidden

    using PyPlot
    α = β + 2


----

Testing math
============

Inline math: :math:`\alpha + \beta`

Display math:

.. math::

    e^{i \pi} + 1 = 0

    dS = \frac{dQ}{T}

That's good.