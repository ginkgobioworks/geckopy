============
Installation
============

Basic installation
==================

.. include:: ../README.rst
    :start-after: installation-start
    :end-before: installation-end

Setting up a virtual environment first
======================================

We highly recommended installing cameo inside a virtual environment (virtualenv_).
virtualenvwrapper_ tremendously simplifies using virtualenv_. Once you
installed virtualenv_ and virtualenvwrapper_, run

.. code-block:: shell

    mkvirtualenv gecko  # or whatever you'd like to call your virtual environment
    workon gecko

and then continue with the installation instructions described below.

.. code-block:: shell

    pip install geckopy

Soft dependencies
=================

Pytfa is required for running thermodynamics simulations.

.. code-block:: shell

    pip install geckopy[pytfa]

The required dependeciens for development can be installed with

.. code-block:: shell

    pip install geckopy[dev]


.. _virtualenvwrapper: https://pypi.python.org/pypi/virtualenvwrapper
.. _virtualenv: https://pypi.python.org/pypi/virtualenv
