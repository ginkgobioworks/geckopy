Contributing
============

Anyone should feel free to make a pull request and contribute to the geckopy codebase. This document gives an overview of the repository structure, how the test harness works, rules for contributing, and information about the CI/CD pipeline in place.

Opening Issues
--------------

To open an issue, go to the `GitHub repo`_ and go to the *Issues* tab. There, click the **New Issue** button and post the feature or bug to be worked on. We ask that, for bug fixes, you enter in the title "[Bug] <Brief description goes here...>", and for features, you title the issues "[Title] <Brief description goes here...>". If an issue falls under neither category, do not add the "[bracket]" tag, and they'll be handled as miscellaneous issues.

Pull requests
-------------
A dev guide is available at the documentation https://geckopy.readthedocs.io/.

1. Fork/branch the master branch.

.. code:: bash

   git clone https://github.com/ginkgobioworks/geckopy.git
   cd geckopy
   git checkout -b 'feat-incredible-new-feature'
   
2. Install the required dependencies for geckopy development:

.. code:: bash

  pip install ".[dev]"

3. Please use incremental `semantic commits`_. When possible, commits should 
   also be atomic; i.e., each commit should pass the tests.

4. The pull requests must pass the CI/CD pipeline, which is equivalent to

.. code:: bash

  # PEP style
  flake8 geckopy tests
  # import sorting
  isort geckopy tests
  # code format
  black geckopy tests
  # CI test suite (at tests/)
  python -m pytest

For pull requests, go to the *Pull Request* tab on the same repo. Select the 
branch that you created after clicking *New Pull Request* and provide a title 
and description for the request. The status of your PR in tests will also be 
available when you go through this process. Approval must be given from one of 
the administrators for a PR to be accepted.

Testing
-------
The ``tests`` directory contains the tests run by pytest_. Briefly, ``pytest``
gathers the different inputs (a.k.a. fixtures_) for the tests at ``tests/conftest.py``. The input files involved should be at ``tests/data/``. After gathering the
fixtures, pytest will execute every function that start with ``tests_`` for each
file that start with ``test_``.


Documentation
-------------

Documentation is stored in the ``docs`` directory, and it is used to generate documentation on https://geckopy.readthedocs.io/. Builds for the `Read the Docs`_ site are triggered with every release.

Release (for maintainers)
-------------------------

Pypi releases are synchronized with git tags. For instance:

.. code:: shell

  git tag -a 0.2.5 -m 'Relevant message for this release'
  git push origin 0.2.5

For a tag to trigger a release, the CI must pass too!

.. _semantic commits: https://sparkbox.com/foundry/semantic_commit_messages
.. _GitHub repo: https://github.com/ginkgobioworks/geckopy
.. _pytest: https://docs.pytest.org/ 
.. _fixtures: https://docs.pytest.org/en/stable/fixture.html 
.. _Read the Docs: https://readthedocs.org/
