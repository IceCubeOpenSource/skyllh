.. _installation:

************
Installation
************


Prerequisites
=============

The SkyLLH framework has several dependencies. They are listed in `requirements.txt <https://github.com/icecube/skyllh/blob/master/requirements.txt>`_ file:

::

    astropy
    numpy
    scipy
    iminuit
    matplotlib

They can be installed from `skyllh` directory with:

.. code:: bash

    pip install -r requirements.txt

On cobalt and NPX servers we can use CVMFS Python 3 virtual environment with all necessary packages already installed. In order to activate it run:

.. code:: bash

    eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`


Setup
=====

Using pip
---------

The latest `skyllh` release can be installed from the
`PyPI <https://pypi.org/project/skyllh/>`_ repository:

.. code:: bash

    pip install skyllh

The current development version can be installed using pip:

.. code:: bash

    pip install git+https://github.com/icecube/skyllh.git#egg=skyllh

Optionally, the editable package version with a specified reference can be
installed by:

.. code:: bash

    pip install -e git+https://github.com/icecube/skyllh.git@[ref]#egg=skyllh

where

* `-e` is the editable flag
* `[ref]` is an optional argument containing a specific commit hash, branch name
  or tag

Cloning from GitHub
-------------------

The framework is split into two packages:

1. `github.com/icecube/skyllh <https://github.com/icecube/skyllh>`_

  * Contains open source code with classes defining the detector independent
    likelihood framework.

2. `github.com/icecube/i3skyllh <https://github.com/icecube/i3skyllh>`_

  * Contains collections of pre-defined SkyLLH IceCube analyses and pre-defined
    IceCube datasets.

In order to set it up, we have to clone git repositories and add them to the
`PYTHONPATH`:

.. code:: bash

    git clone git@github.com:icecube/skyllh.git /path/to/skyllh
    git clone git@github.com:icecube/i3skyllh.git /path/to/i3skyllh
    export PYTHONPATH=$PYTHONPATH:/path/to/skyllh
    export PYTHONPATH=$PYTHONPATH:/path/to/i3skyllh

Alternatively, we can add them inside the python script:

.. code:: python

    import sys

    # Add the skyllh and i3skyllh packages to the PYTHONPATH.
    sys.path.insert(0, '/path/to/skyllh')
    sys.path.insert(0, '/path/to/i3skyllh')
