======================
Oncotator
======================

The name of the directory, oncotator, is also the name of the **distribution**.
This distribution contains the oncotator package.

For more information:
http://www.broadinstitute.org/cancer/cga/oncotator

This distribution is the standalone version of Oncotator.  If you wish to use the web interface:
http://www.broadinstitute.org/oncotator

Please note that the web interface uses an older version of Oncotator and older datasources.

If you are a researcher/employee of the Broad Institute, you can find Broad-specific help in the Cancer Genome Analysis section of Confluence.

Installation
------------

Currently, Windows is unsupported, though this is due to a dependency, pysam, being unsupported in Windows.

IMPORTANT:  You will need root access to your python interpreter or a python virtual environment.  More information about virtual environments can be found on the following site:
https://pypi.python.org/pypi/virtualenv

As a reminder, virtualenv.py can be run as a standalone script, thereby bypassing superuser requirements.  Please see the above link for more details.

Before installing, we recommend installing pyvcf and numpy manually, before attempting the Oncotator install.  You may need to prepend each of the following commands with sudo::

    $ pip install numpy
    $ pip install pyvcf

This distribution is installable through the standard ``setup.py`` method.  Note that Distribute will be installed as part of the setup process if it isn't already::

    $ python setup.py install

Because the setup.py specifies an entry point as a console script, ``oncotator``  and ``initializeDatasource`` will be installed into your Python's ``bin/`` directory


Unit Tests
----------

NOTE: Unit tests require a minimum of 4GB to run.

Before running the unit tests, please perform the following steps:

1) Execute the following three lines in the same directory as setup.py::

    $ mkdir -p out
    $ ln -s test/configs configs
    $ ln -s test/testdata testdata

2) Many unit tests rely on having the standard set of hg19 datasources, which are in a separate download.  To point the unit testing framework to your datasources, you must create a personal test config::

    $ cp configs/personal-test.config.template configs/personal-test.config
    In configs/personal-test.config, replace ```dbDir=MY_DB_DIR/``` with ```dbDir=``` the appropriate path to you oncotator datasource directory.

Unit tests are run from the ``test`` target in setup.py::

    $ python setup.py test


Version Information
-------------------

Once Oncotator is installed, run it with the -V flag to get version information::

    $ oncotator -V
