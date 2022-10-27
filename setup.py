# -*- coding: utf-8 -*-
from setuptools import setup

package_dir = \
{'': 'src'}

packages = \
['ccs2mut']

package_data = \
{'': ['*']}

install_requires = \
['argparse>=1.4.0,<2.0.0',
 'biopython>=1.79,<2.0',
 'click>=7.0,<8.0',
 'cyvcf2>=0.30.16,<0.31.0',
 'natsort>=8.1.0,<9.0.0',
 'pysam>=0.19.1,<0.20.0']

extras_require = \
{':python_version >= "3.8" and python_full_version < "4.0.0"': ['numpy>=1.23.1,<2.0.0']}

entry_points = \
{'console_scripts': ['ccs2mut = ccs2mut.__main__:main']}

setup_kwargs = {
    'name': 'ccs2mut',
    'version': '0.0.1',
    'description': 'ccs2mut',
    'long_description': "ccs2mut\n=======\n\n|PyPI| |Python Version| |License|\n\n|Read the Docs| |Tests| |Codecov|\n\n|pre-commit| |Black|\n\n.. |PyPI| image:: https://img.shields.io/pypi/v/ccs2mut.svg\n   :target: https://pypi.org/project/ccs2mut/\n   :alt: PyPI\n.. |Python Version| image:: https://img.shields.io/pypi/pyversions/ccs2mut\n   :target: https://pypi.org/project/ccs2mut\n   :alt: Python Version\n.. |License| image:: https://img.shields.io/pypi/l/ccs2mut\n   :target: https://opensource.org/licenses/MIT\n   :alt: License\n.. |Read the Docs| image:: https://img.shields.io/readthedocs/ccs2mut/latest.svg?label=Read%20the%20Docs\n   :target: https://ccs2mut.readthedocs.io/\n   :alt: Read the documentation at https://ccs2mut.readthedocs.io/\n.. |Tests| image:: https://github.com/sjin09/ccs2mut/workflows/Tests/badge.svg\n   :target: https://github.com/sjin09/ccs2mut/actions?workflow=Tests\n   :alt: Tests\n.. |Codecov| image:: https://codecov.io/gh/sjin09/ccs2mut/branch/master/graph/badge.svg\n   :target: https://codecov.io/gh/sjin09/ccs2mut\n   :alt: Codecov\n.. |pre-commit| image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white\n   :target: https://github.com/pre-commit/pre-commit\n   :alt: pre-commit\n.. |Black| image:: https://img.shields.io/badge/code%20style-black-000000.svg\n   :target: https://github.com/psf/black\n   :alt: Black\n\n\nFeatures\n--------\n\n* TODO\n\n\nRequirements\n------------\n\n* TODO\n\n\nInstallation\n------------\n\nYou can install *ccs2mut* via pip_ from PyPI_:\n\n.. code:: console\n\n   $ pip install ccs2mut\n\n\nUsage\n-----\n\nPlease see the `Command-line Reference <Usage_>`_ for details.\n\n\nContributing\n------------\n\nContributions are very welcome.\nTo learn more, see the `Contributor Guide`_.\n\n\nLicense\n-------\n\nDistributed under the terms of the MIT_ license,\n*ccs2mut* is free and open source software.\n\n\nIssues\n------\n\nIf you encounter any problems,\nplease `file an issue`_ along with a detailed description.\n\n\nCredits\n-------\n\nThis project was generated from `@cjolowicz`_'s `Hypermodern Python Cookiecutter`_ template.\n\n\n.. _@cjolowicz: https://github.com/cjolowicz\n.. _Cookiecutter: https://github.com/audreyr/cookiecutter\n.. _MIT: http://opensource.org/licenses/MIT\n.. _PyPI: https://pypi.org/\n.. _Hypermodern Python Cookiecutter: https://github.com/cjolowicz/cookiecutter-hypermodern-python\n.. _file an issue: https://github.com/sjin09/ccs2mut/issues\n.. _pip: https://pip.pypa.io/\n.. github-only\n.. _Contributor Guide: CONTRIBUTING.rst\n.. _Usage: https://ccs2mut.readthedocs.io/en/latest/usage.html\n",
    'author': 'Sangjin Lee',
    'author_email': 'sl17@sanger.ac.uk',
    'maintainer': None,
    'maintainer_email': None,
    'url': 'https://github.com/sjin09/ccs2mut',
    'package_dir': package_dir,
    'packages': packages,
    'package_data': package_data,
    'install_requires': install_requires,
    'extras_require': extras_require,
    'entry_points': entry_points,
    'python_requires': '>=3.6.1,<4.0.0',
}


setup(**setup_kwargs)

