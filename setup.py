#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
from shutil import rmtree

from setuptools import setup, Command

NAME = 'm4'
DESCRIPTION = 'A suite of tools to calibrate ELT M4 deformable mirror',
URL = 'https://github.com/ArcetriAdaptiveOptics/M4'
EMAIL = 'marco.xompero@inaf.it',
AUTHOR = 'Marco Xompero, Chiara Selmi',
LICENSE = 'MIT'
KEYWORDS = 'Adaptive Optics, Astrophysics, INAF, Arcetri',

here = os.path.abspath(os.path.dirname(__file__))
# Load the package's __version__.py module as a dictionary.
about = {}
with open(os.path.join(here, NAME, '__version__.py')) as f:
    exec(f.read(), about)


class UploadCommand(Command):
    """Support setup.py upload."""

    description = 'Build and publish the package.'
    user_options = []

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print('\033[1m{0}\033[0m'.format(s))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status('Removing previous builds…')
            rmtree(os.path.join(here, 'dist'))
        except OSError:
            pass

        self.status('Building Source and Wheel (universal) distribution…')
        os.system('{0} setup.py sdist bdist_wheel --universal'.format(sys.executable))

        self.status('Uploading the package to PyPI via Twine…')
        os.system('twine upload dist/*')

        self.status('Pushing git tags…')
        os.system('git tag v{0}'.format(about['__version__']))
        os.system('git push --tags')

        sys.exit()


setup(name=NAME,
      description=DESCRIPTION,
      version=about['__version__'],
      classifiers=['Programming Language :: Python :: 3',
                   ],
      long_description=open('README.md').read(),
      url=URL,
      author_email=EMAIL,
      author=AUTHOR,
      license=LICENSE,
      keywords=KEYWORDS,
      packages=['m4',
                'm4.analyzers',
                'm4.configuration',
                'm4.data',
                'm4.devices',
                'm4.ground',
                'm4.gui',
                'm4.mini_OTT',
                'm4.ott_sim',
                'm4.type',
                'm4.utils'
                ],
      install_requires=["numpy",
                        "arte",
                        "astropy",
                        "h5py",
                        "scikit-image>=0.19",
                        "IPython",
                        "opcua",
                        "mock",
                        "zmq",
                        "pyzabbix",
                        "guietta",
                        "pyyaml",
                        "serial",
                        "pandas",
                        "playsound==1.2.2",
                        "jdcal",
                        "matplotlib"
                        ],
      package_data={
          'm4': ['data/*'],
      },
      include_package_data=True,
            entry_points={
          'gui_scripts': [
              'ott_geometry=m4.ground.GUI:main',
          ],
      },
      test_suite='test',
      cmdclass={'upload': UploadCommand, },
      )
