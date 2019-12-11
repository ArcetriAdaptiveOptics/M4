#!/usr/bin/env python
from setuptools import setup


setup(name='m4',
      description='A suite of tools to calibrate ELT M4 deformable mirror',
      version='0.1',
      classifiers=['Programming Language :: Python :: 3',
                   ],
      long_description=open('README.md').read(),
      url='',
      author_email='chiara.selmi@inaf.it',
      author='Chiara Selmi',
      license='',
      keywords='adaptive optics',
      packages=['m4',
                'm4.utils'],
      install_requires=["numpy",
                        "astropy"],
      test_suite='test',
      )
