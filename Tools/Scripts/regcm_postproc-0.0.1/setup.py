#!/usr/bin/env python

from distutils.core import setup

setup(name='regcm_postproc',
      version='0.0.1',
      author='Niane Aliou, Pattnayak Kanhu Charan, Semie Addisu, Sobolev Andrey',
      author_email='',
      packages=['regcm_postproc'],
      package_dir = {'regcm_postproc': 'src'},
      scripts=['src/app.py',]
      )