# -*- coding: utf-8 -*-
""" setup.py ; setuptools control. """

from setuptools import setup

with open("readme","rb") as f:
   long_descr=f.read().decode("utf-8")

setup(
  name="bpb3",
  packages=[
   "bpb3",
   "bpb3.script",
   "bpb3.script.script_high" ,
   "bpb3.script.script_high.other"],
  entry_points={
    "console_scripts":  ['bpb3 = bpb3.bpb3:main'] },
  version= 1.0,
  install_requires=[
                     'setuptools',
                     'pandas',
                      'numpy',
                     'argparse',
                      'matplotlib',
                     'seaborn','datetime','scipy','numba'
                     ],
  description = "BayesPI-BAR in Python3 : bpb3 ",
  long_description=long_descr,
  include_package_data=True,
  author= "Kirill Batmanov and Junbai Wang",
  author_email= "junbai@gmail.com"
 )
