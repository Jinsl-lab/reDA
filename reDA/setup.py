import setuptools
from distutils.core import  setup

setuptools.setup(
    name='reDA',
    version='1.0.0',
    description='differential abundance testing on scATAC-seq data using random walk with restart',
    author='Zirui Chen',
    author_email='22S112105@stu.hit.edu.cn',
    packages=setuptools.find_packages(),
    license='GPL-3.0'
   )
