import os
from setuptools import setup, find_packages
import distutils.core

with open(os.path.join(os.path.dirname(__file__), 'README.txt')) as readme:
    long_description = readme.read()

distutils.core.setup(name='Distutils', long_description=long_description)

# Allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

setup(
      name='staNMF',
      version='0.1',
      packages=find_packages(),
      include_package_data=True,
      license='LICENSE.txt',
      description='python 2.7 implementation of stability NMF (Siqi Wu 2016)',
      long_description=long_description,
      url='https://github.com/greenelab/staNMF',
      author='Greene Lab',
      author_email='team@greenelab.com',
      install_requires=['numpy', 'pandas', 'scipy', 'matplotlib', 'spams'],
      classifiers=[
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
                  ],
      )
