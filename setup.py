import os
from setuptools import setup, find_packages


with open(os.path.join(os.path.dirname(__file__), 'README.rst')) as readme:
    long_description = readme.read()

# Allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

setup(
    name='staNMF',
    version='1.1',
    packages=['staNMF'],
    include_package_data=True,
    license='LICENSE.txt',
    description='python 2.7 implementation of stability NMF (Siqi Wu 2016)',
    long_description=long_description,
    url='https://github.com/greenelab/staNMF',
    author='Greene Lab',
    author_email='team@greenelab.com',
    classifiers=[
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    install_requires=['numpy', 'pandas', 'scipy', 'matplotlib'],
)
