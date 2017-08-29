#from distutils.core import setup
from setuptools import setup

setup(
    name='AbPyTools',
    version='0.1',
    package_dir={'abpytools': 'abpytools'},
    packages=['abpytools',
              'abpytools.utils',
              'abpytools.core',
              'abpytools.analysis',
              'abpytools.features'],
    package_data={'abpytools': ['data/*.json']},
    url='https://github.com/gf712/AbPyTools',
    license='MIT',
    author='Gil Ferreira Hoben',
    author_email='gil.hoben.16@ucl.ac.uk',
    description='Python package for antibody analysis',
    install_requires=['numpy>=1.13',
                      'joblib>=0.11',
                      'pandas>=0.19',
                      'tqdm>=4.15',
                      'seaborn>=0.7',
                      'matplotlib>=1.5'],
    test_suite="tests"
)
