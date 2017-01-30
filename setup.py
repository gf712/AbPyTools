from distutils.core import setup
# from setuptools import setup

setup(
    name='AbPyTools',
    version='0.1',
    package_dir={'abpytools': 'abpytools'},
    packages=['abpytools', 'abpytools.utils', 'abpytools.antibody', 'abpytools.analysis'],
    package_data={'abpytools': ['data/*.json']},
    url='https://github.com/gf712/AbPyTools',
    license='MIT',
    author='Gil Ferreira Hoben',
    author_email='gil.hoben.16@ucl.ac.uk',
    description='Python package for antibody analysis',
    # install_requires=['numpy', 'tqdm', 'joblib', 'seaborn', 'matplotlib']
)
