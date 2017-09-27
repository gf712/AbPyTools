from setuptools import setup

about = {}
with open('./abpytools/__about__.py', 'r') as f:
    exec(f.read(), about)

setup(
    name='AbPyTools',
    classifiers=['Development Status :: 2 - Pre-Alpha',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: MIT License',
                 'Natural Language :: English',
                 'Programming Language :: Python :: 3 :: Only',
                 'Topic :: Scientific/Engineering :: Bio-Informatics'],
    keywords='antibody-analysis bioinformatics data-processing',
    version=about['__version__'],
    package_dir={'abpytools': 'abpytools'},
    packages=['abpytools',
              'abpytools.utils',
              'abpytools.core',
              'abpytools.analysis',
              'abpytools.features'],
    package_data={'abpytools': ['data/*.json']},
    url='https://github.com/gf712/AbPyTools',
    license=about['__license__'],
    author=about['__author__'],
    author_email=about['__author_email__'],
    description='Python package for antibody analysis',
    install_requires=['numpy',
                      'joblib',
                      'pandas',
                      'tqdm>=4.15',
                      'seaborn',
                      'matplotlib',
                      'beautifulsoup4',
                      'lxml',
                      'scipy'],
    test_suite="tests"
)
