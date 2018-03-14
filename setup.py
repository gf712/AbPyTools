from setuptools import setup #, Extension

# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
import distutils.sysconfig

# cfg_vars = distutils.sysconfig.get_config_vars()
# for key, value in cfg_vars.items():
#     if type(value) == str:
#         cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

about = {}
with open('abpytools/__about__.py', 'r') as f:
    exec(f.read(), about)

# cython_extension = Extension("abpytools.Cextensions",
#                              ["./abpytools/cython_extensions/Cextensions.pyx"],
#                              extra_compile_args=['-O3'],
#                              language='c++')

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
    test_suite="tests",
    # ext_modules=[cython_extension],
    # cmdclass={'build_ext': build_ext}
)
