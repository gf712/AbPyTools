from setuptools import setup
from distutils.extension import Extension


# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
import distutils.sysconfig

# cfg_vars = distutils.sysconfig.get_config_vars()
# for key, value in cfg_vars.items():
#     if type(value) == str:
#         cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

about = {}
with open('abpytools/__about__.py', 'r') as f:
    exec(f.read(), about)

try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext

    use_cython = True

except ImportError:
    use_cython = False

    raise NotImplementedError("Currenly Cython is always required, but in future versions both "
                              "Python and Cython implementations will be available.")

if use_cython:

    cython_extensions = [Extension("abpytools.cython_extensions.convert_py_2_C",
                                   ["abpytools/cython_extensions/convert_py_2_C.pyx"]),
                         Extension("abpytools.utils.math_utils_",
                                   ["abpytools/utils/math_utils_.pyx"]),
                         Extension("abpytools.analysis.distance_metrics_",
                                   ["abpytools/analysis/distance_metrics_.pyx"])
                         ]
    cython_extensions_ = cythonize(cython_extensions)
else:
    cython_extensions_ = None

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
    ext_modules=cython_extensions_,
    cmdclass={'build_ext': build_ext}
)
