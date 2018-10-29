from setuptools import setup, Extension
import os
import sys
import shutil
from distutils.command.build_py import build_py as _build_py
from distutils.command.clean import clean as _clean
import re
import configparser

# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
import distutils.sysconfig
import subprocess

cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

config = configparser.ConfigParser()

####################################################################
#                           VERSION
####################################################################
about = {}
with open('abpytools/__about__.py', 'r') as f:
    exec(f.read(), about)

####################################################################
#                            CYTHON
####################################################################

try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext

    use_cython = True

    compile_time_env = {}


except ImportError:
    use_cython = False

    raise NotImplementedError("Currenly Cython is always required, but in future versions both "
                              "Python and Cython implementations will be available.")

if use_cython:

    cython_extensions = [Extension("abpytools.cython_extensions.convert_py_2_C",
                                   ["abpytools/cython_extensions/convert_py_2_C.pyx"],
                                   language='c++'),
                         Extension("abpytools.utils.math_utils",
                                   ["abpytools/utils/math_utils.pyx", "abpytools/utils/ops.cpp"],
                                   extra_compile_args=['-march=native', '-fopenmp', '-std=c++11'],
                                   extra_link_args=['-fopenmp'],
                                   language='c++'),
                         Extension("abpytools.analysis.distance_metrics_",
                                   ["abpytools/analysis/distance_metrics_.pyx"],
                                   language='c++')
                         ]

    # compile with SSE instructions?
    # TODO: check precise minimum SSE instruction set required
    try:
        call = subprocess.check_output("gcc -mavx2 -dM -E - < /dev/null | egrep 'SSE4_2'",
                                       stderr=subprocess.STDOUT, shell=True)
        compile_time_env['SSE4_2'] = 1
    except subprocess.CalledProcessError:
        compile_time_env['SSE4_2'] = 0

    cython_extensions_ = cythonize(cython_extensions,
                                   compile_time_env=compile_time_env,
                                   compiler_directives={'language_level': 3})


else:
    cython_extensions_ = None

####################################################################
#                           PROTOBUF
####################################################################
# adapted from https://github.com/protocolbuffers/protobuf/blob/master/python/setup.py
try:
    import google.protobuf as protobuf

    HAS_PROTOBUF = True
except:
    HAS_PROTOBUF = False

protoc = None


def get_protoc_version(path):
    """
    Get the protoc version
    Args:
        path: path to protoc

    Returns:
        protobuf version with the format "%d.%d.%d"
    """

    with subprocess.Popen(f"{path} --version", shell=True, stdout=subprocess.PIPE) as proc:
        # get stdout
        protobuf_version_b = proc.communicate()[0]
        protobuf_version_out = protobuf_version_b.decode("utf-8")
        protobuf_version_re = re.findall(r"\d.\d.\d", protobuf_version_out)
        if len(protobuf_version_re) == 1:
            protobuf_version = protobuf_version_re[0]
        else:
            protobuf_version = None
    return protobuf_version


if HAS_PROTOBUF:
    # Find the Protocol Compiler.
    if 'PROTOC' in os.environ and os.path.exists(os.environ['PROTOC']):
        candidate_protoc = os.environ['PROTOC']
        if get_protoc_version(candidate_protoc) == protobuf.__version__:
            protoc = candidate_protoc

    # look for protoc in anaconda first
    if protoc is None and os.path.exists('/'.join(sys.executable.split('/')[:-1] + ["protoc"])):
        candidate_protoc = '/'.join(sys.executable.split('/')[:-1] + ["protoc"])
        if get_protoc_version(candidate_protoc) == protobuf.__version__:
            protoc = candidate_protoc
    # look for protoc in /usr/local/bin
    if protoc is None and os.path.exists("/usr/local/bin/protoc"):
        candidate_protoc = "/usr/local/bin/protoc"
        if get_protoc_version(candidate_protoc) == protobuf.__version__:
            protoc = candidate_protoc
    if protoc is None:
        candidate_protoc = shutil.which("protoc")
        if get_protoc_version(candidate_protoc) == protobuf.__version__:
            protoc = candidate_protoc
    if protoc is None:
        print("-- Failed to find a matching proto compiler and protobuf python library!")
        print("-- Disabling protobuf support.")
        print(f"-- Protobuf python library version: {protobuf.__version__}")
        HAS_PROTOBUF = False
    else:
        print(f"-- Protoc  : {protoc} (Found {get_protoc_version(protoc)})\n"
              f"-- Protobuf: {protobuf.__file__} (Found {protobuf.__version__})")
        print("-- Found a matching proto compiler  and protobuf python library")
        HAS_PROTOBUF = True

if HAS_PROTOBUF:
    def generate_proto(source, require=True):
        """Invokes the Protocol Compiler to generate a _pb2.py from the given
        .proto file.  Does nothing if the output already exists and is newer than
        the input."""

        if not require and not os.path.exists(source):
            return

        output = source.replace(".proto", "_pb2.py").replace("../src/", "")

        if (not os.path.exists(output) or
                (os.path.exists(source) and
                 os.path.getmtime(source) > os.path.getmtime(output))):
            print("Generating %s..." % output)

            if not os.path.exists(source):
                sys.stderr.write("Can't find required file: %s\n" % source)
                sys.exit(-1)

            if protoc is None:
                sys.stderr.write(
                    "protoc is not installed nor found in ../src.  Please compile it "
                    "or install the binary package.\n")
                sys.exit(-1)

            protoc_command = [protoc, "-I.", "--python_out=.", source]
            if subprocess.call(protoc_command) != 0:
                sys.exit(-1)


    class build_py(_build_py):
        def run(self):
            # Generate necessary .proto file if it doesn't exist.
            generate_proto("abpytools/formats/chain.proto")
            generate_proto("abpytools/formats/fab.proto")

            _build_py.run(self)

else:
    build_py = _build_py

if HAS_PROTOBUF:
    config['PROTOBUF'] = {'PROTOBUF': 1,
                          'PROTOC_VERSION': get_protoc_version(protoc),
                          'PROTOBUF_LIB_VERSION': protobuf.__version__}
else:
    config['PROTOBUF'] = {'PROTOBUF': 0,
                          'PROTOC_VERSION': 0,
                          'PROTOBUF_LIB_VERSION': 0}

with open("abpytools/config.ini", 'w') as configfile:
    config.write(configfile)


####################################################################
#                           CLEANUP
####################################################################
# adapted from https://stackoverflow.com/questions/27843481/python-project-using-protocol-buffers-deployment-issues
class clean(_clean):
    def run(self):
        # Delete generated files in the code tree.
        for (dirpath, dirnames, filenames) in os.walk("."):
            for filename in filenames:
                filepath = os.path.join(dirpath, filename)
                if filepath.endswith("_pb2.py") or filepath.endswith("so") or filepath.endswith("pyc"):
                    os.remove(filepath)
                if filepath.endswith(".pyx"):
                    # print(filepath.replace('pyx', 'cpp'))
                    try:
                        os.remove(filepath.replace('pyx', 'cpp'))
                    except:
                        # cython generated file no longer exists
                        pass
        # _clean is an old-style class, so super() doesn't work.
        _clean.run(self)


####################################################################
#                           SETUP
####################################################################
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
              'abpytools.core.flags',
              'abpytools.formats',
              'abpytools.analysis',
              'abpytools.features',
              'abpytools.cython_extensions'],
    package_data={'abpytools': ['data/*.json',
                                'config.ini']},
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
    cmdclass={'build_ext': build_ext,
              'build_py': build_py,
              'clean': clean}
)
