#!/usr/bin/env bash
source activate test-environment
pip install sphinx recommonmark
python setup.py build_ext -i
python setup.py build
cd docs && make html && cd ..