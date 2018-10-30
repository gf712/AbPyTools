#!/usr/bin/env bash
pip install sphinx recommonmark
python setup.py build_ext -i
python setup.py build
cd docs && make html && cd ..