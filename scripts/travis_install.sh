#!/usr/bin/env bash
sudo apt-get update
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda info -a
conda create -q -n test-environment python=${PYTHON_VERSION} scipy numpy pip matplotlib pandas libprotobuf=3.6.0
source activate test-environment
pip install coveralls coverage cython protobuf==3.6.0
python setup.py install