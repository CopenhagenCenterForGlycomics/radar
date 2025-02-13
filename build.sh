#!/bin/bash

apt-get install -y build-essential autoconf python3-setuptools python3-setuptools-whl cython3 python3-wheel wget python2-dev python2-pip python2-setuptools

pip2 install Cython

(
	cd fasta2
	make lfasta
	mv lfasta /usr/local/bin/
)

autoreconf -i
./configure
python3 setup.py sdist bdist_wheel
python2 setup.py install