#!/bin/bash

apt-get install -y buildtools-essential autoconf python3-setuptools-whl cython3 python3-wheel wget

wget 'https://fasta.bioch.virginia.edu/wrpearson/fasta/fasta2/fasta2.shar.Z'
gunzip fasta2.shar.Z

(
	cd fasta2
	/bin/sh fasta2.shar
	make lfasta
	mv lfasta /usr/local/bin/
)

autoreconf -i
./configure
python3 setup.py sdist bdist_wheel