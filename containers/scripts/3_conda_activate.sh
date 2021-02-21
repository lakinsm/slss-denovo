#!/usr/bin/env bash
# written by oleg osipenko
# 10012020

# ----- activate Conda -----
source /usr/local/miniconda3/etc/profile.d/conda.sh  # This must be invoked to use conda activate in Singularity
conda update -y -n base -c defaults conda
conda create -y -n DENOVO python=3 numpy matplotlib
source /usr/local/miniconda3/etc/profile.d/conda.sh  # This must be invoked to use conda activate in Singularity
source activate DENOVO

# ----- install Python packages -----

cd /opt/biotools || exit
wget https://downloads.sourceforge.net/project/quast/quast-5.0.2.tar.gz
tar -xzf quast-5.0.2.tar.gz
cd quast-5.0.2 || exit
./setup.py install
cd ../

git clone https://github.com/fenderglass/Flye
cd Flye || exit
python setup.py install
cd ../
