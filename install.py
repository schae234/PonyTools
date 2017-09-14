#!/usr/bin/env bash
# Script to setup stuff in an anaconda environment
# Written by Joe Jeffers and Rob Schaefer
# Email: jeffe174@umn.edu, schae234@umn.edu

NAME="ponytools"

function usage(){
cat <<EOF
    Usage: $0 [flags]

    Flags
    -----
    -h | --help
        print help message
    -b | --base
        Base installation directory for $NAME (default: ~/.$NAME).
EOF
exit 0
}

# Configurable variables
BASE=$HOME/.$NAME

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

function red(){
    printf "${RED}$1${NC}\n"
}
function green(){
    printf "${GREEN}$1${NC}\n"
}

while [[ $# > 0 ]]
do
key="$1"
case $key in 
    -h|--help)
    usage
    shift
    ;;
    -b|--base)
    BASE=$2
    shift
    ;;
    *)  
        #unknown options
    ;;
esac
shift
done


#===================================================
#----------Check for internet connection------------
#===================================================
wget -q --tries=10 --timeout=20 --spider http://github.com
if [[ $? -ne 0 ]]; then
        echo "Check your internet connection and try again."
		exit
fi

export CWD=$(pwd)
#===================================================
#----------Setup the build Environment--------------
#===================================================
echo "Setting up the build environment"
source $HOME/.bashrc
mkdir -p $BASE
mkdir -p $BASE/conda
mkdir -p $BASE/bin 
mkdir -p $BASE/lib 

cd $BASE

export LD_LIBRARY_PATH=$BASE/lib:$LD_LIBRARY_PATH
export PATH=$PATH:$BASE/bin:$BASE/conda/bin

#===================================================
#----------------Install conda ---------------------
#===================================================
# if ! hash conda 2>/dev/null
if [ ! -e $BASE/conda/bin/conda ]
then
    cd $BASE
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
	bash miniconda.sh -b -f -p $BASE/conda
	rm -f miniconda.sh
else
    green "Conda Already Installed" 
fi

#===================================================
#----------Build the Conda Environment--------------
#===================================================
if [ ! -d $BASE/conda/envs/$NAME ]
then
    echo "Making the conda virtual environment named $NAME in $BASE"
    conda remove -y --name $NAME --all
    conda config --add envs_dirs $BASE/conda/envs
    conda config --append channels conda-forge
    conda config --append channels blaze
    conda create -y -n $NAME python=3 setuptools pip cython numpy scipy pandas \
        matplotlib ipython ipdb pytest-cov tqdm  
else
    green 'conda already installed'
fi

#===================================================
#----------Activate the Conda Environment-----------
#===================================================
green "activating $NAME"
source $BASE/conda/bin/activate $NAME
green 'checking python'
which python

pip install pysam

#==================================================
#-----------------Install Camoco-------------------
#=================================================
green "Installing $NAME"
cd $CWD
python setup.py install
python -c 'import ponytools'
if [ $? -eq 1 ]
then
    red '$NAME failed to install!'
    exit 1
else
    green '$NAME installed!'
fi
source deactivate 


