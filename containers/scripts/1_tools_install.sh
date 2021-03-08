#!/usr/bin/env bash
# written by oleg osipenko
# 10012020
#

# ----- array of needed apps -----
packages[0]='build-essential'
packages[1]='git'
packages[2]='libz-dev'
packages[3]='libbz2-dev'
packages[4]='zlib1g-dev'
packages[5]='libncurses5-dev'
packages[6]='libncursesw5-dev'
packages[7]='liblzma-dev'
packages[8]='wget'
packages[9]='gfortran'
packages[10]='pkg-config'
packages[11]='libfreetype6-dev'
packages[12]='libpng-dev'

# ----- function to check if package exists -----
function package_exists(){
	dpkg -s "$1" &> /dev/null
	return $?
}

# ----- check if package exists for each element in array -----
for pkg in "${packages[@]}"
do
	if ! package_exists $pkg ; then
		echo "$pkg is not isntalled, installing now" 
		apt-get install -y $pkg
	else
		echo "$pkg is installed ... good to go"
	fi
done

# -----  build directory for cmake  -----
[ ! -d /opt/cmake ]
  mkdir /opt/cmake
cd /opt/cmake || exit

wget https://github.com/Kitware/CMake/releases/download/v3.19.6/cmake-3.19.6-Linux-x86_64.sh
sh ./cmake-3.19.6-Linux-x86_64.sh --prefix=/opt/cmake --skip-license
ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake

# -----  set directory for biotools and start cloning -----
[ ! -d /opt/biotools ]
	mkdir /opt/biotools
cd /opt/biotools || exit


git clone https://github.com/lh3/bwa.git
cd bwa || exit
make
cd ../

wget https://cab.spbu.ru/files/release3.15.1/SPAdes-3.15.1.tar.gz
tar -xzf SPAdes-3.15.1.tar.gz
cd SPAdes-3.15.1 || exit
PREFIX=/usr/local ./spades_compile.sh
cd /opt/biotools || exit

cd /usr/bin || exit
ln -s /opt/biotools/bwa/bwa
