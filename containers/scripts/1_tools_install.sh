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
packages[10]='cmake'
packages[11]='pkg-config'
packages[12]='libfreetype6-dev'
packages[13]='libpng-dev'

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

# -----  set directory for biotools and start cloning -----
[ ! -d /opt/biotools ]
	mkdir /opt/biotools
cd /opt/biotools || exit


git clone https://github.com/lh3/bwa.git
cd bwa || exit
make
cd ../

wget http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0.tar.gz
tar -xzf SPAdes-3.13.0.tar.gz
cd SPAdes-3.13.0 || exit
PREFIX=/usr/local ./spades_compile.sh
cd /opt/biotools || exit

cd /usr/bin || exit
ln -s /opt/biotools/bwa/bwa
