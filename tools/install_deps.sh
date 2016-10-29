#!/bin/bash -e

# Install all dependencies under the provided directory, if they're not
# there already. This is mainly intended for use by Travis CI
# (see ../.travis.yml)

if [ $# -ne 3 ]; then
  echo "Usage: $0 top_directory python_version modeller_license_file"
  exit 1
fi

top_dir=$1
python_version=$2
modeller_license_file=$3
bin_dir=${top_dir}/bin
lib_dir=${top_dir}/lib
conda_dir=${top_dir}/miniconda
temp_dir=`mktemp -d`

cd ${temp_dir}

mkdir -p ${bin_dir} ${lib_dir}

# Use miniconda Python rather than the Travis environment (we do this because
# the latter has terrible support for scipy, while the former trivially supports
# both it and Sali lab packages like Modeller)
if [ ! -e ${conda_dir}/envs/python${python_version} ]; then
  # Clean up after a potential previous install failure
  rm -rf ${conda_dir}
  # Save on some downloading if the version is the same
  if [ "${python_version}" == "2.7" ]; then
    wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh
  else
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
  fi
  bash miniconda.sh -b -p ${conda_dir}
  export PATH=${conda_dir}/bin:$PATH
  conda update --yes -q conda
  conda create --yes -q -n python${python_version} python=${python_version}
fi

# Make sure that our conda environment is up to date
export PATH=${conda_dir}/bin:$PATH
conda update --yes -q conda
source activate python${python_version}

export KEY_MODELLER=`cat ${modeller_license_file}`
conda install --yes -c salilab pip biopython scikit-learn scipy modeller nose
pip install coverage

# MUSCLE
if [ ! -e ${bin_dir}/muscle ]; then
  wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
  tar -xvzf muscle3.8.31_i86linux64.tar.gz
  mv muscle3.8.31_i86linux64 ${bin_dir}/muscle
fi

# DSSP
if [ ! -e ${bin_dir}/mkdssp ]; then
  wget ftp://ftp.cmbi.ru.nl/pub/software/dssp/dssp-2.0.4-linux-amd64
  chmod a+x dssp-2.0.4-linux-amd64
  mv dssp-2.0.4-linux-amd64 ${bin_dir}/mkdssp
fi

# fpocket2
if [ ! -e ${bin_dir}/fpocket ]; then
  wget http://downloads.sourceforge.net/project/fpocket/fpocket2.tar.gz
  tar -xzf fpocket2.tar.gz
  (cd fpocket2 && sed -e 's/\$(LFLAGS) \$^/\$^ \$(LFLAGS)/' makefile > makefile.new && mv makefile.new makefile && make bin/fpocket && cp bin/fpocket ${bin_dir})
fi

cd /
rm -rf ${temp_dir}
