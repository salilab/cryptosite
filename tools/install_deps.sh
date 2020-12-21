#!/bin/bash -e

# Install all dependencies under the provided directory, if they're not
# there already. This is mainly intended for use by Travis CI
# (see ../.travis.yml)

if [ $# -ne 2 ]; then
  echo "Usage: $0 top_directory python_version"
  exit 1
fi

top_dir=$1
python_version=$2
bin_dir=${top_dir}/bin
lib_dir=${top_dir}/lib
temp_dir=`mktemp -d`

cd ${temp_dir}

mkdir -p ${bin_dir} ${lib_dir}

conda update --yes -q conda
conda create --yes -q -n python${python_version} -c salilab python=${python_version} pip biopython scikit-learn scipy modeller
eval "$(conda shell.bash hook)"
conda activate python${python_version}
pip install coverage pytest-cov pytest-flake8

# MUSCLE
if [ ! -e ${bin_dir}/muscle ]; then
  wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
  tar -xvzf muscle3.8.31_i86linux64.tar.gz
  mv muscle3.8.31_i86linux64 ${bin_dir}/muscle
fi

# DSSP
if [ ! -e ${bin_dir}/mkdssp ]; then
  # This is a Sali-lab-maintained mirror, since the "official" DSSP FTP
  # server is often down (which will cause our Travis builds to fail)
  wget https://salilab.org/dssp/dssp-2.0.4-linux-amd64
# wget ftp://ftp.cmbi.ru.nl/pub/software/dssp/dssp-2.0.4-linux-amd64
  chmod a+x dssp-2.0.4-linux-amd64
  mv dssp-2.0.4-linux-amd64 ${bin_dir}/mkdssp
fi

# fpocket2
if [ ! -e ${bin_dir}/fpocket ]; then
  wget http://downloads.sourceforge.net/project/fpocket/fpocket2.tar.gz
  tar -xzf fpocket2.tar.gz
  (cd fpocket2 && sed -e 's/\$(LFLAGS) \$^/\$^ \$(LFLAGS)/' makefile > makefile.new && mv makefile.new makefile && make bin/fpocket && cp bin/fpocket ${bin_dir})
fi

# PatchDock
if [ ! -e ${bin_dir}/patch_dock.Linux ]; then
  wget http://bioinfo3d.cs.tau.ac.il/PatchDock/download/patch_dock_download.zip
  unzip patch_dock_download.zip
  (cd PatchDock && cp -r * ${bin_dir})
fi

cd /
rm -rf ${temp_dir}
