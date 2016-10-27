#!/bin/sh -e

# Install all dependencies under the provided directory, if they're not
# there already. This is mainly intended for use by Travis CI
# (see ../.travis.yml)

if [ $# -ne 2 ]; then
  echo "Usage: $0 top_directory modeller_license_file"
  exit 1
fi

top_dir=$1
modeller_license_file=$2
bin_dir=${top_dir}/bin
lib_dir=${top_dir}/lib
python_dir=${top_dir}/python
temp_dir=`mktemp -d`

cd ${temp_dir}

mkdir -p ${bin_dir} ${lib_dir} ${python_dir}

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

# MODELLER
if [ ! -e ${python_dir}/_modeller.so ]; then
  MODVER=9.17
  mod_dir=${top_dir}/modeller-${MODVER}
  wget https://salilab.org/modeller/${MODVER}/modeller-${MODVER}.tar.gz
  tar -xzf modeller-${MODVER}.tar.gz
  echo -e "\n${mod_dir}" > inst-pre
  echo -e "\n\n\n" > inst-post
  cat inst-pre ${modeller_license_file} inst-post > inst-cmd
  (cd modeller-${MODVER} && ./Install < ../inst-cmd > /dev/null)
  ln -sf ${mod_dir}/bin/mod${MODVER} ${bin_dir}
  ln -sf ${mod_dir}/lib/x86_64-intel8/lib*so* ${lib_dir}
  rm -f ${python_dir}/modeller
  ln -sf ${mod_dir}/modlib/modeller ${python_dir}
  ln -sf ${mod_dir}/lib/x86_64-intel8/python2.5/_modeller.so ${python_dir}
fi

cd /
rm -rf ${temp_dir}
