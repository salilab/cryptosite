include ../../../Makefile.include

.PHONY: install
PY=${PYTHONDIR}/cryptosite/config

CONFIG=${PY}/__init__.py

install: ${CONFIG}

${CONFIG}::
	@if [ ! -d ${PY} ]; then mkdir -p ${PY}; fi
	@if [ ! -e ${UNIPROT} ]; then echo "Cannot find UniProt database as ${UNIPROT}; set the make UNIPROT variable to its location"; exit 1; fi
	echo "datadir=r'${DATADIR}'" > $@
	echo "uniprot=r'${UNIPROT}'" >> $@
	${PYTHON} -m py_compile $@
