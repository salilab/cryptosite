include ../../../Makefile.include

.PHONY: install
PY=${PYTHONDIR}/cryptosite/config

CONFIG=${PY}/__init__.py

install: ${CONFIG}

${CONFIG}::
	@if [ ! -d ${PY} ]; then mkdir -p ${PY}; fi
	echo "datadir=r'${DATADIR}'" > $@
	python -m py_compile $@
