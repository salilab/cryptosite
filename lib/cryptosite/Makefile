include ../../Makefile.include

.PHONY: install
PY=${PYTHONDIR}/cryptosite

FILES=${PY}/__init__.py

install: ${FILES}

${PY}/%: %
	@if [ ! -d ${PY} ]; then mkdir -p ${PY}; fi
	install $< $@
	python -m py_compile $@
