include ../../Makefile.include

.PHONY: install
PY=${PYTHONDIR}/cryptosite

FILES=${PY}/__init__.py ${PY}/cleaning.py \
      ${PY}/bmi_feature_parser.py ${PY}/res_parser_bmi.py \
      ${PY}/hyd_chr_sse.py ${PY}/patch_mapper.py ${PY}/seq_conservation.py \
      ${PY}/predict.py ${PY}/pockets.py ${PY}/am_bmi.py \
      ${PY}/soap.py ${PY}/setup.py \
      ${PY}/gather.py ${PY}/analysis.py ${PY}/chimera.py

install: ${FILES}

${PY}/%: %
	@if [ ! -d ${PY} ]; then mkdir -p ${PY}; fi
	install $< $@
	${PYTHON} -m py_compile $@
