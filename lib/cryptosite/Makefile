include ../../Makefile.include

.PHONY: install
PY=${PYTHONDIR}/cryptosite

FILES=${PY}/__init__.py ${PY}/chasa.py ${PY}/cleaning.py \
      ${PY}/bmi_feature_parser.py ${PY}/res_parser_bmi.py \
      ${PY}/hyd_chr_sse.py ${PY}/patch_mapper.py ${PY}/seq_conservation.py \
      ${PY}/predicter.py ${PY}/pocket_parser.py ${PY}/am_bmi.py \
      ${PY}/soap_clean.py

install: ${FILES}

${PY}/%: %
	@if [ ! -d ${PY} ]; then mkdir -p ${PY}; fi
	install $< $@
	python -m py_compile $@
