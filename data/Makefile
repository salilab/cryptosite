include ../Makefile.include

.PHONY: install

FILES=${DATADIR}/Scaler_Final_Final.pkl ${DATADIR}/SVM_Final_Final.pkl

install: ${FILES}

${DATADIR}/%: %
	@if [ ! -d ${DATADIR} ]; then mkdir -p ${DATADIR}; fi
	install $< $@
