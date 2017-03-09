include ../Makefile.include

.PHONY: install

FILES=${DATADIR}/Scaler_Final_Final.pkl ${DATADIR}/SVM_Final_Final.pkl \
      ${DATADIR}/PolyScaler_Final.pkl ${DATADIR}/PolySVC_FinalModel.pkl \
      ${DATADIR}/script.chimerax

install: ${FILES}

${DATADIR}/%: %
	@if [ ! -d ${DATADIR} ]; then mkdir -p ${DATADIR}; fi
	install $< $@
