include ../Makefile.include

.PHONY: install

FILES=${BINDIR}/cryptosite

install: ${FILES}

${BINDIR}/%: %
	@if [ ! -d ${BINDIR} ]; then mkdir -p ${BINDIR}; fi
	install $< $@
