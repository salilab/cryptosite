include ../Makefile.include

.PHONY: install clean

cryptosite: cryptosite.in
	FULL_PYTHON=`which ${PYTHON}`; sed -e "s^#!.*^#!$$FULL_PYTHON^" $< > $@ && chmod a+x $@

FILES=${BINDIR}/cryptosite

install: ${FILES}

clean:
	rm -f cryptosite

${BINDIR}/%: %
	@if [ ! -d ${BINDIR} ]; then mkdir -p ${BINDIR}; fi
	install $< $@
