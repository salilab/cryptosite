include Makefile.include

.PHONY: install test

install:
	${MAKE} -C data install
	${MAKE} -C data/ligands install
	${MAKE} -C bin install
	${MAKE} -C lib/cryptosite install
	${MAKE} -C lib/cryptosite/config install

test:
	cd test && nosetests
