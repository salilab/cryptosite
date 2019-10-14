include Makefile.include

.PHONY: install test bin clean

bin:
	${MAKE} -C bin

clean:
	${MAKE} -C bin clean

install: bin
	${MAKE} -C data install
	${MAKE} -C data/ligands install
	${MAKE} -C bin install
	${MAKE} -C lib/cryptosite install
	${MAKE} -C lib/cryptosite/config install

test: bin
	nosetests --processes=-1 test
