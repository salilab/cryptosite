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
	pytest -v .
	flake8 --ignore=E402,W503,W504 --exclude=doc/conf.py
