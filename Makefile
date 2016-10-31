include Makefile.include

.PHONY: install

install:
	${MAKE} -C data install
	${MAKE} -C bin install
	${MAKE} -C lib/cryptosite install
	${MAKE} -C lib/cryptosite/config install
