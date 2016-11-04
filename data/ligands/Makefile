include ../../Makefile.include

.PHONY: install

D=${DATADIR}/ligands

FILES=${D}/2F2.mol2 ${D}/ACM.mol2 ${D}/BNZ.mol2 ${D}/CHX.mol2 \
      ${D}/EHN.mol2 ${D}/HBX.mol2 ${D}/IPH.mol2 ${D}/NME.mol2 \
      ${D}/2F2.pdb ${D}/ACM.pdb ${D}/BNZ.pdb ${D}/CHX.pdb \
      ${D}/EHN.pdb ${D}/HBX.pdb ${D}/IPH.pdb ${D}/NME.pdb \
      ${D}/ACE.mol2 ${D}/ACN.mol2 ${D}/CCN.mol2 ${D}/DMF.mol2 \
      ${D}/EOH.mol2 ${D}/IPA.mol2 ${D}/ISB.mol2 ${D}/URE.mol2 \
      ${D}/ACE.pdb ${D}/ACN.pdb ${D}/CCN.pdb ${D}/DMF.pdb \
      ${D}/EOH.pdb ${D}/IPA.pdb ${D}/ISB.pdb ${D}/URE.pdb

install: ${FILES}

${D}/%: %
	@if [ ! -d ${D} ]; then mkdir -p ${D}; fi
	install $< $@
