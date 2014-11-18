XPI = cfacdb.xpi

INSTALL_RDF = install.rdf
MANIFEST    = chrome.manifest

CONTENT_DIR = chrome/content
SKIN_DIR    = chrome/skin/classic
ICON_DIR    = chrome/icons/default

CHROME_SRC = \
	$(CONTENT_DIR)/overlay.xul \
	$(CONTENT_DIR)/overlay.js \
	$(CONTENT_DIR)/cfacdb.xul \
	$(CONTENT_DIR)/cfacdb.js

all : $(XPI)

$(XPI) : $(INSTALL_RDF) $(MANIFEST) $(CHROME_SRC)
	rm -f $@ && zip -r $@ $^

release : $(XPI)
	./scripts/makerelease.sh

upload :
	rsync -av releases/*.xpi releases/update.rdf plasma-gate.weizmann.ac.il:cfacdb/xul/releases

clean :
	rm -f $(XPI) $(CASES_RDF)

changes :
	cvs2cl
