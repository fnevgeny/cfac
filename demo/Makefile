DEMOS = structure pirr aidr ionization excitation excitationM \
	structureUTA pirrUTA aidrUTA ionizationUTA excitationUTA 

all: 

check: 
	@set -e; for i in $(DEMOS); do (cd $$i; $(MAKE) check) || exit 1; done

clean: 
	@set -e; for i in $(DEMOS); do (cd $$i; $(MAKE) clean) || exit 1; done

install : 
