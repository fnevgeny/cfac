DEMOS = structure pirr aidr ionization excitation \
	structureUTA pirrUTA aidrUTA ionizationUTA excitationUTA \
	aidrM excitationM

all: 

check: 
	@set -e; for i in $(DEMOS); do (cd $$i; $(MAKE) check) || exit 1; done

clean: 
	@set -e; for i in $(DEMOS); do (cd $$i; $(MAKE) clean) || exit 1; done

install : 
