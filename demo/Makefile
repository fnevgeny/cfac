DEMOS = structure excitation excitationM pirr aidr ionization

all: 

check: 
	@set -e; for i in $(DEMOS); do (cd $$i; $(MAKE) check) || exit 1; done

clean: 
	@set -e; for i in $(DEMOS); do (cd $$i; $(MAKE) clean) || exit 1; done

install : 
