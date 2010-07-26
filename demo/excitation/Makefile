SFAC = ../../sfac/sfac
all: fe17_excitation.py fe17_excitationM.py
	python fe17_excitation.py
	python fe17_excitationM.py

test: fe17_excitation.sf fe17_excitationM.sf
	$(SFAC) fe17_excitation.sf
	$(SFAC) fe17_excitationM.sf

clean:
	rm -f ne.*
