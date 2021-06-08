all: octave.log
	$(eval OCTTABS=$(shell grep _TAB $< | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(OCTTABS), sed -n '/^$i_TAB/,/^$i_END_nodal1/{p;/^$i_END/q}' $< | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/@//g' | sed 's/$$/\\\\ \\hline/g' > $i_TAB.tex;)
	

octave.log: t5.m
	octave $< > $@
	
clean:
	@rm -f octave.log octave-workspace *.eps *~ *.png *_TAB.tex *.txt *.odg

.PHONY: clean


