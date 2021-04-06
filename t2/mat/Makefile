all: octave.log
	$(eval OCTTABS=$(shell grep _TAB_nodal1 $< | sed 's/_TAB_nodal1//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(OCTTABS), sed -n '/^$i_TAB_nodal1/,/^$i_END_nodal1/{p;/^$i_END_nodal1/q}' $< | grep -v $i_TAB_nodal1 | grep -v $i_END_nodal1 | grep -v '#' | sed 's/\=/\&/g' | sed 's/@//g' | sed 's/$$/\\\\ \\hline/g' > $i_TAB_nodal1.tex;)

octave.log: t2.m
	octave $< > $@
	
clean:
	@rm -f octave.log octave-workspace *.eps *~ *.png *_TAB_nodal1.tex

.PHONY: clean
