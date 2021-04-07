all: octave.log
	$(eval OCTTABS=$(shell grep _TAB_nodal1 $< | sed 's/_TAB_nodal1//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(OCTTABS), sed -n '/^$i_TAB_nodal1/,/^$i_END_nodal1/{p;/^$i_END_nodal1/q}' $< | grep -v $i_TAB_nodal1 | grep -v $i_END_nodal1 | grep -v '#' | sed 's/\=/\&/g' | sed 's/@//g' | sed 's/$$/\\\\ \\hline/g' > $i_TAB_nodal1.tex;)
	$(eval OCTTABS=$(shell grep _TAB_nodal2 $< | sed 's/_TAB_nodal2//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(OCTTABS), sed -n '/^$i_TAB_nodal2/,/^$i_END_nodal2/{p;/^$i_END_nodal2/q}' $< | grep -v $i_TAB_nodal2 | grep -v $i_END_nodal2 | grep -v '#' | sed 's/\=/\&/g' | sed 's/@//g' | sed 's/$$/\\\\ \\hline/g' > $i_TAB_nodal2.tex;)
	$(eval OCTTABS=$(shell grep _TAB_nodal4 $< | sed 's/_TAB_nodal4//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(OCTTABS), sed -n '/^$i_TAB_nodal4/,/^$i_END_nodal4/{p;/^$i_END_nodal4/q}' $< | grep -v $i_TAB_nodal4 | grep -v $i_END_nodal4 | grep -v '#' | sed 's/\=/\&/g' | sed 's/@//g' | sed 's/$$/\\\\ \\hline/g' > $i_TAB_nodal4.tex;)
	$(eval OCTTABS=$(shell grep _TAB_nodal4ph $< | sed 's/_TAB_nodal4ph//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(OCTTABS), sed -n '/^$i_TAB_nodal4ph/,/^$i_END_nodal4ph/{p;/^$i_END_nodal4ph/q}' $< | grep -v $i_TAB_nodal4ph | grep -v $i_END_nodal4 | grep -v '#' | sed 's/\=/\&/g' | sed 's/@//g' | sed 's/$$/\\\\ \\hline/g' > $i_TAB_nodal4ph.tex;)


octave.log: t2.m
	octave $< > $@
	
clean:
	@rm -f octave.log octave-workspace *.eps *~ *.png *_TAB_nodal1.tex *_TAB_nodal2.tex *_TAB_nodal4.tex *_TAB_nodal4ph.tex

.PHONY: clean


