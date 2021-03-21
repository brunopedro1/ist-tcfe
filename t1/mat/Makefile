all: octave.log
	$(eval OCTTABS=$(shell grep _TAB $< | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(OCTTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $< | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/@//g' | sed 's/$$/\\\\ \\hline/g' > $i_TAB.tex;)
	$(eval OCTTAB1S=$(shell grep _TAB1 $< | sed 's/_TAB1//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(OCTTAB1S), sed -n '/^$i_TAB1/,/^$i_END/{p;/^$i_END/q}' $< | grep -v $i_TAB1 | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/@//g' | sed 's/$$/\\\\ \\hline/g' > $i_TAB1.tex;)
	$(eval OCTTAB2S=$(shell grep _TAB2 $< | sed 's/_TAB2//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(OCTTAB2S), sed -n '/^$i_TAB2/,/^$i_END/{p;/^$i_END/q}' $< | grep -v $i_TAB2 | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/@//g' | sed 's/$$/\\\\ \\hline/g' > $i_TAB2.tex;)

octave.log: t1.m
	octave $< > $@

clean:
	@rm -f octave.log octave-workspace *.tex *~

.PHONY: clean
