all: ngspice.log 
	$(eval NGSFIGS=$(shell grep _FIG ngspice.log | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	$(eval NGSTABS=$(shell grep _TAB ngspice.log | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $< | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > op_tab.tex;)
	

ngspice.log: ../sim/t3_fullwave.net
	ngspice -b $< -o $@

clean:
	@rm -f *.log *.tex *.pdf *.ps *~ *.txt *.eps

.PHONY: all clean
