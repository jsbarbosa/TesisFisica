thesis.pdf : thesis.tex Chapters/*.tex
	pdflatex $<

thesis_lowres.pdf : thesis.pdf
	gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/ebook -dNOPAUSE -dBATCH -sOutputFile=$@ $<

clean : 
	rm -r *.aux *.out *.toc *.gz *.xml *.lof *.loc *.log *.lot *.blg *.bbl Chapters/*.aux

