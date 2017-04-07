PROG= pmfiberbiref

all:
	echo "Make clean, test, dist, windist, realclean, distclean, update-version"

dist:
	make distclean
	pyinstaller $(PROG).py
	cd dist; zip -9rv $(PROG) $(PROG)

update-version:
	cp version.py /home/dba/mesdb/coherent_data/utils/$(PROG)_version

windist:
	make distclean
	C:/python27/scripts/pyinstaller.exe $(PROG).py
	cd dist; zip -9rv $(PROG) $(PROG)

test:
	./$(PROG).py &

clean:
	rm -f Makefile.bak *.png *.dat *~ *.txt *.csv

realclean: clean distclean
	rm -f *.pyc

distclean:
	rm -rf dist build $(PROG).spec
