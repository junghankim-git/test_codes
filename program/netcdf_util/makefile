
build:
	make -C src $@

gen:
	make -C src $@

clean: clean_test
	make -C src $@

nc4:
	nc_util -m 1 input/pnc.nc -o nc4.nc

pnc:
	nc_util -m 2 input/nc4.nc -o pnc.nc

dnc4:
	ncdump nc4.nc > new.txt
	diff input/org.txt new.txt 2>&1

dpnc:
	ncmpidump pnc.nc > new.txt
	diff input/org_pnc.txt new.txt 2>&1

clean_test:
	rm -f *.txt *.nc
