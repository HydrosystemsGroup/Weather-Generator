.PHONY: mon db rda clean all

fetch:
	mkdir -p data
	curl -o data/east_monthly_met.tgz http://hydro.engr.scu.edu/files/gridded_obs/monthly/ascii/east_monthly_met.tgz
	tar -zxvf data/east_monthly_met.tgz -C data/
	rm data/east_monthly_met.tgz

db:
	mkdir -p db
	RScript --slave r/make_db.R $(abspath ./data/east) $(abspath ./db) maurer_mon.db

rda:
	RScript --slave r/make_rda.R $(abspath ./db/maurer_mon.db) $(abspath ../../data/maurer.rda)

clean:
	rm -r ./data
	rm -r ./db

all: fetch db rda