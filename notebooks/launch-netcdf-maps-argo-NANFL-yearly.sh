for year in $(seq 2009 2013); do
		python 2020-11-02-AA-netcdf-maps-all-profiles-ARGO-NANFL-date-cal1.py --datemin '01-01-'${year} --datemax '31-12-'${year}
		python 2020-11-03-AA-maps-all-profiles-ARGO-NANFL-date-cal1.py --datemin '01-01-'${year} --datemax '31-12-'${year}

done

