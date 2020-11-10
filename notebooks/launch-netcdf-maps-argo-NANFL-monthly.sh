for year in $(seq 2009 2013); do
	for month in $(seq 1 12); do

		mm=$(printf "%02d" $month)
		case $month in
			1|3|5|7|8|10|12) end=31 ;;
			4|6|9|11) end=30;;
			2) end=28;;
		esac
		python 2020-11-02-AA-netcdf-maps-all-profiles-ARGO-NANFL-date-cal1.py --datemin '01-'${mm}'-'${year} --datemax ${end}'-'${mm}'-'${year}
		python 2020-11-03-AA-maps-all-profiles-ARGO-NANFL-date-cal1.py --datemin '01-'${mm}'-'${year} --datemax ${end}'-'${mm}'-'${year}

	done
done

