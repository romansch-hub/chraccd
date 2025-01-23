timestamp(){
	date +"%T"
}

timestamp

while IFS= read -r line
	do 
		awk '{if ($11 == "nan") $11 = 0;
        		if ($12 == "nan") $12 = 0;
        		if ($13 == "nan") $13 = 0;
        		print $1"\t"$2"\t"$3"\t"$11"\t"$13"\t"$12}' > out.bed
	done < mod-counts.bed


timestamp
