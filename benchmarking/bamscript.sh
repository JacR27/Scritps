#
# A script for finding the time for each stage of BAM creation in Isaac and outputing it as a CSV file for plotting
#
# Paul Smith 26 August 2014
#
# 
# Find the first entry and then the entries that correspond to each bin
grep -P "Making sure all bins fit in memory done|\.dat" | \
#
# Remove all the unwatend lines
grep -v Break | grep -v offset | grep -v isaac-align |grep -v WAR| grep -v gnuplot | \
# 
# Parse the strings for the bin name so they are all equivalent
sed 's/)//g' | sed 's/"//g'|\
#
# Get the time stamp
awk '{ cmd="date \"+%s\" -d \""$1" "$2"\""; cmd | getline time;
if (NR==1) {start=time; print "a-a- 000000000bin,load wait,load, re-align wait,re-align, serialise, save wait,save" }}
/Saving/{if (NF==12) startsave[$12$3]=time-start; else {endsave[$12$3]=time-start;}}
/Reading align/{if (NF==16) {startread[$16$3]=time-start;} else { endread[$17$3]=time-start; }} 
/Reading unalign/{if (NF==16) {startread[$16$3]=time-start;} else { endread[$17$3]=time-start;startrel[$17$3]=time-start }} 
/Serializing records/{if (NF==19) startser[$19$3]=time-start; else {endser[$20$3]=time-start;}} 
/Realigning against/{ startrel[$17$3]=time-start} 
END{for (var in startsave) print var,",",startread[var],",",endread[var]-startread[var],",",startrel[var]-endread[var],",",startser[var]-startrel[var],",",endser[var]-startser[var],",",startsave[var]-endser[var],",",endsave[var]-startsave[var]}' | \
cut -f 3- -d- | \
sed 's/00000000.dat/100000000.dat/' | \
sort | \
sed 's/100000000.dat/00000000.dat/' | \
sed 's/000000000bin/bin/'
