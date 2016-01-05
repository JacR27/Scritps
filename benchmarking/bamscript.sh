#!/bin/bash

# Find the first entry and then the entries that correspond to each bin
grep -P "Making sure all bins fit in memory done|\.dat" | \
#
# Remove all the unwatend lines
grep -v Break | grep -v isaac-align |grep -v WAR| grep -v gnuplot | \
# 
# Parse the strings for the bin name so they are all equivalent
sed 's/)//g' | sed 's/"//g'|\
#
# Get the time stamp
gawk '

{cmd="date \"+%s\" -d \""$1" "$2"\""; cmd | getline time; close(cmd);
#
# Print column headings and make sure the titles come first in a sort 
if (NR==1) {print "a-a- 000000000bin,load wait,load, dup wait, dup, re-align wait, re-align, sortbam wait, sortbam,ser wait, serialise, save wait,save" }}

/Saving/{if (NF==12) startsave[$12]=time-start; else {endsave[$12]=time-start;}}

/Reading alignment records from/{if (length(start) == 0) {start = time;} startread[$16]=time-start;}

/Reading alignment records done/{endread[$17]=time-start;}

/Reading unaligned records from/{ if (length(startread[$16]) == 0) startread[$16]=time-start;}

/Reading unaligned records done/{
        endread[$17]=time-start;
        startrel[$17]=time-start;endrel[$17]=time-start;
        startdup[$17]=time-start;enddup[$17]=time-start;
}

/Resolving duplicates for/{ if (length(startdup[$16]) == 0) startdup[$16]=time-start;}

/Resolving duplicates done/{if (length(enddup[$17]) == 0) enddup[$17]=time-start;}

/Sorting offsets for bam BinMetadata/{startsortbam[$16]=time-start;}

/Sorting offsets for bam done/ {endsortbam[$17]=time-start;}

/Serializing records/{if (NF==19) startser[$19]=time-start; else {endser[$20]=time-start;}}

/Realigning against/{ startrel[$17]=time-start}

/Realigning gaps/{ endrel[$15]=time-start}


END{for (var in startsave) print var,",",startread[var],",",endread[var]-startread[var],",",startdup[var]-endread[var],",",enddup[var]-startdup[var],",",startrel[var]-enddup[var],",",endrel[var]-startrel[var],",",startsortbam[var]-endrel[var],",",endsortbam[var]-startsortbam[var],",",startser[var]-endsortbam[var],",",endser[var]-startser[var],",",startsave[var]-endser[var],",",endsave[var]-startsave[var]}' | \



# Cut down the name of the bin to just the last part of the name
cut -f 3- -d- | \
#
# Get rid of empty bins
#grep -v ",  ,"
#
# Change the name of the zero unaligned bin so it comes last in a sort.
sed 's/00000000.dat/100000000.dat/' | \
#
# Sort bins in numerical order
sort | \
#
# remove sort tag on zero bin so it now has the correct name
sed 's/100000000.dat/00000000.dat/' | \
#
# remove sort tag from column heading so it has correct name
sed 's/000000000bin/bin/' | \
# Get rid of empty bins
grep -v ",  ,"



