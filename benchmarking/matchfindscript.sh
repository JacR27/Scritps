grep -P 'Selecting matches|Loading Bcl'| grep -v expectedBinSize | sed 's/,//g' | \
awk '{ cmd="date \"+%s\" -d \""$1" "$2"\""; cmd | getline time; if (NR==1) {start=time;print "0 , tile, load wait, load, select wait, select"}}\
/Selecting matches/{ if (NF==15) { startsel[$12,$11]=time-start } else {endsel[$16,$15]=time-start}}
/Loading Bcl/{if (NF==13) startread[$10,$9]=time-start; else { endread[$11,$10]=time-start}}
END{ for  (var in startread) { split(var,xx,SUBSEP); print xx[1],",",xx[2]," ",startread[var],"\011"endread[var]-startread[var]," ",startsel[var]-endread[var]," ",endsel[var]-startsel[var]}}' |\
sort -k1,1 -k3,3nr | sed 's/ , //' | sed 's/0tile/tile/'

