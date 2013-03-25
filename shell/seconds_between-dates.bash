#!/bin/bash
if test $1 = '-h'
then
echo "*****************************"
echo "*** Difference in seconds ***"
echo "***   between two dates   ***"
echo "***     (up to days)      ***"
echo "***     DATE2 - DATE1     ***"
echo "*****************************"
echo "seconds_between-dates.bash 'DATE1'(date 1, [YYYY][MM][DD][HH][MI][SS]) \
'DATE2' (data 2, [YYYY][MM][DD][HH][MI][SS])"
else
file=$0
llfile=`(expr length $file)`
path=`(expr $llfile - 26)`
foresthome=`(expr substr $file 1 $path)`
foresthome=${foresthome}/..

# Date1
##
YMD1=`(expr substr $1 1 8)`
hours1=`(expr substr $1 9 2)`
minutes1=`(expr substr $1 11 2)`
seconds1=`(expr substr $1 13 2)`
hours1=`(expr $hours1 '*' 3600)`
minutes1=`(expr $minutes1 '*' 60)`
seconds1=`(expr $hours1 + $minutes1 + $seconds1)`
sec1=`(date +%s -u -d"$YMD1 $seconds1 seconds")`

# Date2
##
YMD2=`(expr substr $2 1 8)`
hours2=`(expr substr $2 9 2)`
minutes2=`(expr substr $2 11 2)`
seconds2=`(expr substr $2 13 2)`
hours2=`(expr $hours2 '*' 3600)`
minutes2=`(expr $minutes2 '*' 60)`
seconds2=`(expr $hours2 + $minutes2 + $seconds2)`
sec2=`(date +%s -u -d"$YMD2 $seconds2 seconds")`
seconds=`(expr $sec2 - $sec1)`

echo $seconds
fi
