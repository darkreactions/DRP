featureset=`cat $1 | egrep -o "Selected attributes: [[:digit:]]+(,[[:digit:]]+)*" | egrep -o "[[:digit:]]+(,[[:digit:]]+)*"`
featureset=$featureset",last"
echo $featureset