# Get File lists
# $1 : dataset name
# $2 : process name

das_client.py --query="file dataset="$1" instance=prod/phys03" --limit 1000 > list_$2.txt

grep DQMIO list_$2.txt > list_invalid_$2.txt

python $CRABPYTHON/DBS3SetFileStatus.py --url=https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter --status=invalid --recursive=False  --files=list_invalid_$2.txt
