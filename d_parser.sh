if [ $# -eq 0 ]; then
  echo -n ""
else
  list=( "$@" )
  IFS=:
  deplist="${list[*]}"
  echo -n "-W depend=afterok:$deplist"
fi;
