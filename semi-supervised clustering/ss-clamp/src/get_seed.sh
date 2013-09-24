#!/bin/bash
#
# Parameters:
#   $1: dabatabase name
#   $2: main name of the pwc file
#
database="$1.db"
pwc="%$2.pwc"
sqlite3 $database <<!
select a.seed
from algorithm a, input i
where a.sha1 = i.algorithm_id
and i.file like '$pwc';
!
