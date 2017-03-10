#! /bin/sh

if [ $# -lt 1  ]; then
    echo "Usage: $0 <program_executable>"
    exit 1
fi

program=$1
test_files=`ls -1 test/*.data`

for test in ${test_files}
do
    >&2 echo "./$program $test"
    ./$program $test
done
