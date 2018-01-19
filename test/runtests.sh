#!/bin/bash
LOGFILE="$(pwd)/test.log"

failed=0

print() {
  echo $@ | tee -a $LOGFILE
}
buildtest () {
  test_name=$(basename $1)
  print
  print -------------------------
  print Building Test: $test_name
  print -------------------------
  
  make build 2>&1 | sed -e 's/^/'"$test_name"': /' | tee -a $LOGFILE
  rc=${PIPESTATUS[0]}
  
  #    print " Returned: $rc " 
  
  if  [ $rc -eq 0 ];
  then
    print Building TEST $1 PASSED !
  else
    print Building TEST $1 FAILED !
    failed=$[$failed + 1]
  fi
}

runtest () {
  test_name=$(basename $1)
  print
  print ------------------------
  print Running Test: $test_name
  print -------------------------
  
  make run 2>&1 | sed -e 's/^/'"$test_name"': /' | tee -a $LOGFILE
  rc=${PIPESTATUS[0]}
  
  #    print " Returned: $rc " 
  
  if  [ $rc -eq 0 ];
  then
    print Test $1 PASSED !
  else
    print Test $1 FAILED !
    failed=$[$failed + 1]
  fi
}
initlog()
{
  echo `date`: Starting tests > $LOGFILE
}
initlog

# Build all tests

for test in `ls -d */`
do
  if test -x $test
  then 
    cd $test
    buildtest $test
    cd ..
  fi
done

# Run all tests

for test in `ls -d */`
do
  if test -x $test
  then 
    cd $test
    runtest $test
    cd ..
  fi
done

# Cleanup all tests

for test in `ls -d */`
do
  if test -x $test
  then 
    cd $test
    make clean
    cd ..
  fi
done
if [ $failed -ne 0 ]
then
  echo
  echo $failed test\(s\) FAILED see $LOGFILE for more details.
fi
