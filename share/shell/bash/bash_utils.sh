#!/bin/bash

# functions 
function check_err {
: '
  1) abstract: check the shell error code.
  2) history log:
    2017-08-09  junghan kim    initial setup
  3) usage:
    check_err $?
    check_err $? "some message"
'
    if [ $# -lt 1 ]; then
        echo "need error code... ex) check_err \$? message "
        exit 
    fi
    err=$1
  
    if [ $# -ge 2 ]; then
        message=$2
    else
        message='message none'
    fi
  
    if [ $err -ne 0 ]; then
        echo "error: " $message
        exit 
    fi
}
