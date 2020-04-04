#!/bin/bash

for file in $(ls *.gjf)
    do
    prefix=`basename $file .gjf`
    g09 $file
    mv fort.7 $prefix.mo
    done
