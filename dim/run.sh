#!/bin/bash

for file in $(ls *.gjf)
    do
    prefix=`basename $file .gjf`
    g09 $file
    mv FILE.43 $prefix.om
    mv FILE.44 $prefix.fm
    done
