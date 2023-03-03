#!/bin/bash

for filename in /Users/yashravipati/Downloads/DockedFiles/*.pdbqt; do
    # get the name of the file without the extension
    base=`basename $filename .pdbqt`
    # convert the pdbqt file to pdb format
    obabel -ipdbqt $filename -opdb -O /Users/yashravipati/Downloads/DockedFiles/$base.pdb
done