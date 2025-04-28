#!/bin/bash

: '
Preparation for Intel MKL compilation
argument will be stored as compiler option
If no argument is passed, `SDL(Single Dynamic Linking)` will be applied.
More information of default compiler option is as follows:
'

mklflags="$@"
prefix="MKLFLAGS="
repeat=false
timeout=15
file="Makefile.inc"

while IFS= read -r line
do
    sline=$(sed "s/ //g" <<< $line)
    if [[ "$sline" == *"$prefix"* ]]; then
        # echo "Found"
        repeat=true
        break
    fi
done < "$file"

if [[ $repeat == true ]]; then
    echo "MKLFLAGS already exists in Makefile.inc file."
    read -t "$timeout" -p "Do you want to overwrite MKLFLAGS? [yes|no]: " user_input

    if [[ $? -ne 0 ]]; then
        echo
        echo "Timeout. No change applied."
        exit 0
    else
        user_input=$(echo "$user_input" | tr '[:upper:]' '[:lower:]')

        if [[ "$user_input" == "yes" || "$user_input" == "y" ]]; then
            sed -i -E \
            "s|^[[:space:]]*MKLFLAGS[[:space:]]*=[[:space:]]*.*|MKLFLAGS = $mklflags|" $file
            echo "Flags overwritten"
        fi
    fi
else
    # Write MKLFLAGS into `Makefile.inc` file
    sed -i "\$a MKLFLAGS = $mklflags" $file
fi
