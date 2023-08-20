#!/bin/bash

# define colors
YEL='\x1b[1;33m'
GRE='\x1b[1;32m'
RED='\x1b[1;31m'
BLU='\x1b[1;34m'
NC='\e[0m'

script_dir="$(dirname "$(readlink -f "$0")")"
script_dir="${script_dir%'/utils'}"

# error to terminate if failed
error_cheker(){
    lastexit=$1
    if [[ $(( lastexit )) -ne 0 ]];then
        echo -e "${RED}ERROR in exit value of last process${NC}"
        exit
    fi
}

destination=$1
echo "$destination"
if ( -z "$destination"  ) then 
    cmd=""
else
    cmd="--directory-prefix $destination"
    mkdir -p "$destination"
fi

# create env from yml file 
echo -e "${YEL}Creating wes_gatk ENV${NC}"
mamba env create -f ${script_dir}/workflow/env/wes_gatk.yml
error_cheker $?

# download kraken2 database
echo -e "${YEL}Downloading kraken2 database${NC}"
wget $cmd -c "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230605.tar.gz"
error_cheker $?
