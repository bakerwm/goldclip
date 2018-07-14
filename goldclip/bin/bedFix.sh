#!/usr/bin/env bash
#
# Date: 2017-08-04
#
# parse: parameter or STDIN
#[ $# -ge 1 -a -f "$1" ] && input="$1" || input="/dev/stdin"
[[ $# -lt 1 ]] && echo "Usage: $(basename $0) <in.bed> <rename:0|1>" && exit 0
input=$1
flag_tag=$2
[[ ! -a ${input} ]] && echo "${input} - file not exists" && exit 1
[[ -z ${flag_tag} ]] && flag_tag=0 # default not change id
case ${flag_tag} in
    0|n|no)    flag_tag=0 ;;
    1|y|yes)   flag_tag=1 ;;
    *)         echo "${flag_tag} - unknown value" && exit 1 ;;
esac

## fix BED
# 0. chr (chr)
# 1. start <= end
# 2. id (width < 40)
# 3. score (integer)
# 4. strand (+, -, *)
## add tag to <name> field
num=1
#flag_tag=1
while IFS=$'\t' read -r -a va ; do
    [[ ${#va[@]} -lt 6 ]] && continue # skip 
    chr=${va[0]}
    start=${va[1]}
    end=${va[2]}
    id=${va[3]}
    score=${va[4]}
    strand=${va[5]}
    ## 0. chr
    [[ ${va[0]} != chr* ]] && chr="chr${va[0]}"
    ## 1. start <= end
    [[ ${va[1]} -gt ${va[2]} ]] && start=${va[2]} && end=${va[1]}
    [[ ${start} -lt 0 ]] && start=0
    [[ ${end} -lt 1 ]] && end=1
    ## 2. id
    id="${va[3]:0:40}"
    [[ ${flag_tag} -eq 1 ]] && id="${id}_${num}"
    ## 3. score
    #re='^[0-9]+$' # integer
    re='^[0-9]+([.][0-9]+)?$' # handle float
    #re='^-?[0-9]+([.][0-9]+)?$' # handle negative numbers
    [[ ${va[4]} =~ $re ]] || score=255 # default
    score=$(printf "%.0f" ${score})
    ## 4. strand
    st='^[+-]$'
    [[ ${va[5]} =~ $st ]] || strand="*" # default *
    # output
    #echo ${chr} ${start} ${end} ${id} ${score} "${strand}"
    printf "%s\t%d\t%d\t%s\t%d\t%s\n" ${chr} ${start} ${end} ${id} ${score} ${strand}
#    break
    num=$((num+1))
done < ${input}

## version-2017-01-18
##
## 1. fix the coordinates of BED file
## 2. chop the id_field of BED (limit: 20 characters)
## save original file to *.bak
#
#[[ $# -lt 1 ]] && echo "Usage: $0 <in.bed>" && exit 0
#in=$1
#bak="${in}.bak"
#cp -rf ${in} ${bak}
#cat ${bak} | awk '{FS=OFS="\t"}{$4=substr($4,0,20); if($2>$3){t=$2;$2=$3;$3=t} print $0}' > ${in}


## EOF
