#!/usr/bin/env bash
#
# Date: 2017-08-05
# Author: WM
#
# Calling peaks using clipper (version 1.1)
# Refer to: https://github.com/YeoLab/clipper/wiki for instructions 
#

function echo2() {
    date_fmt='+%Y-%m-%d %H:%M:%S'
    case $2 in 
        error)     echo -e "[$(date "${date_fmt}")] error: $1" && exit 1 ;;
        warning)   echo -e "[$(date "${date_fmt}")] warning: $1" && exit 1 ;;
        *)         echo -e "[$(date "${date_fmt}")] $1" ;;
    esac
}

function get_python_ver() {
    # return the version of python
    # like: 2.7.12
    python -c 'import sys; version=sys.version_info[:3]; print("{0}.{1}.{2}".format(*version))'
}


function env_creater() {
    env_in=$1
    py_ver=$2 
    [[ -z $env_in ]] && echo2 "require env_path" "error"
    [[ -z $py_ver ]] && py_ver=2.7
    if [[ -d ${env_in} ]] ; then
        echo2 "env_path exists: ${env_in}"
    else
        virtualenv --python=python${py_ver} ${env_in}
    fi
    # env_checker ${env_in} ${py_ver}
}



function env_checker() {
    env_in=$1
    py_ver=$2
    create=$3 # 0=no, 1=yes
    [[ -z $env_in ]] && echo2 "require env_path" "error"
    [[ -z $py_ver ]] && py_ver=2.7
    [[ -z ${create} ]] && create=1 # yes
    env_name=$(basename ${env_in}) # name of env
    if [[ ${VIRTUAL_ENV} = *${env_name} ]] ; then
        py_cur=$(get_python_ver)
        a=${py_ver:0:1} # (2).7
        b=${py_cur:0:1} # (2).7
        if [[ $a = $b ]] ; then
            # echo "0" # success
            echo2 "env: ok, python: ok, ${env_in}"
        else
            echo2 "env: ok, python: failed, ${env_in}" "error"
            # echo "1" # env:exists, python:not_correct
        fi
    elif [[ -f ${env_in}/bin/activate ]] ; then
        source ${env_in}/bin/activate # switch to
        py_cur=$(get_python_ver)
        a=${py_ver:0:1} # (2).7
        b=${py_cur:0:1} # (2).7
        if [[ $a = $b ]] ; then
            # echo "0" # success
            echo2 "env: ok, python: ok, ${env_in}"
        else
            echo2 "env: ok, python: failed, ${env_in}" "error"
            # echo "1" # env:exists, python:not_correct
        fi
        deactivate # exit env
    else
        if [[ ${create} = 1 ]] ; then
            env_creater ${env_in} ${py_ver}
            env_checker ${env_in} ${py_ver}
        else
            echo2 "env: failed, python: failed, ${env_in}" "error"
        fi
    fi
}


## main
[[ $# -lt 3 ]] && echo "Usage: $(basename $0) <genome> <bam> <path_out>" && exit 0
genome=$1
inbam=$(readlink -f $2)
path_out=$(readlink -f $3)
[[ -z $(which clipper) ]] && echo2 "clipper - command not found, exiting ..." "error"
[[ -f ${inbam} ]] || echo2 "[${inbam}] file not exists" "error"
[[ ${inbam} = *.bam ]] || echo2 "not a [BAM] file" "error"
#check_genome ${genome} >/dev/null #
[[ ${genome} != hg19 && ${genome} != mm9 ]] && echo2 "${genome} - not supported" "error"
prefix=$(basename ${2%.bam}) # original bam_name
path_sub="${path_out}/${prefix}"
mkdir -p ${path_sub} || echo2 "${path_sub} - cannot create directory" "error"


## switch to py27 env
## ~/envs/clipper_env using python2.7
## using virtualenv to manage env
echo2 "checking virtualenv"
env_checker ${HOME}/envs/clipper_env 2.7
source ${HOME}/envs/clipper_env/bin/activate 

## start
echo2 "start:clipper ${prefix}"
pushd ${path_sub} >/dev/null 2>&1 #

## run clipper
mkdir -p logs
clipperLog="${prefix}.clipper.log"
peakbed="${prefix}.bed"
clipper_para="--bonferroni --superlocal --threshold-method binomial --save-pickle"
if [[ -f ${peakbed} ]] ; then
    echo2 "........ peak exists, clipper skipped ..."
else 
    clipper -b ${inbam} -s ${genome} -o ${peakbed} ${clipper_para} 1>${clipperLog} 2>&1
fi
echo2 "finish:clipper ${prefix}"

popd >/dev/null 2>&1 #

## exit env
deactivate # exit clipper env


## EOF
