#!/bin/bash

############################################################
# Script to download data from S3, run an arbitrary command, then upload the result back to S3
#
## Options
# -i Parameter string defining the source (S3) and destination (path) locations in the form -i <S3 SOURCE URL>|<LOCAL DESTINATION PATH>
# -o Parameter string defining the source (path) and destination (s3) locations in the form -o <LOCAL SOURCE PATH>|<S3 DESTINATION URL>
#
# Example:
# bash run_s3_batch_job.sh \
#   -i "s3://my-bucket/job42/input/input1.txt|/data/input/" \
#   -i "s3://my-bucket/job42/input/input2.txt|/data/input/" \
#   -i "s3://my-bucket/job42/input/more_input/|/data/input/" \
#   -o "/data/output|s3://my-bucket/job/42/output/" \
#   "python3 my_python_script.py --input-path /data/input --output-path /data/output"

# make the script stop when error (non-true exit code) occurs
set -e
START="$(date +%s)"

if ! command -v aws &> /dev/null
then
    cat <<EOF
The AWS CLI could not be found. Please refer to\n
https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html
for more information.
EOF
    exit
fi

unset -v INPUT_PARAMS OUTPUT_PARAMS COMMAND

while getopts "i:o:" option
do
    case $option in
    i) INPUT_PARAMS+=($OPTARG) ;;
    o) OUTPUT_PARAMS=($OPTARG) ;;
    *) exit 1 ;;
    esac
done

shift "$((OPTIND-1))"
COMMAND=$1

if [ -n "$INPUT_PARAMS" ]
then
    IFS="|"
    for val in "${INPUT_PARAMS[@]}";
    do
        read -a strarr <<< "$val"
        if [[ ${strarr[0]} =~ \/$ ]]
        then
            aws s3 cp --recursive ${strarr[0]} ${strarr[1]}
        else
            aws s3 cp ${strarr[0]} ${strarr[1]}
        fi
    done  
fi

eval $COMMAND

if [ -n "$OUTPUT_PARAMS" ]
then
    IFS="|"
    for val in "${OUTPUT_PARAMS[@]}";
    do
        read -a strarr <<< "$val"
        if [[ ${strarr[0]} =~ \/$ ]]
        then
            aws s3 cp --recursive ${strarr[0]} ${strarr[1]}
        else
            aws s3 cp ${strarr[0]} ${strarr[1]}
        fi
    done  
fi

DURATION=$[ $(date +%s) - ${START} ]
echo "Job completed in ${DURATION} seconds."