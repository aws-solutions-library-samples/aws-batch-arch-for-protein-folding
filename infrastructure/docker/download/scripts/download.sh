#!/bin/bash
# Download and (optionally) extract files
# Usage: bash download.sh -x myfile.tar.gz destination_folder

if [[ $# -eq 0 ]]; then
    echo "Error: download directory must be provided as an input argument."
    exit 1
fi

if ! command -v aria2c &> /dev/null ; then
    echo "Error: aria2c could not be found. Please install aria2c (sudo apt install aria2)."
    exit 1
fi

if ! command -v aws &> /dev/null ; then
    echo "Error: awscli could not be found. Please refer to https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html"
    exit 1
fi


# make the script stop when error (non-true exit code) is occuredcd
set -e
START="$(date +%s)"

EXTRACT=0
OPTIND=1
while getopts ":x" option
do
    case $option in
    x) EXTRACT=$((EXTRACT+1)) ;; # Extract files after download?
    *) exit 1 ;;
    esac
done

shift "$((OPTIND-1))"

URL=$1
SCHEME=${URL%%:*}
BASENAME=$(basename "${URL}")
DEST="output"

echo "Downloading $URL"
if [[ $SCHEME == "http" ]] ; then
    aria2c "${URL}" --dir="${DEST}"
elif [[ $SCHEME == "s3" ]] ; then
    aws s3 cp "${URL}" "${DEST}/${BASENAME}"
fi

if [[ $EXTRACT -eq 1 ]] && [[ $BASENAME == *"tar"* ]] ; then
    tar --extract --verbose --file="${DEST}/${BASENAME}"  --directory="${DEST}"
    rm "${DEST}/${BASENAME}"
fi
