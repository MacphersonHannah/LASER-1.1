#!/bin/bash
# Download and extract the LASER tutorial data (mouse chr17 subset)
set -euo pipefail

# run from the tutorial/ folder, wherever the script is called from
cd "$(dirname "$0")/.."

URL="https://github.com/MacphersonHannah/LASER-1.1/releases/download/tutorial/laser_tutorial_data.zip"

echo "Downloading tutorial data..."
wget -O laser_tutorial_data.zip "$URL"   # or: curl -L -o laser_tutorial_data.zip "$URL"

echo "Extracting..."
unzip -o laser_tutorial_data.zip
rm laser_tutorial_data.zip

if [ -f data/checksums.md5 ]; then
    echo "Verifying files..."
    (cd data && md5sum -c checksums.md5)
fi

echo "Done. Tutorial data is in $(pwd)/data/"
