#!/usr/bin/env bash

set -e

VERSION=$(git rev-parse --short HEAD)

mkdir -p dist
cp LICENSE README.md style.css usage.js visualizer.js dist/

sed \
  -e "s|src=\"./|src=\"https://cdn.jsdelivr.net/gh/yuanfeifuzzy/epv@${VERSION}/src/epv/web/dist/|g" \
  -e "s|href=\"./|href=\"https://cdn.jsdelivr.net/gh/yuanfeifuzzy/epv@${VERSION}/src/epv/web/dist/|g" \
  index.html > dist/index.html

git add dist
git commit -m "Update dist for ${VERSION}"
git push
