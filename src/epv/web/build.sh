#!/usr/bin/env bash

mkdir -p dist
cp LICENSE README.md style.css usage.js visualizer.js dist/
sed -e 's|src="./|src="https://cdn.jsdelivr.net/gh/yuanfeifuzzy/epv@main/src/epv/web/dist/|g' \
    -e 's|href="./|href="https://cdn.jsdelivr.net/gh/yuanfeifuzzy/epv@main/src/epv/web/dist/|g' \
    index.html > dist/index.html
