#!/usr/bin/env bash

NAME="pierre-24"
REPO="osmipy"
HTML_DOC_DIR="documentation/build/html"

# config GIT
git config --global user.email "travis@travis-ci.org"
git config --global user.name "travis-ci"

# make doc
make doc

# push doc
git add $HTML_DOC_DIR -f
git commit -m "Deploy"
git subtree split --branch build_doc --prefix $HTML_DOC_DIR
git push https://$NAME:$GITHUB_API_KEY@github.com/$NAME/$REPO build_doc:gh-pages -fq > /dev/null 2>&1