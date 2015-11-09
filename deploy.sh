#!/bin/bash

set -o errexit -o nounset

rev=$(git rev-parse --short HEAD)

npm run build
npm run render

git init
git config user.name "Bruno Vieira"
git config user.email "mail@bmpvieira.com"

git remote add upstream "https://$GH_TOKEN@github.com/bionode/get-bionode.git"
git fetch upstream
git reset upstream/gh-pages

touch .

git add -A .
git commit -m "rebuild pages at ${rev}"
git push -q upstream HEAD:gh-pages
