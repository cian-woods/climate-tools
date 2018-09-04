#!/bin/bash

git add .
git commit -m "First commit"
git remote add origin https://github.com/cian-woods/climate-tools.git
git remote -v
git push -u origin master
