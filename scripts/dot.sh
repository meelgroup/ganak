#!/bin/bash
set -e
set -x
./ganak --compile my.ddnnf a.cnf
./bdnnf dot my.ddnnf my.dot
dot -Tpdf my.dot -o my.pdf
dot -Tpng my.dot -o my.png
okular my.pdf
