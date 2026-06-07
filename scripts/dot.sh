#!/bin/bash
./ganak --compile my.ddnnf a.cnf
./ddnf-cleanup my.ddnnf my-cleaned.ddnnf
./ddnnf2dot my-cleaned.ddnnf my.dot
dot -Tpdf my.dot -o my.pdf
okular my.pdf

