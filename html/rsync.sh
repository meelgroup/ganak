#!/bin/bash
set -euxo pipefail

rsync -vaP *.js *.wasm  index.html msoos.org:/var/www/ganak/
