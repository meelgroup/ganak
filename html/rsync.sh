#!/bin/bash
set -euxo pipefail

rsync -vaP *.js *.wasm index.html version.txt msoos.org:/var/www/ganak/
