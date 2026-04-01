#!/bin/bash
set -euxo pipefail

rsync -vaP ganak.js ganak.wasm  index.html msoos.org:/var/www/ganak/
