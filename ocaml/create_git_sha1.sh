#!/bin/bash

SHA1=$(git log -1  | head -1 | cut -d ' ' -f 2)
DATE=$(git log -1  | grep Date | cut -d ':' -f 2-)
MESSAGE=$(git log -1  | tail -1 | sed 's/"/\\"/g')
cat << EOF > Git.ml
open Core
let sha1 = "$SHA1" |> String_ext.strip
let date = "$DATE" |> String_ext.strip
let message = "$MESSAGE" |> String_ext.strip
EOF

