#!/bin/bash

SHA1=$(git log -1  | head -1 | cut -d ' ' -f 2)
DATE=$(git log -1  | grep Date | cut -d ':' -f 2-)
MESSAGE=$(git log -1  | tail -1 | sed 's/"/\\"/g')
cat << EOF > Git.ml
let sha1 = "$SHA1" |> String.trim
let date = "$DATE" |> String.trim
let message = "$MESSAGE" |> String.trim
EOF

