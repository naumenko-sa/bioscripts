#!/bin/bash

#https://www.quora.com/How-do-I-download-a-very-large-file-from-Google-Drive
# get link in browser
# https://drive.google.com/file/d/XXXXX/view?usp=sharing
id=XXXXX
# get token from https://developers.google.com/oauthplayground/
# Drive API v3, and select https://www.googleapis.com/auth/drive.readonly, “Authorize APIs” and then “Exchange authorization code for tokens
token=YYYYY
outfile=ZZZZ.tar.gz

curl -H "Authorization: Bearer "$token "https://www.googleapis.com/drive/v3/files/"${id}"?alt=media" -o $outfile

# some more methods
# https://gist.github.com/iamtekeste/3cdfd0366ebfd2c0d805
