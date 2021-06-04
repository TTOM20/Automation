#!/bin/bash
# #OCR
find "./img/" -type f -name "*.png" | sed 's/\.png$//' | xargs -P8 -n1 -I% tesseract %.png % -l eng+jpn pdf

# #sum
pdftk ./img/*.pdf cat output "./result/a.pdf"


