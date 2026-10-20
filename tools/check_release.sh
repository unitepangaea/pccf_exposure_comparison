#!/bin/sh
# Checks the repository before a commit or release. Exits with an error if:
#   1. any file is larger than 50 MB (GitHub warns at 50 MB and refuses 100 MB)
#   2. any data file would be committed (data/ holds pipeline outputs only)
#   3. any R data file (.rds, .RData) would be committed (record-level intermediate files)
#   4. any file contains a machine-specific path
#   5. any published CSV or Markdown table has a record-level identifier column or more
#      than 50,000 rows (aggregate tables are short)
# Run from the repository root: sh tools/check_release.sh
cd "$(dirname "$0")/.." || exit 1
fail=0
if [ -d .git ]; then FILES=$(git ls-files --cached --others --exclude-standard); else FILES=$(find . -type f ! -path './.git/*' | sed 's|^\./||'); fi

big=$(echo "$FILES" | tr '\n' '\0' | xargs -0 -I{} find {} -size +50M 2>/dev/null)
[ -n "$big" ] && { echo "FAIL files over 50 MB:"; echo "$big"; fail=1; }

tracked=$(echo "$FILES" | grep '^data/' | grep -v '^data/README.md$')
[ -n "$tracked" ] && { echo "FAIL data files would be committed:"; echo "$tracked"; fail=1; }

rds=$(echo "$FILES" | grep -i -E '\.(rds|rdata|rda)$')
[ -n "$rds" ] && { echo "FAIL R data files would be committed:"; echo "$rds"; fail=1; }

paths=$(echo "$FILES" | grep -v '^tools/check_release.sh$' | tr '\n' '\0' | xargs -0 grep -Il -e '/Users/' -e 'Dropbox' -e 'C:\\\\Users' 2>/dev/null)
[ -n "$paths" ] && { echo "FAIL machine-specific paths in:"; echo "$paths"; fail=1; }

ID='(^|,)"?(ADDR_GUID|LOC_GUID|PCODE|PCODE5|PCODE4|DAuid|DAUID|B_DAuid|ref_DAUID|pccf_DAUID|cand_id|x|y|lat|lon|B_x|B_y|pccf_x|pccf_y|pccf_lat|pccf_lon)"?(,|$)'
for f in $(echo "$FILES" | grep -E '^outputs/.*\.csv$'); do
  head -1 "$f" | grep -q -E "$ID" && { echo "FAIL record-level columns in $f"; fail=1; }
  [ "$(wc -l < "$f")" -gt 50000 ] && { echo "FAIL more than 50,000 rows (not an aggregate table?) in $f"; fail=1; }
done
gz=$(echo "$FILES" | grep -E '\.gz$')
[ -n "$gz" ] && { echo "FAIL compressed files would be committed (check by hand):"; echo "$gz"; fail=1; }

[ $fail -eq 0 ] && echo "check_release: OK" || echo "check_release: FAILED"
exit $fail
