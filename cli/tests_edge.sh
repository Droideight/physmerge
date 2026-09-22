#!/bin/bash
# Edge cases for the physmerge executable: malformed input, boundary values, and
# the guards that stop a failed or misdirected run from leaving bad files behind.
# Run from the cli directory after `make`. Exits non-zero on the first failure.
set -u
PM=./physmerge
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
fail=0
ok()   { printf '  ok    %s\n' "$1"; }
bad()  { printf '  FAIL  %s\n' "$1"; fail=1; }
check(){ if [ "$2" = "$3" ]; then ok "$1"; else bad "$1 (got '$2', want '$3')"; fi; }
CUST="-f custom --chrom-col #CHROM --pos-col POS --id-col ID --value-col P"

printf '#CHROM\tPOS\tID\tP\r\n1\t1000\trsA\t1e-9\r\n1\t9000\trsB\t1e-9\r\n' > "$T/crlf.tsv"
check "CRLF line endings" "$($PM -i "$T/crlf.tsv" $CUST -w 500 -q | tail -n +2 | wc -l | tr -d ' ')" "2"

printf '#CHROM\tPOS\tID\tP\n1\t1000\trsA\t1e-9\n1\t2000\n1\t9000\trsB\t1e-9\n' > "$T/short.tsv"
check "row with missing fields is dropped" \
  "$($PM -i "$T/short.tsv" $CUST -w 500 -q | tail -n +2 | wc -l | tr -d ' ')" "2"

printf '#CHROM\tPOS\tID\tP\n' > "$T/hdr.tsv"
$PM -i "$T/hdr.tsv" $CUST -q >/dev/null 2>&1
check "header-only file exits 0" "$?" "0"

printf '#CHROM,POS,ID,P\n1,1000,rsA,1e-9\n' > "$T/c.csv"
check "comma separator auto-detected" \
  "$($PM -i "$T/c.csv" $CUST -w 500 -q | tail -n +2 | wc -l | tr -d ' ')" "1"

python3 -c "
print('#CHROM\tPOS\tID\tP\tJUNK')
print('1\t1000\trsA\t1e-9\t'+'X'*1200000)
print('1\t9000\trsB\t1e-9\tshort')" > "$T/long.tsv"
check "line longer than the read buffer" \
  "$($PM -i "$T/long.tsv" $CUST -w 500 -q | tail -n +2 | wc -l | tr -d ' ')" "2"

printf '#CHROM\tPOS\tID\tP\n1\t1000\ta\t1e-9\n2\t1000\tb\t1e-9\n1\t2000\tc\t1e-9\n' > "$T/inter.tsv"
$PM -i "$T/inter.tsv" $CUST -w 500 -q --out "$T/partial.tsv" >/dev/null 2>&1
check "interleaved chromosomes are rejected" "$?" "2"
if [ -e "$T/partial.tsv" ]; then bad "failed run left a partial --out file"; else ok "failed run leaves no --out file"; fi
check "--sort accepts the same file" \
  "$($PM -i "$T/inter.tsv" $CUST -w 500 -q --sort | tail -n +2 | wc -l | tr -d ' ')" "3"

$PM -i "$T/crlf.tsv" $CUST -w 500 -q --out "$T/crlf.tsv" >/dev/null 2>&1
check "--out equal to --input is refused" "$?" "2"
check "input survived that attempt" "$(head -1 "$T/crlf.tsv" | cut -f1)" "#CHROM"

printf '#CHROM\tPOS\tID\tP\n../../escape\t1000\tx\t1e-9\n' > "$T/evil.tsv"
mkdir -p "$T/d"
$PM -i "$T/evil.tsv" $CUST -w 500 -q --snp-list-dir "$T/d" --out /dev/null >/dev/null 2>&1
if [ -e "$T/escape.txt" ] || [ -e "$T/../escape.txt" ]; then bad "chromosome name escaped the output directory"
else ok "chromosome name cannot escape --snp-list-dir"; fi

printf '#CHROM\tPOS\tID\tP\n1\t1000\tx\t5e-8\n' > "$T/eq.tsv"
check "value equal to the threshold is not significant" \
  "$($PM -i "$T/eq.tsv" $CUST -s 5e-8 -w 500 -q | tail -n +2 | wc -l | tr -d ' ')" "0"

check "stdin" "$(cat "$T/crlf.tsv" | $PM -i - $CUST -w 500 -q | tail -n +2 | wc -l | tr -d ' ')" "2"

$PM -i "$T/crlf.tsv" $CUST -w 0 -q >/dev/null 2>&1;   check "--window 0 is refused" "$?" "2"

# --- regressions from tests/qa/FINDINGS.txt (2026-09-22) -----------------------

# C1: a quoted field that contains the separator
printf '#CHROM,POS,ID,P\n1,1000,"rs1,alt",1e-9\n1,9000,rs2,1e-9\n' > "$T/quoted.csv"
check "quoted field containing the separator" \
  "$($PM -i "$T/quoted.csv" $CUST -w 500 -q | awk -F'\t' 'NR==2{print $6}')" "rs1,alt"
printf '#CHROM,POS,ID,P\n1,1000,"say ""hi""",1e-9\n' > "$T/quoted2.csv"
check "doubled quote inside a quoted field" \
  "$($PM -i "$T/quoted2.csv" $CUST -w 500 -q | awk -F'\t' 'NR==2{print $6}')" 'say "hi"'

# C2: a UTF-8 BOM in front of the header
printf '\xef\xbb\xbf#CHROM\tPOS\tID\tP\n1\t1000\trsA\t1e-9\n' > "$T/bom.tsv"
check "UTF-8 BOM in the header" \
  "$($PM -i "$T/bom.tsv" $CUST -w 500 -q | tail -n +2 | wc -l | tr -d ' ')" "1"

# C3: runs of spaces in a whitespace-delimited file
printf '#CHROM POS ID P\n1  1000  rsA 1e-9\n1   9000 rsB 1e-9\n' > "$T/runs.txt"
check "runs of spaces count as one separator" \
  "$($PM -i "$T/runs.txt" $CUST --sep space -w 500 -q | tail -n +2 | wc -l | tr -d ' ')" "2"
printf '#CHROM POS ID P\n   1 1000 rsA 1e-9\n' > "$T/lead.txt"
check "leading whitespace is skipped" \
  "$($PM -i "$T/lead.txt" $CUST --sep space -w 500 -q | tail -n +2 | wc -l | tr -d ' ')" "1"

# C4: two long contig names must not collide on one output file
L=$(python3 -c "print('A'*70)")
printf '#CHROM\tPOS\tID\tP\n%sX\t1000\ta\t1e-9\n%sY\t1000\tb\t1e-9\n' "$L" "$L" > "$T/long2.tsv"
mkdir -p "$T/cd"
$PM -i "$T/long2.tsv" $CUST -w 500 -q --snp-list-dir "$T/cd" -o /dev/null >/dev/null 2>&1
check "long chromosome names do not collide" "$(ls "$T/cd" | wc -l | tr -d ' ')" "2"

# a block's end must never fall below its start
printf '#CHROM\tPOS\tID\tP\n1\t-1000\ta\t1e-9\n1\t-200\tb\t1e-9\n' > "$T/neg.tsv"
check "negative positions give no inverted interval" \
  "$($PM -i "$T/neg.tsv" $CUST -w 500 -q | awk -F'\t' 'NR>1 && $4 < $3 {n++} END{print n+0}')" "0"

# tab files must be untouched by the whitespace rule
printf '#CHROM\tPOS\tID\tP\n1\t1000\t\t1e-9\n' > "$T/empty.tsv"
check "an empty tab-separated field stays empty" \
  "$($PM -i "$T/empty.tsv" $CUST -w 500 -q | awk -F'\t' 'NR==2{print "["$6"]"}')" "[]"

$PM -i "$T/crlf.tsv" $CUST -w abc -q >/dev/null 2>&1; check "non-numeric --window is refused" "$?" "2"
$PM -i "$T/missing_file" $CUST -q >/dev/null 2>&1;    check "missing input is refused" "$?" "2"

[ $fail -eq 0 ] && echo "all edge cases passed" || echo "SOME EDGE CASES FAILED"
exit $fail
