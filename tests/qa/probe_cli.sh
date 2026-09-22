set -u
PM="${PHYSMERGE_PKG:-.}/cli/physmerge"
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
C="-f custom --chrom-col #CHROM --pos-col POS --id-col ID --value-col P"
p() { printf '%s\n' "=== $1 ==="; }

p "C1 quoted CSV field (fread handles quotes, C does not)"
printf '#CHROM,POS,ID,P\n1,1000,"rs1,alt",1e-9\n1,9000,rs2,1e-9\n' > "$T/q.csv"
"$PM" -i "$T/q.csv" $C -w 500 -q 2>&1

p "C2 BOM in header"
printf '\xef\xbb\xbf#CHROM\tPOS\tID\tP\n1\t1000\trs1\t1e-9\n' > "$T/bom.tsv"
"$PM" -i "$T/bom.tsv" $C -w 500 -q 2>&1 | head -3

p "C3 blank lines inside the file"
printf '#CHROM\tPOS\tID\tP\n1\t1000\trs1\t1e-9\n\n1\t9000\trs2\t1e-9\n' > "$T/bl.tsv"
"$PM" -i "$T/bl.tsv" $C -w 500 -q 2>&1

p "C4 big positions + snp-list falls back to rps_BP"
printf '#CHROM\tPOS\tP\n1\t900000\t1e-9\n1\t248000000\t1e-9\n' > "$T/big.tsv"
"$PM" -i "$T/big.tsv" -f custom --chrom-col '#CHROM' --pos-col POS --id-col NA --value-col P -w 500 -q --snp-list "$T/s.txt" ; echo "-- snp list:"; cat "$T/s.txt"

p "C5 p-value underflow 1e-400 and literal 0"
printf '#CHROM\tPOS\tID\tP\n1\t1000\ta\t1e-400\n1\t9000\tb\t0\n' > "$T/uf.tsv"
"$PM" -i "$T/uf.tsv" $C -w 500 -q 2>&1

p "C6 non-numeric --sig-th nan / inf"
"$PM" -i "$T/bl.tsv" $C -w 500 -s nan -q 2>&1 | head -3; echo "exit=$?"
"$PM" -i "$T/bl.tsv" $C -w 500 -s inf -q 2>&1 | head -3

p "C7 fractional positions"
printf '#CHROM\tPOS\tID\tP\n1\t1000.5\ta\t1e-9\n1\t1600.5\tb\t1e-9\n' > "$T/fr.tsv"
"$PM" -i "$T/fr.tsv" $C -w 500 -q 2>&1

p "C8 --no-chrom on a file that has chromosomes"
printf '#CHROM\tPOS\tID\tP\n1\t1000\ta\t1e-9\n2\t1100\tb\t1e-9\n' > "$T/nc.tsv"
"$PM" -i "$T/nc.tsv" -f custom --pos-col POS --id-col ID --value-col P --no-chrom -w 500 -q 2>&1

p "C9 --chrom filter"
"$PM" -i "$T/nc.tsv" $C --chrom 2 -w 500 -q 2>&1

p "C10 long chromosome names collide in --snp-list-dir"
LONG1=$(python3 -c "print('A'*70+'X')"); LONG2=$(python3 -c "print('A'*70+'Y')")
printf '#CHROM\tPOS\tID\tP\n%s\t1000\ta\t1e-9\n%s\t1000\tb\t1e-9\n' "$LONG1" "$LONG2" > "$T/long.tsv"
mkdir -p "$T/dir"; "$PM" -i "$T/long.tsv" $C -w 500 -q --snp-list-dir "$T/dir" -o /dev/null 2>&1
echo "files:"; ls "$T/dir"; echo "contents:"; cat "$T/dir"/*

p "C11 --snp-list-dir that does not exist"
"$PM" -i "$T/nc.tsv" $C -w 500 -q --snp-list-dir "$T/missing" -o /dev/null 2>&1; echo "exit=$?"

p "C12 duplicate column names in header"
printf '#CHROM\tPOS\tID\tP\tP\n1\t1000\ta\t1e-9\t0.9\n' > "$T/dup.tsv"
"$PM" -i "$T/dup.tsv" $C -w 500 -q 2>&1

p "C13 --annotate-full"
"$PM" -i "$T/nc.tsv" $C -w 500 -q --annotate-full 2>&1

p "C14 trailing whitespace / space-separated with runs of spaces"
printf '#CHROM POS ID P\n1  1000 a 1e-9\n' > "$T/sp.txt"
"$PM" -i "$T/sp.txt" $C -w 500 -q 2>&1

p "C15 --out into a directory"
"$PM" -i "$T/nc.tsv" $C -w 500 -q -o "$T" 2>&1; echo "exit=$?"

p "C16 negative positions"
printf '#CHROM\tPOS\tID\tP\n1\t-1000\ta\t1e-9\n1\t-200\tb\t1e-9\n' > "$T/neg.tsv"
"$PM" -i "$T/neg.tsv" $C -w 500 -q 2>&1

p "C17 --sort with 1e6 rows memory sanity (timing)"
python3 -c "
import random
print('#CHROM\tPOS\tID\tP')
for i in range(300000): print('1\t%d\trs%d\t%g'%(i*10, i, random.random()))" > "$T/perf.tsv"
/usr/bin/time -l "$PM" -i "$T/perf.tsv" $C -w 500000 -s 1e-4 -q -o /dev/null 2>&1 | grep -E "real|maximum resident" | head -3
