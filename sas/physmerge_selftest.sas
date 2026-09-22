/*===========================================================================
  physmerge for SAS -- self-test

  Runs %physmerge over the vectors in sas/testdata/ and compares the block
  table against sas/testdata/vectors_expected.tsv, which was produced by the
  R package.  A clean run prints "SELFTEST: PASS" and every PROC COMPARE
  reports NOTE: No unequal values were found.

  Usage
    %let PMDIR = /path/to/physmerge_pkg/sas;   * <- edit this one line ;
    %include "&PMDIR/physmerge_selftest.sas";
===========================================================================*/

%if %symexist(PMDIR) = 0 %then %do;
  %put ERROR: set %nrstr(%let PMDIR = ...;) to the sas/ directory first.;
%end;

%include "&PMDIR/physmerge.sas";
%pm_version;

/*--------------------------------------------------------------- vectors --*/
proc import datafile="&PMDIR/testdata/vectors_input.tsv" out=vin dbms=dlm replace;
  delimiter='09'x; getnames=yes; guessingrows=max;
run;
proc import datafile="&PMDIR/testdata/vectors_params.tsv" out=vpar dbms=dlm replace;
  delimiter='09'x; getnames=yes; guessingrows=max;
run;
proc import datafile="&PMDIR/testdata/vectors_expected.tsv" out=vexp_raw dbms=dlm replace;
  delimiter='09'x; getnames=yes; guessingrows=max;
run;

/* P and rps_value arrive as text in full precision; parse them explicitly so
   nothing depends on how PROC IMPORT guessed the type */
data vin;
  set vin;
  length case $ 32 CHROM $ 32 ID $ 64;
  case  = cats(case);
  CHROM = cats(CHROM);
  ID    = cats(ID);
  POSITION = input(cats(POS), ?? best32.);
  VALUE    = input(cats(P),   ?? best32.);
  keep case CHROM POS ID P POSITION VALUE;
run;

data vexp;
  length case $ 32 CHROM $ 32 rps_ID $ 64;
  set vexp_raw;
  case   = cats(case);
  CHROM  = cats(CHROM);
  rps_ID = cats(rps_ID);
  start  = input(cats(start),  ?? best32.);
  end    = input(cats(end),    ?? best32.);
  rps_BP = input(cats(rps_BP), ?? best32.);
  rps_value = input(cats(rps_value), ?? best32.);
  keep case serial CHROM start end rps_BP rps_ID rps_value;
run;

/*----------------------------------------------------------- run each case */
proc datasets lib=work nolist nowarn; delete allblocks; quit;

%macro pm_selftest;
  %local ncase i cs sg wn rw ro;
  proc sql noprint;
    select count(*) into :ncase trimmed from vpar;
  quit;

  %do i = 1 %to &ncase;
    proc sql noprint;
      select cats(case), cats(sig_th), cats(window), cats(reward), cats(reset_on)
        into :cs trimmed, :sg trimmed, :wn trimmed, :rw trimmed, :ro trimmed
        from vpar(firstobs=&i obs=&i);
    quit;

    data _one;
      set vin;
      where case = "&cs";
    run;

    %physmerge(data=_one, out=_blk, sig_th=&sg, window=&wn,
               reward=&rw, reset_on=&ro,
               chrom=CHROM, pos=POSITION, value=VALUE, id=ID, quiet=1);

    data _blk;
      length case $ 32;
      set _blk;
      case = "&cs";
    run;

    proc append base=allblocks data=_blk force; run;
  %end;

  proc sort data=allblocks; by case serial; run;
  proc sort data=vexp;      by case serial; run;

  /* exact on the integers and the id, 1e-12 relative on the value */
  proc compare base=vexp compare=allblocks out=_cmp outnoequal noprint
               criterion=1e-12 method=relative;
    id case serial;
    var CHROM start end rps_BP rps_ID rps_value;
  run;

  %local ndiff nb nc;
  %let ndiff = %pm_nobs(_cmp);
  %let nb    = %pm_nobs(vexp);
  %let nc    = %pm_nobs(allblocks);

  %put NOTE: expected blocks = &nb, produced blocks = &nc, unequal rows = &ndiff;
  %if &nb = &nc and &ndiff = 0 %then %put SELFTEST: PASS (&ncase cases, &nb blocks).;
  %else %do;
    %put SELFTEST: FAIL -- see the tables below.;
    proc print data=_cmp(obs=40) noobs; title "rows that differ"; run;
    title;
  %end;
%mend pm_selftest;

%pm_selftest;

/*------------------------------------------------- guards that should fail */
%put NOTE: --- the next three lines must each print an ERROR ---;
%physmerge(data=vin, out=_x, window=0,       chrom=CHROM, id=ID, quiet=1);
%physmerge(data=vin, out=_x, reward=sideways,chrom=CHROM, id=ID, quiet=1);
%physmerge(data=vin, out=_x, reset_on=maybe, chrom=CHROM, id=ID, quiet=1);

/*---------------------------------------------- missing values must not win */
data _miss;
  length CHROM $ 32 ID $ 64;
  input CHROM $ POSITION VALUE ID $;
datalines;
1 1000 . rsNA
1 2000 1e-9 rsHIT
;
run;
%physmerge(data=_miss, out=_mblk, sig_th=5e-8, window=500, chrom=CHROM, id=ID, quiet=1);
%put NOTE: missing-value test: expect exactly 1 block with rps_ID=rsHIT;
proc print data=_mblk noobs; run;

/*------------------------------------------------ export must not use e+06 */
%pm_export(blocks=_mblk, file="&PMDIR/testdata/_selftest_ids.txt");
data _null_;
  infile "&PMDIR/testdata/_selftest_ids.txt" truncover;
  input line $200.;
  put "NOTE: exported id: " line;
run;
