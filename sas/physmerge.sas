/*===========================================================================
  physmerge for SAS  --  panel-free physical locus merging
  ---------------------------------------------------------------------------
  Port of the R package (R/read_sumstat.R, R/physical_merge.R, R/export.R) and
  of the C command-line tool (cli/physmerge.c).  The three implementations are
  held to the same block table: sas/physmerge_selftest.sas checks this one
  against vectors produced by the R code.

  Macros
    %pm_read      read a summary-statistics file into a SAS data set
    %physmerge    forward sliding-window merge  -> block table
    %pm_annotate  attach the lead SNP's original columns to the block table
    %pm_export    write representative SNP ids, one file or one file per chrom
    %pm_version   print the version

  Helpers (called by the above, usable on their own)
    %pm_col       resolve a file header name to its SAS variable name
    %pm_vtype     variable type, C or N
    %pm_prep      normalise / filter / sort an input data set

  Base SAS only: no SAS/STAT, no SAS/ACCESS, no PROC FCMP.  Written against
  SAS 9.4; nothing here is 9.4-specific, so Viya and OnDemand should also run
  it.

  STATUS.  This has not yet been executed on a SAS installation.  The algorithm
  is a line-by-line transliteration of R/physical_merge.R, and the same
  transliteration into C agrees with R on 400 randomised inputs; what is
  unverified is the SAS syntax.  Run physmerge_selftest.sas first: it merges 13
  cases and compares them against sas/testdata/, which the R package produced,
  and prints "SELFTEST: PASS" when the port is faithful.

  THREE THINGS THAT DIFFER FROM R, ON PURPOSE
    Missing values.  A SAS missing numeric compares below every number, so
      `value < 5e-8` is TRUE for `.`.  %pm_read and %pm_prep drop rows with a
      missing position or value before the scan.  Do not hand %physmerge a data
      set you assembled without that filter.
    Chromosome order.  Chromosomes are ranked by first appearance, not
      alphabetically, so the block serial numbers match physical_merge(), which
      walks unique(data$CHROM).  %pm_prep sorts by (rank, position, input row),
      reproducing R's stable order().
    Exported ids.  A numeric id is written with BEST32., so position 900000
      comes out as 900000 rather than R's "9e+05", which PLINK cannot use.
      %pm_export(by_chrom=1) writes plain snp_ch<CHR>.txt files into dir=; it
      does not zip them, because Base SAS has no portable zip.

  NUMERIC RANGE -- READ THIS BEFORE MERGING ON P-VALUES
    The BEST32. informat does not read subnormals reliably, and on z/OS the
    native floating-point format underflows around 1e-78.  GWAS p-values
    routinely go below 1e-300.  A p-value that underflows to 0 is still
    "significant", so the merge is usually still correct, but rps_value reads 0
    and any downstream ranking on it breaks.  If your summary statistics reach
    that range, merge on LOG10_P instead:

      %physmerge(data=ss, out=blocks, value=LOG10_P,
                 sig_th=7.3, reward=max, window=500000, ...);

    -log10(5e-8) = 7.301.  This is the same advice read_sumstat() gives in R.

  Compressed input.  %pm_read cannot read .gz; gunzip first, or use the C tool,
  which reads gzip natively.
===========================================================================*/

%macro pm_version;
  %put NOTE: physmerge for SAS 0.3.0 (sas-1), matching physmerge R 0.3.0.;
%mend pm_version;


/*---------------------------------------------------------------------------
  %pm_vtype(ds, var) -> C or N (empty if the variable does not exist)
---------------------------------------------------------------------------*/
%macro pm_vtype(ds, var);
  %local dsid vnum t rc;
  %let dsid = %sysfunc(open(&ds));
  %if &dsid %then %do;
    %let vnum = %sysfunc(varnum(&dsid, &var));
    %if &vnum %then %let t = %sysfunc(vartype(&dsid, &vnum));
    %let rc = %sysfunc(close(&dsid));
  %end;
  &t
%mend pm_vtype;


/*---------------------------------------------------------------------------
  %pm_col(ds, want) -> the SAS variable holding the file column named `want`.

  PROC IMPORT under VALIDVARNAME=V7 rewrites a header such as "#CHROM" to the
  variable _CHROM and keeps "#CHROM" as the label.  This looks at the label
  first, then the name, then the V7 rewrite of `want`, so both spellings work.
  Returns empty when nothing matches.
---------------------------------------------------------------------------*/
%macro pm_col(ds, want);
  %local lib mem hit;
  %if %index(&ds, .) %then %do;
    %let lib = %upcase(%scan(&ds, 1, .));
    %let mem = %upcase(%scan(&ds, 2, .));
  %end;
  %else %do;
    %let lib = WORK;
    %let mem = %upcase(&ds);
  %end;

  %global pm__want;
  %let pm__want = %superq(want);

  proc sql noprint;
    select name into :hit trimmed from dictionary.columns
     where libname = "&lib" and memname = "&mem"
       and (upcase(label) = upcase(symget('pm__want'))
         or upcase(name)  = upcase(symget('pm__want'))
         /* PROC IMPORT under VALIDVARNAME=V7 rewrites every character that is
            not A-Z 0-9 _ to _, and prefixes _ when the result does not start
            with a letter or _ : "#CHROM" becomes _CHROM */
         or upcase(name)  = upcase(prxchange('s/[^A-Za-z0-9_]/_/', -1,
                                             strip(symget('pm__want'))))
         or upcase(name)  = upcase(cats('_', prxchange('s/[^A-Za-z0-9_]/_/', -1,
                                             strip(symget('pm__want'))))))
     order by case when upcase(label) = upcase(symget('pm__want')) then 1
                   when upcase(name)  = upcase(symget('pm__want')) then 2
                   else 3 end;
  quit;
  &hit
%mend pm_col;


/*---------------------------------------------------------------------------
  %pm_read  --  the %pm_prep-ready reader, mirroring read_sumstat()

    path=         delimited text file (plain text; gunzip .gz first)
    out=          output data set                      (default pm_sumstat)
    format=       plink2 | gpcm | custom               (default plink2)
    chrom_col=    default #CHROM for plink2 and gpcm
    pos_col=      default POS
    id_col=       default ID; pass id_col=_none_ when the file has no id
    value_col=    default P (plink2) or P_HPI (gpcm)
    test_filter=  1 | 0; default 1 for plink2, 0 otherwise
    test_col=     default TEST
    test_val=     default ADD
    chrom=        comma-separated list to keep, e.g. chrom=%str('1','2')
    dlm=          delimiter; default '09'x (tab).  Use ',' for csv.
    guessingrows= PROC IMPORT GUESSINGROWS (default MAX)

  The output data set keeps every original column and adds POSITION and VALUE,
  the two numeric columns %physmerge works from.  Rows with a missing POSITION
  or VALUE are dropped, as read_sumstat() drops them.
---------------------------------------------------------------------------*/
%macro pm_read(path=, out=pm_sumstat, format=plink2,
               chrom_col=, pos_col=, id_col=, value_col=,
               test_filter=, test_col=TEST, test_val=ADD, chrom=,
               dlm='09'x, guessingrows=MAX, quiet=0);

  %local fmt vc vp vv vi vt n0 n1 save_vvn;
  %let fmt = %upcase(&format);
  %let save_vvn = %sysfunc(getoption(validvarname));

  /* ---- format defaults, as in read_sumstat() ---- */
  %if &fmt = PLINK2 %then %do;
    %if %length(&chrom_col) = 0 %then %let chrom_col = %str(#CHROM);
    %if %length(&pos_col)   = 0 %then %let pos_col   = POS;
    %if %length(&id_col)    = 0 %then %let id_col    = ID;
    %if %length(&value_col) = 0 %then %let value_col = P;
    %if %length(&test_filter) = 0 %then %let test_filter = 1;
  %end;
  %else %if &fmt = GPCM %then %do;
    %if %length(&chrom_col) = 0 %then %let chrom_col = %str(#CHROM);
    %if %length(&pos_col)   = 0 %then %let pos_col   = POS;
    %if %length(&id_col)    = 0 %then %let id_col    = ID;
    %if %length(&value_col) = 0 %then %let value_col = P_HPI;
    %if %length(&test_filter) = 0 %then %let test_filter = 0;
  %end;
  %else %if &fmt = CUSTOM %then %do;
    %if %length(&test_filter) = 0 %then %let test_filter = 0;
    %if %length(&chrom_col) = 0 or %length(&pos_col) = 0 or %length(&value_col) = 0 %then %do;
      %put ERROR: For format = 'custom', you must supply chrom_col=, pos_col= and value_col=.;
      %return;
    %end;
  %end;
  %else %do;
    %put ERROR: format= must be plink2, gpcm or custom.;
    %return;
  %end;

  options validvarname=v7;
  proc import datafile="&path" out=_pm_imp dbms=dlm replace;
    delimiter = &dlm;
    getnames  = yes;
    guessingrows = &guessingrows;
  run;
  %if &syserr > 1 %then %do;
    %put ERROR: Failed to read file: &path;
    options validvarname=&save_vvn;
    %return;
  %end;

  %let vc = %pm_col(_pm_imp, &chrom_col);
  %let vp = %pm_col(_pm_imp, &pos_col);
  %let vv = %pm_col(_pm_imp, &value_col);
  %let vt = %pm_col(_pm_imp, &test_col);
  %if %upcase(&id_col) = _NONE_ %then %let vi = ;
  %else %let vi = %pm_col(_pm_imp, &id_col);

  %local miss;
  %let miss = ;
  %if %length(&vc) = 0 %then %let miss = &miss &chrom_col;
  %if %length(&vp) = 0 %then %let miss = &miss &pos_col;
  %if %length(&vv) = 0 %then %let miss = &miss &value_col;
  %if %upcase(&id_col) ne _NONE_ and %length(&vi) = 0 %then %let miss = &miss &id_col;
  %if %length(&miss) %then %do;
    %put ERROR: Column(s) not found:&miss;
    options validvarname=&save_vvn;
    %return;
  %end;

  /* ---- TEST filter ---- */
  %if &test_filter = 1 %then %do;
    %if %length(&vt) = 0 %then %do;
      %put WARNING: test_col '&test_col' not found; TEST filter skipped.;
    %end;
    %else %do;
      %let n0 = %pm_nobs(_pm_imp);
      data _pm_imp;
        set _pm_imp;
        where upcase(cats(&vt)) = %upcase("&test_val");
      run;
      %let n1 = %pm_nobs(_pm_imp);
      %if &quiet = 0 %then
        %put NOTE: TEST filter: kept &n1 of &n0 rows where &test_col = "&test_val".;
      %if &n1 = 0 %then %do;
        %put ERROR: No rows remain after TEST filter.;
        options validvarname=&save_vvn;
        %return;
      %end;
    %end;
  %end;

  /* ---- interface columns.  PROC IMPORT turns a column that contains "NA"
          into a character variable, so always route through input(). ---- */
  data &out;
    set _pm_imp;
    %if %pm_vtype(_pm_imp, &vp) = C %then %do;
      POSITION = input(cats(&vp), ?? best32.);
    %end;
    %else %do;
      POSITION = &vp;
    %end;
    %if %pm_vtype(_pm_imp, &vv) = C %then %do;
      VALUE = input(cats(&vv), ?? best32.);
    %end;
    %else %do;
      VALUE = &vv;
    %end;
    %if %length(&chrom) %then %do;
      if cats(&vc) not in (&chrom) then delete;
    %end;
    if missing(POSITION) or missing(VALUE) then do; _pm_na + 1; delete; end;
    drop _pm_na;
  run;

  %if %pm_nobs(&out) = 0 %then %put ERROR: No usable rows after filtering.;

  options validvarname=&save_vvn;

  /* remember the resolved names so %physmerge can default to them */
  %global pm_chrom_var pm_pos_var pm_id_var pm_value_var pm_reward;
  %let pm_chrom_var = &vc;
  %let pm_pos_var   = &vp;
  %let pm_id_var    = &vi;
  %let pm_value_var = &vv;
  %if %index(%str( LOG10_P T_STAT Z_STAT CHISQ F_STAT T_STAT_DIRECT T_STAT_TE HPI ),
             %str( )%upcase(&value_col)%str( )) %then %let pm_reward = max;
  %else %let pm_reward = min;
  %if &quiet = 0 %then
    %put NOTE: physmerge: suggested reward = &pm_reward for value column &value_col..;

  proc datasets lib=work nolist nowarn; delete _pm_imp; quit;
%mend pm_read;


/*---------------------------------------------------------------------------
  %pm_nobs(ds) -> observation count
---------------------------------------------------------------------------*/
%macro pm_nobs(ds);
  %local dsid n rc;
  %let n = 0;
  %let dsid = %sysfunc(open(&ds));
  %if &dsid %then %do;
    %let n  = %sysfunc(attrn(&dsid, nlobs));
    %let rc = %sysfunc(close(&dsid));
  %end;
  &n
%mend pm_nobs;


/*---------------------------------------------------------------------------
  %pm_prep  --  normalise, drop missing rows, and sort the way R's order() does

  Chromosomes are ranked by first appearance, not alphabetically, so that the
  block serial numbers line up with physical_merge()'s, which walks
  unique(data$CHROM).  Within a chromosome the sort is (position, input row),
  reproducing R's stable order().
---------------------------------------------------------------------------*/
%macro pm_prep(data=, out=_pm_prep, chrom=, pos=POSITION, value=VALUE, id=,
               idlen=200, chromlen=32);
  data &out(keep=_chrord _chrom _pos _val _id _seq);
    length _chrom $ &chromlen _id $ &idlen _chrord 8;
    if _n_ = 1 then do;
      declare hash _h();
      _h.defineKey('_chrom');
      _h.defineData('_chrord');
      _h.defineDone();
    end;
    retain _nchr 0;
    set &data;
    %if %length(&chrom) %then %do; _chrom = cats(&chrom); %end;
    %else %do;                     _chrom = '';           %end;
    %if %length(&id) %then %do;    _id    = cats(&id);    %end;
    %else %do;                     _id    = '';           %end;
    _pos = &pos;
    _val = &value;
    /* a SAS missing value compares low, so it would look significant under
       reward=min: drop those rows before they reach the scan */
    if missing(_pos) or missing(_val) then delete;
    _seq + 1;
    if _h.find() ne 0 then do;
      _nchr + 1;
      _chrord = _nchr;
      _h.add();
    end;
    drop _nchr;
  run;
  proc sort data=&out; by _chrord _pos _seq; run;
%mend pm_prep;


/*---------------------------------------------------------------------------
  %physmerge  --  forward sliding-window merge

    data=       input data set
    out=        block table                            (default pm_blocks)
    sig_th=     significance threshold                 (default 5e-8)
    window=     window in bp                           (default 500000)
    reward=     min for p-values (default) | max for test statistics
    reset_on=   best (default) | any
    chrom=      chromosome variable; leave blank for a single-chromosome run
    pos=        position variable                      (default POSITION)
    value=      value variable                         (default VALUE)
    id=         lead-SNP id variable, optional
    idlen=      character length for the id            (default 200)

  Output columns: serial CHROM start end rps_BP rps_ID rps_value, the same set
  and the same order as annotate_blocks() in R and as the C tool's table.
---------------------------------------------------------------------------*/
%macro physmerge(data=, out=pm_blocks, sig_th=5e-8, window=500000,
                 reward=min, reset_on=best,
                 chrom=, pos=POSITION, value=VALUE, id=,
                 idlen=200, chromlen=32, quiet=0);

  %local c rw ro nin nout;
  %let rw = %upcase(&reward);
  %let ro = %upcase(&reset_on);
  %if %length(&rw) = 0 %then %let rw = .;
  %if %length(&ro) = 0 %then %let ro = .;

  %if &rw = MAX %then %let c = >;
  %else %if &rw = MIN %then %let c = <;
  %else %do; %put ERROR: `reward` must be either 'min' or 'max'.; %return; %end;

  %if &ro ne BEST and &ro ne ANY %then %do;
    %put ERROR: `reset_on` must be either 'best' or 'any'.; %return;
  %end;
  %if %sysevalf(&window <= 0) %then %do;
    %put ERROR: `window` must be a single positive numeric value.; %return;
  %end;
  %if %length(%superq(sig_th)) = 0 %then %do;
    %put ERROR: `sig_th` must be a single numeric value.; %return;
  %end;

  %pm_prep(data=&data, out=_pm_prep, chrom=&chrom, pos=&pos, value=&value,
           id=&id, idlen=&idlen, chromlen=&chromlen);
  %let nin = %pm_nobs(_pm_prep);

  /* ---- pass 1: the forward scan (physical_merge's for loop) -------------- */
  data _pm_raw(keep=_chrord _chrom start end rps_BP rps_value rps_ID);
    length _chrom $ &chromlen rps_ID $ &idlen;
    retain in_block steps sig_this last_pos start end rps_BP rps_value rps_ID;
    /* _chrom is read from the input record; blocks never span a chromosome */
    set _pm_prep;
    by _chrord;

    if first._chrord then do;
      in_block = 0;
      steps    = &window;
      sig_this = &sig_th;
      last_pos = _pos;
    end;

    if not in_block then do;
      if _val &c &sig_th then link pm_open;
    end;
    else do;
      remaining = steps - (_pos - last_pos);
      if remaining <= 0 then do;
        link pm_close;
        if _val &c &sig_th then link pm_open;
      end;
      else do;
        steps = remaining;
        %if &ro = ANY %then %do;
          /* locusDefiner style: any significant SNP resets the window */
          if _val &c &sig_th then do;
            steps = &window;
            if _val &c sig_this then link pm_rep;
          end;
          else if _val &c sig_this then do;
            steps = &window;
            link pm_rep;
          end;
        %end;
        %else %do;
          /* reset_on = best: reset only on improvement */
          if _val &c sig_this then do;
            steps = &window;
            link pm_rep;
          end;
        %end;
      end;
    end;

    last_pos = _pos;
    if last._chrord and in_block then link pm_close;
    return;

  pm_open:
    start     = max(0, _pos - &window);
    end       = .;
    rps_BP    = _pos;
    rps_value = _val;
    rps_ID    = _id;
    in_block  = 1;
    steps     = &window;
    sig_this  = _val;
    return;

  pm_rep:
    sig_this  = _val;
    rps_BP    = _pos;
    rps_value = _val;
    rps_ID    = _id;
    return;

  pm_close:
    end      = last_pos + steps;
    output;
    in_block = 0;
    steps    = &window;
    sig_this = &sig_th;
    return;
  run;

  /* ---- pass 2: collapse, then trim, with a one-block lookahead ----------- */
  data _pm_col(keep=_chrord h_chrom h_start h_end h_bp h_val h_id);
    length h_chrom $ &chromlen h_id $ &idlen;
    retain has_held h_chrom h_start h_end h_bp h_val h_id;
    set _pm_raw(rename=(_chrom=i_chrom start=i_start end=i_end
                        rps_BP=i_bp rps_value=i_val rps_ID=i_id));
    by _chrord;

    if first._chrord then has_held = 0;

    if has_held then do;
      if (i_bp - h_bp) < &window then do;          /* .collapse_blocks() */
        if i_end > h_end then h_end = i_end;
        if i_val &c h_val then do;
          h_bp = i_bp; h_val = i_val; h_id = i_id;
        end;
      end;
      else do;
        if h_end > i_start then h_end = i_start;   /* trim pass */
        output;
        link pm_hold;
      end;
    end;
    else link pm_hold;

    if last._chrord and has_held then output;
    return;

  pm_hold:
    h_chrom = i_chrom; h_start = i_start; h_end = i_end;
    h_bp    = i_bp;    h_val   = i_val;   h_id  = i_id;
    has_held = 1;
    return;
  run;

  /* ---- pass 3: serial numbers and the public column set ------------------ */
  data &out(keep=serial CHROM start end rps_BP rps_ID rps_value);
    length serial 8 CHROM $ &chromlen start 8 end 8 rps_BP 8
           rps_ID $ &idlen rps_value 8;
    set _pm_col;
    serial + 1;
    CHROM     = h_chrom;
    start     = h_start;
    /* start is clamped at 0 by max(0, pos - window) but end is not, so a
       negative position would otherwise invert the interval.  For a
       non-negative position this is a no-op. */
    end       = max(h_end, h_start);
    rps_BP    = h_bp;
    rps_ID    = h_id;
    rps_value = h_val;
  run;

  %let nout = %pm_nobs(&out);
  %if &quiet = 0 %then
    %put NOTE: physmerge: &nin SNPs -> &nout blocks (window=&window, sig_th=&sig_th, reward=&reward, reset_on=&reset_on).;

  proc datasets lib=work nolist nowarn; delete _pm_prep _pm_raw _pm_col; quit;
%mend physmerge;


/*---------------------------------------------------------------------------
  %pm_annotate  --  attach the lead SNP's original columns (annotate_blocks)

  Joins the block table back to the input on (CHROM, rps_BP).  When several
  variants share one base pair the join would be ambiguous, so it also matches
  on the lead SNP id when the block table carries one.
---------------------------------------------------------------------------*/
%macro pm_annotate(blocks=pm_blocks, data=, out=pm_blocks_annot,
                   chrom=, pos=POSITION, id=, keep=);

  %local addcols;
  /* every input column except the join keys and the interface columns */
  proc sql noprint;
    select strip(name) into :addcols separated by ' '
      from dictionary.columns
     where libname = 'WORK' and memname = %upcase("&data")
       and upcase(name) not in ('POSITION', 'VALUE'
            %if %length(&chrom) %then , %upcase("&chrom") ;
            %if %length(&pos)   %then , %upcase("&pos")   ;
            %if %length(&id)    %then , %upcase("&id")    ; );
  quit;
  %if %length(&keep) %then %let addcols = &keep;

  proc sql;
    create table &out as
    select b.*
           %if %length(&addcols) %then %do;
             %local i col;
             %do i = 1 %to %sysfunc(countw(&addcols));
               %let col = %scan(&addcols, &i);
               , d.&col
             %end;
           %end;
      from &blocks as b
      left join &data as d
        on d.&pos = b.rps_BP
           %if %length(&chrom) %then and cats(d.&chrom) = b.CHROM ;
           %if %length(&id)    %then and cats(d.&id)    = b.rps_ID ;
     order by b.serial;
  quit;
%mend pm_annotate;


/*---------------------------------------------------------------------------
  %pm_export  --  representative SNP ids, one per line (export_snp_list)

    blocks=    block table
    file=      output file, when by_chrom=0
    dir=       output directory, when by_chrom=1; writes snp_ch<CHR>.txt
    id=        column to write; rps_ID by default, rps_BP when there is no id
    by_chrom=  0 (default) | 1

  A numeric id is written with BEST32., never in scientific notation, so the
  file stays usable as a PLINK --extract list.  As in the C tool, a chromosome
  name is sanitised before it becomes part of a file name.
---------------------------------------------------------------------------*/
%macro pm_export(blocks=pm_blocks, file=, dir=, id=, by_chrom=0);
  %local v t;
  %if %length(&id) %then %let v = &id;
  %else %if %length(%pm_vtype(&blocks, rps_ID)) %then %let v = rps_ID;
  %else %let v = rps_BP;

  %let t = %pm_vtype(&blocks, &v);
  %if %length(&t) = 0 %then %do;
    %put ERROR: Column '&v' not found in &blocks..;
    %return;
  %end;
  %if %pm_nobs(&blocks) = 0 %then %do;
    %put WARNING: No blocks to export.;
    %return;
  %end;

  %if &by_chrom = 0 %then %do;
    %if %length(&file) = 0 %then %do;
      %put ERROR: by_chrom=0 needs file=.; %return;
    %end;
    data _null_;
      set &blocks;
      length _pmline $ 32767;
      %if &t = C %then %do; _pmline = &v; %end;
      %else %do;             _pmline = strip(put(&v, best32.)); %end;
      file "&file" lrecl=32767;
      put _pmline;
    run;
    %put NOTE: Wrote %pm_nobs(&blocks) IDs to &file..;
  %end;
  %else %do;
    %if %length(&dir) = 0 %then %do;
      %put ERROR: by_chrom=1 needs dir=.; %return;
    %end;
    %if %length(%pm_vtype(&blocks, CHROM)) = 0 %then %do;
      %put ERROR: by_chrom=1 requires a 'CHROM' column in &blocks..; %return;
    %end;
    data _null_;
      set &blocks;
      by CHROM notsorted;
      length _path $ 1024 _safe $ 64 _pmline $ 32767;
      /* the chromosome comes from the input file, so keep it out of the path */
      _safe = translate(cats(CHROM), '___', '/\.');
      _path = cats("&dir", "/snp_ch", _safe, ".txt");
      %if &t = C %then %do; _pmline = &v; %end;
      %else %do;             _pmline = strip(put(&v, best32.)); %end;
      file _pmout filevar=_path lrecl=32767;
      put _pmline;
    run;
    %put NOTE: Wrote per-chromosome SNP lists to &dir..;
  %end;
%mend pm_export;
