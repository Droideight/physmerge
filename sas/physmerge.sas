%macro pm_version;
  %put NOTE: physmerge for SAS 0.4.0 (sas-1), matching physmerge R 0.4.0.;
%mend pm_version;

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
      or upcase(name) = upcase(symget('pm__want'))
      or upcase(name) = upcase(prxchange('s/[^A-Za-z0-9_]/_/', -1,
      strip(symget('pm__want'))))
      or upcase(name) = upcase(cats('_', prxchange('s/[^A-Za-z0-9_]/_/', -1,
      strip(symget('pm__want'))))))
      order by case when upcase(label) = upcase(symget('pm__want')) then 1
      when upcase(name) = upcase(symget('pm__want')) then 2
      else 3 end;
  quit;
  &hit
    %mend pm_col;

%macro pm_read(path=, out=pm_sumstat, format=plink2,
  chrom_col=, pos_col=, id_col=, value_col=,
  test_filter=, test_col=TEST, test_val=ADD, chrom=,
  dlm='09'x, grows=MAX, quiet=0);

  %local fmt vc vp vv vi vt n0 n1 vvn;
  %let fmt = %upcase(&format);
  %let vvn = %sysfunc(getoption(validvarname));

  %if &fmt = PLINK2 %then %do;
    %if %length(&chrom_col) = 0 %then %let chrom_col = %str(#CHROM);
    %if %length(&pos_col) = 0 %then %let pos_col = POS;
    %if %length(&id_col) = 0 %then %let id_col = ID;
    %if %length(&value_col) = 0 %then %let value_col = P;
    %if %length(&test_filter) = 0 %then %let test_filter = 1;
  %end;
  %else %if &fmt = GPCM %then %do;
    %if %length(&chrom_col) = 0 %then %let chrom_col = %str(#CHROM);
    %if %length(&pos_col) = 0 %then %let pos_col = POS;
    %if %length(&id_col) = 0 %then %let id_col = ID;
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
    getnames = yes;
    grows = &grows;
  run;
  %if &syserr > 1 %then %do;
    %put ERROR: Failed to read file: &path;
    options validvarname=&vvn;
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
    options validvarname=&vvn;
    %return;
  %end;

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
        options validvarname=&vvn;
        %return;
      %end;
    %end;
  %end;

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

  options validvarname=&vvn;

  %global pm_chrom_var pm_pos_var pm_id_var pm_value_var pm_reward;
  %let pm_chrom_var = &vc;
  %let pm_pos_var = &vp;
  %let pm_id_var = &vi;
  %let pm_value_var = &vv;
  %if %index(%str( LOG10_P T_STAT Z_STAT CHISQ F_STAT T_STAT_DIRECT T_STAT_TE HPI ),
    %str( )%upcase(&value_col)%str( )) %then %let pm_reward = max;
  %else %let pm_reward = min;
  %if &quiet = 0 %then
    %put NOTE: physmerge: suggested reward = &pm_reward for value column &value_col..;

  proc datasets lib=work nolist nowarn; delete _pm_imp; quit;
%mend pm_read;

%macro pm_nobs(ds);
  %local dsid n rc;
  %let n = 0;
  %let dsid = %sysfunc(open(&ds));
  %if &dsid %then %do;
    %let n = %sysfunc(attrn(&dsid, nlobs));
    %let rc = %sysfunc(close(&dsid));
  %end;
  &n
    %mend pm_nobs;

%macro pm_prep(data=, out=_pm_prep, chrom=, pos=POSITION, value=VALUE, id=,
  idlen=200, chromlen=32);
  data &out(keep=_cord _chrom _pos _val _id _seq);
    length _chrom $ &chromlen _id $ &idlen _cord 8;
    if _n_ = 1 then do;
      declare hash _h();
      _h.defineKey('_chrom');
      _h.defineData('_chrord');
      _h.defineDone();
    end;
    retain _nc 0;
    set &data;
    %if %length(&chrom) %then %do; _chrom = cats(&chrom); %end;
    %else %do; _chrom = ''; %end;
    %if %length(&id) %then %do; _id = cats(&id); %end;
    %else %do; _id = ''; %end;
    _pos = &pos;
    _val = &value;

    if missing(_pos) or missing(_val) then delete;
    _seq + 1;
    if _h.find() ne 0 then do;
      _nc + 1;
      _cord = _nc;
      _h.add();
    end;
    drop _nc;
  run;
  proc sort data=&out; by _cord _pos _seq; run;
%mend pm_prep;

%macro physmerge(data=, out=pm_blocks, sig_th=5e-8, window=500000,
  reward=min, reset_on=any,
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

  data _pm_raw(keep=_cord _chrom start end rps_BP rps_value rps_ID);
    length _chrom $ &chromlen rps_ID $ &idlen;
    retain inblk steps best lpos start end rps_BP rps_value rps_ID;

    set _pm_prep;
    by _cord;

    if first._cord then do;
      inblk = 0;
      steps = &window;
      best = &sig_th;
      lpos = _pos;
    end;

    if not inblk then do;
      if _val &c &sig_th then link pm_open;
    end;
    else do;
      rem = steps - (_pos - lpos);
      if rem <= 0 then do;
        link pm_close;
        if _val &c &sig_th then link pm_open;
      end;
      else do;
        steps = rem;
        %if &ro = ANY %then %do;

          if _val &c &sig_th then do;
            steps = &window;
            if _val &c best then link pm_rep;
          end;
          else if _val &c best then do;
            steps = &window;
            link pm_rep;
          end;
        %end;
        %else %do;

          if _val &c best then do;
            steps = &window;
            link pm_rep;
          end;
        %end;
      end;
    end;

    lpos = _pos;
    if last._cord and inblk then link pm_close;
    return;

  pm_open:
    start = max(0, _pos - &window);
    end = .;
    rps_BP = _pos;
    rps_value = _val;
    rps_ID = _id;
    inblk = 1;
    steps = &window;
    best = _val;
    return;

  pm_rep:
    best = _val;
    rps_BP = _pos;
    rps_value = _val;
    rps_ID = _id;
    return;

  pm_close:
    end = lpos + steps;
    output;
    inblk = 0;
    steps = &window;
    best = &sig_th;
    return;
  run;

  data _pm_col(keep=_cord h_chrom h_start h_end h_bp h_val h_id);
    length h_chrom $ &chromlen h_id $ &idlen;
    retain held h_chrom h_start h_end h_bp h_val h_id;
    set _pm_raw(rename=(_chrom=i_chrom start=i_start end=i_end
      rps_BP=i_bp rps_value=i_val rps_ID=i_id));
    by _cord;

    if first._cord then held = 0;

    if held then do;
      if (i_bp - h_bp) < &window then do;
        if i_end > h_end then h_end = i_end;
        if i_val &c h_val then do;
          h_bp = i_bp; h_val = i_val; h_id = i_id;
        end;
      end;
      else do;
        if h_end > i_start then h_end = i_start;
        output;
        link pm_hold;
      end;
    end;
    else link pm_hold;

    if last._cord and held then output;
    return;

  pm_hold:
    h_chrom = i_chrom; h_start = i_start; h_end = i_end;
    h_bp = i_bp; h_val = i_val; h_id = i_id;
    held = 1;
    return;
  run;

  data &out(keep=serial CHROM start end rps_BP rps_ID rps_value);
    length serial 8 CHROM $ &chromlen start 8 end 8 rps_BP 8
      rps_ID $ &idlen rps_value 8;
    set _pm_col;
    serial + 1;
    CHROM = h_chrom;
    start = h_start;

    end = max(h_end, h_start);
    rps_BP = h_bp;
    rps_ID = h_id;
    rps_value = h_val;
  run;

  %let nout = %pm_nobs(&out);
  %if &quiet = 0 %then
    %put NOTE: physmerge: &nin SNPs -> &nout blocks (window=&window, sig_th=&sig_th, reward=&reward, reset_on=&reset_on).;

  proc datasets lib=work nolist nowarn; delete _pm_prep _pm_raw _pm_col; quit;
%mend physmerge;

%macro pm_annotate(blocks=pm_blocks, data=, out=pm_blocks_annot,
  chrom=, pos=POSITION, id=, keep=);

  %local addcols;

  proc sql noprint;
    select strip(name) into :addcols separated by ' '
      from dictionary.columns
      where libname = 'WORK' and memname = %upcase("&data")
      and upcase(name) not in ('POSITION', 'VALUE'
      %if %length(&chrom) %then , %upcase("&chrom") ;
            %if %length(&pos) %then , %upcase("&pos") ;
            %if %length(&id) %then , %upcase("&id") ; );
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
           %if %length(&id) %then and cats(d.&id) = b.rps_ID ;
     order by b.serial;
  quit;
%mend pm_annotate;

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
      %else %do; _pmline = strip(put(&v, best32.)); %end;
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

      _safe = translate(cats(CHROM), '___', '/\.');
      _path = cats("&dir", "/snp_ch", _safe, ".txt");
      %if &t = C %then %do; _pmline = &v; %end;
      %else %do; _pmline = strip(put(&v, best32.)); %end;
      file _pmout filevar=_path lrecl=32767;
      put _pmline;
    run;
    %put NOTE: Wrote per-chromosome SNP lists to &dir..;
  %end;
%mend pm_export;
