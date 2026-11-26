%macro pm_version;
  %put NOTE: physmerge for SAS 0.5.0.;
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
  %let hit = ;

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

%macro pm_kind(name);
  %local u;
  %let u = %qupcase(%superq(name));
  %if %sysfunc(prxmatch(/^(P|PVAL|PVALUE|P_VALUE|P-VALUE|P_.*)$/, &u)) %then p;
  %else %if %sysfunc(prxmatch(/^(NEG_LOG10_P|LOG10_P)$/, &u)) %then logp;
  %else other;
%mend pm_kind;

%macro pm_read(path=, out=pm_sumstat, format=plink2,
  chrom_col=, pos_col=, id_col=, value_col=,
  test_val=ADD, chrom=, sep=, quiet=0);

  %local fmt vc vp vv vi vt n0 n1 vvn dl sp tvs i q;
  %let fmt = %upcase(&format);
  %let vvn = %sysfunc(getoption(validvarname));

  %if &fmt = PLINK2 %then %do;
    %if %length(&chrom_col) = 0 %then %let chrom_col = %str(#CHROM);
    %if %length(&pos_col) = 0 %then %let pos_col = POS;
    %if %length(&id_col) = 0 %then %let id_col = ID;
    %if %length(&value_col) = 0 %then %let value_col = P;
  %end;
  %else %if &fmt = CUSTOM %then %do;
    %if %length(&chrom_col) = 0 or %length(&pos_col) = 0 or %length(&value_col) = 0 %then %do;
      %put ERROR: For format = 'custom', you must supply chrom_col, pos_col and value_col%str(;) set chrom_col = NA to ignore the chromosome column.;
      %return;
    %end;
  %end;
  %else %do;
    %put ERROR: format= must be plink2 or custom.;
    %return;
  %end;
  %if %length(&id_col) = 0 %then %let id_col = NA;
  %if %length(&chrom) and %upcase(&chrom_col) = NA %then %do;
    %put ERROR: chrom needs a chromosome column%str(;) it cannot be used with chrom_col = NA.;
    %return;
  %end;

  %if %upcase(%qscan(&path, -1, .)) = GZ %then %do;
    filename _pmsrc zip "&path" gzip;
  %end;
  %else %do;
    filename _pmsrc "&path";
  %end;

  %let sp = %lowcase(&sep);
  %if %length(&sp) = 0 %then %do;
    data _null_;
      infile _pmsrc obs=1 lrecl=1000000;
      input;
      length s $ 5;
      if index(_infile_, '09'x) then s = 'tab';
      else if index(_infile_, ',') then s = 'comma';
      else if index(_infile_, ';') then s = 'semi';
      else s = 'space';
      call symputx('pm__sep', s);
    run;
    %let sp = &pm__sep;
  %end;
  %if &sp = tab %then %let dl = '09'x;
  %else %if &sp = space %then %let dl = ' ';
  %else %if &sp = comma %then %let dl = ',';
  %else %if &sp = semi %then %let dl = ';';
  %else %if %length(&sep) = 1 %then %let dl = "&sep";
  %else %do;
    %put ERROR: sep must be tab, space or a single character.;
    %return;
  %end;

  options validvarname=v7;
  proc import datafile=_pmsrc out=_pm_imp dbms=dlm replace;
    delimiter = &dl;
    getnames = yes;
    guessingrows = max;
  run;
  %if &syserr > 1 %then %do;
    %put ERROR: Failed to read file: &path;
    options validvarname=&vvn;
    %return;
  %end;

  %if %upcase(&chrom_col) = NA %then %let vc = ;
  %else %let vc = %pm_col(_pm_imp, &chrom_col);
  %let vp = %pm_col(_pm_imp, &pos_col);
  %let vv = %pm_col(_pm_imp, &value_col);
  %let vt = %pm_col(_pm_imp, TEST);
  %if %upcase(&id_col) = NA %then %let vi = ;
  %else %let vi = %pm_col(_pm_imp, &id_col);

  %local miss;
  %let miss = ;
  %if %upcase(&chrom_col) ne NA and %length(&vc) = 0 %then %let miss = &miss &chrom_col;
  %if %length(&vp) = 0 %then %let miss = &miss &pos_col;
  %if %length(&vv) = 0 %then %let miss = &miss &value_col;
  %if %upcase(&id_col) ne NA and %length(&vi) = 0 %then %let miss = &miss &id_col;
  %if %length(&miss) %then %do;
    %put ERROR: Column(s) not found:&miss;
    options validvarname=&vvn;
    %return;
  %end;

  %if %pm_nobs(_pm_imp) = 0 %then %do;
    %put ERROR: No usable rows after filtering.;
    options validvarname=&vvn;
    %return;
  %end;

  %if &fmt = PLINK2 %then %do;
    %if %length(&vt) = 0 %then %do;
      %put WARNING: physmerge: TEST column not found%str(;) TEST filter skipped.;
    %end;
    %else %do;
      %let n0 = %pm_nobs(_pm_imp);
      data _pm_tst;
        set _pm_imp;
        where cats(&vt) = "&test_val";
      run;
      %let n1 = %pm_nobs(_pm_tst);
      %if &n1 = 0 %then %do;
        proc sql noprint outobs=10;
          select distinct cats(&vt) into :tvs separated by ', ' from _pm_imp
            where not missing(&vt);
        quit;
        %put ERROR: no rows have TEST = %str(%')&test_val%str(%')%str(;) TEST values in this file: %superq(tvs). Set test_val to one of them.;
        options validvarname=&vvn;
        %return;
      %end;
      data _pm_imp; set _pm_tst; run;
      %if &quiet = 0 %then
        %put NOTE: physmerge: TEST filter: kept &n1 of &n0 rows where TEST = %str(%')&test_val%str(%').;
    %end;
  %end;

  %if %length(&chrom) and %length(&vc) %then %do;
    %let q = ;
    %do i = 1 %to %sysfunc(countw(&chrom, %str( )));
      %let q = &q "%scan(&chrom, &i, %str( ))";
    %end;
    %let n0 = %pm_nobs(_pm_imp);
    data _pm_imp;
      set _pm_imp;
      if cats(&vc) in (&q);
    run;
    %let n1 = %pm_nobs(_pm_imp);
    %if &quiet = 0 and &n1 < &n0 %then
      %put NOTE: physmerge: %eval(&n0 - &n1) row(s) dropped by the chromosome filter.;
  %end;

  %let n0 = %pm_nobs(_pm_imp);
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
    if missing(POSITION) or missing(VALUE) then delete;
  run;
  %let n1 = %pm_nobs(&out);
  %if &quiet = 0 and &n1 < &n0 %then
    %put NOTE: physmerge: %eval(&n0 - &n1) row(s) dropped (NA in position or value).;
  options validvarname=&vvn;
  filename _pmsrc clear;
  %if &n1 = 0 %then %do;
    %put ERROR: No usable rows after filtering.;
    %return;
  %end;

  %global pm_out pm_chrom_var pm_pos_var pm_id_var pm_value_var pm_value_name;
  %let pm_out = %upcase(&out);
  %let pm_chrom_var = &vc;
  %let pm_pos_var = &vp;
  %let pm_id_var = &vi;
  %let pm_value_var = &vv;
  %let pm_value_name = %superq(value_col);

  proc datasets lib=work nolist nowarn; delete _pm_imp _pm_tst; quit;
%mend pm_read;

%macro pm_prep(data=, chrom=, pos=, value=, id=, sort=0, quiet=0);
  %global pm__cl pm__il pm__err pm__nin pm__neg pm__span pm__ndrop;
  %let pm__cl = 1;
  %let pm__il = 1;
  %let pm__err = ;
  %let pm__neg = 0;
  %let pm__span = 0;
  %let pm__ndrop = 0;
  proc sql noprint;
    select max(1, max(length(cats(&pos))))
      %if %length(&chrom) %then , max(1, max(length(cats(&chrom))));
      %if %length(&id) %then , max(1, max(length(cats(&id))));
      into :pm__dummy
      %if %length(&chrom) %then , :pm__cl;
      %if %length(&id) %then , :pm__il;
      from &data;
  quit;
  %let pm__cl = &pm__cl;
  %let pm__il = &pm__il;

  data _pm_prep(keep=_cord _chrom _pos _val _id _seq)
    _pm_chk(keep=_neg _min _max _ndrop);
    length _chrom $ &pm__cl _id $ &pm__il _cord 8;
    if _n_ = 1 then do;
      declare hash _h();
      _h.defineKey('_chrom');
      _h.defineData('_cord');
      _h.defineDone();
    end;
    retain _nc 0 _neg 0 _min . _max . _ndrop 0;
    set &data end=_eof;
    _seq = _n_;
    %if %length(&chrom) %then %do;
      if missing(&chrom) then _chrom = '';
      else _chrom = cats(&chrom);
    %end;
    %else %do; _chrom = ''; %end;
    %if %length(&id) %then %do; _id = cats(&id); %end;
    %else %do; _id = ''; %end;
    _pos = &pos;
    _val = &value;

    if missing(_pos) or missing(_val) then do;
      _ndrop + 1;
      if _eof then output _pm_chk;
      delete;
    end;
    if _pos < 0 then _neg = 1;
    _min = min(_min, _pos);
    _max = max(_max, _pos);
    if _h.find() ne 0 then do;
      _nc + 1;
      _cord = _nc;
      _h.add();
    end;
    output _pm_prep;
    if _eof then output _pm_chk;
    drop _nc;
  run;

  data _null_;
    set _pm_chk;
    call symputx('pm__neg', _neg);
    call symputx('pm__span', _max - _min);
    call symputx('pm__ndrop', _ndrop);
  run;
  %if &quiet = 0 and &pm__ndrop > 0 %then
    %put NOTE: physmerge: &pm__ndrop row(s) dropped (NA in position or value).;
  %if &pm__neg = 1 %then
    %put WARNING: physmerge: negative position(s) found%str(;) block boundaries are clamped at 0.;
  %if %length(&chrom) = 0 and %sysevalf(&pm__span > 2.5e8) %then
    %put WARNING: physmerge: positions span more than 250 Mb but there is no chromosome column%str(;) blocks may cross chromosome boundaries.;

  %if &sort = 1 %then %do;
    proc sort data=_pm_prep; by _cord _pos _seq; run;
  %end;
  %else %do;
    data _null_;
      if _n_ = 1 then do;
        declare hash _s();
        _s.defineKey('_chrom');
        _s.defineDone();
      end;
      length _pc $ &pm__cl _msg $ 512;
      retain _pc '' _pp . _first 1;
      set _pm_prep;
      if not _first and _chrom = _pc then do;
        if _pos < _pp then do;
          %if %length(&chrom) %then %do;
            _msg = 'input is not position-sorted (chromosome ' || strip(_chrom) || ': ' ||
              strip(put(_pos, 32.)) || ' after ' || strip(put(_pp, 32.)) ||
              '); set sort to merge it anyway.';
          %end;
          %else %do;
            _msg = 'input is not position-sorted (' || strip(put(_pos, 32.)) || ' after ' ||
              strip(put(_pp, 32.)) || '); set sort to merge it anyway.';
          %end;
          call symputx('pm__err', _msg);
          stop;
        end;
      end;
      else do;
        if _s.check() = 0 then do;
          _msg = catx(' ', 'chromosome', _chrom,
            'appears in more than one stretch of the input; set sort to merge it anyway.');
          call symputx('pm__err', _msg);
          stop;
        end;
        _s.add();
      end;
      _pc = _chrom;
      _pp = _pos;
      _first = 0;
    run;
  %end;
  %let pm__nin = %pm_nobs(_pm_prep);
%mend pm_prep;

%macro physmerge(data=, out=pm_blocks, sig_th=, window=500000,
  reward=, chrom=, pos=POSITION, value=VALUE, id=,
  sort=0, quiet=0);

  %local c rw nout vname kind vout fromread;
  %let fromread = 0;
  %if %symexist(pm_out) %then %do;
    %if %upcase(&data) = &pm_out %then %let fromread = 1;
  %end;
  %if &fromread %then %do;
    %if %length(&chrom) = 0 %then %let chrom = &pm_chrom_var;
    %if %length(&id) = 0 %then %let id = &pm_id_var;
    %let vname = %superq(pm_value_name);
    %let vout = &pm_value_var;
  %end;
  %else %do;
    %if %length(&chrom) = 0 and %length(%pm_vtype(&data, CHROM)) %then %let chrom = CHROM;
    %if %length(&id) = 0 and %length(%pm_vtype(&data, ID)) %then %let id = ID;
    %let vname = &value;
    %let vout = &value;
  %end;
  %let vout = %substr(rps_&vout, 1, %sysfunc(min(32, %length(rps_&vout))));

  %let kind = %pm_kind(&vname);
  %if (%length(&reward) = 0 or %length(&sig_th) = 0) and &kind = other %then %do;
    %put ERROR: cannot tell whether a larger or a smaller %superq(vname) is more significant%str(;) set reward and sig_th.;
    %return;
  %end;
  %if %length(&reward) = 0 %then %do;
    %if &kind = p %then %let reward = min;
    %else %let reward = max;
  %end;
  %if %length(&sig_th) = 0 %then %do;
    %if &kind = p %then %let sig_th = 5e-8;
    %else %let sig_th = %sysevalf(-%sysfunc(log10(5e-8)));
  %end;

  %let rw = %upcase(&reward);
  %if &rw = MAX %then %let c = >;
  %else %if &rw = MIN %then %let c = <;
  %else %do; %put ERROR: `reward` must be either 'min' or 'max'.; %return; %end;
  %if %sysevalf(&window <= 0) %then %do;
    %put ERROR: `window` must be a single positive numeric value.; %return;
  %end;

  %pm_prep(data=&data, chrom=&chrom, pos=&pos, value=&value, id=&id,
    sort=&sort, quiet=&quiet);
  %if %length(&pm__err) %then %do;
    %put ERROR: %superq(pm__err);
    %return;
  %end;

  data _pm_raw(keep=_cord _chrom start end rps_BP rps_value rps_ID rps_seq);
    length _chrom $ &pm__cl rps_ID $ &pm__il;
    retain inblk steps best lpos start end rps_BP rps_value rps_ID rps_seq;

    set _pm_prep;
    by _cord;

    if first._cord then do;
      inblk = 0;
      steps = &window;
      best = &sig_th;
      lpos = _pos;
    end;

    if inblk then do;
      rem = steps - (_pos - lpos);
      if rem <= 0 then link pm_close;
      else do;
        steps = rem;
        if _val &c &sig_th then do;
          steps = &window;
          if _val &c best then link pm_rep;
        end;
      end;
    end;
    if not inblk and _val &c &sig_th then link pm_open;

    lpos = _pos;
    if last._cord and inblk then link pm_close;
    return;

  pm_open:
    start = max(0, _pos - &window);
    end = .;
    rps_BP = _pos;
    rps_value = _val;
    rps_ID = _id;
    rps_seq = _seq;
    inblk = 1;
    steps = &window;
    best = _val;
    return;

  pm_rep:
    best = _val;
    rps_BP = _pos;
    rps_value = _val;
    rps_ID = _id;
    rps_seq = _seq;
    return;

  pm_close:
    end = lpos + steps;
    output;
    inblk = 0;
    steps = &window;
    best = &sig_th;
    return;
  run;

  data _pm_col(keep=_cord h_chrom h_start h_end h_bp h_val h_id h_seq);
    length h_chrom $ &pm__cl h_id $ &pm__il;
    retain held h_chrom h_start h_end h_bp h_val h_id h_seq;
    set _pm_raw(rename=(_chrom=i_chrom start=i_start end=i_end
      rps_BP=i_bp rps_value=i_val rps_ID=i_id rps_seq=i_seq));
    by _cord;

    if first._cord then held = 0;

    if held then do;
      if (i_bp - h_bp) < &window then do;
        if i_end > h_end then h_end = i_end;
        if i_val &c h_val then do;
          h_bp = i_bp; h_val = i_val; h_id = i_id; h_seq = i_seq;
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
    h_bp = i_bp; h_val = i_val; h_id = i_id; h_seq = i_seq;
    held = 1;
    return;
  run;

  data &out(keep=serial %if %length(&chrom) %then CHROM; start end rps_BP
      %if %length(&id) %then rps_ID; &vout)
    _pm_rows(keep=serial _seq);
    length serial 8 %if %length(&chrom) %then CHROM $ &pm__cl; start 8 end 8 rps_BP 8
      %if %length(&id) %then rps_ID $ &pm__il; &vout 8;
    set _pm_col;
    serial + 1;
    %if %length(&chrom) %then CHROM = h_chrom;;
    start = h_start;
    end = max(h_end, h_start);
    rps_BP = h_bp;
    %if %length(&id) %then rps_ID = h_id;;
    &vout = h_val;
    _seq = h_seq;
  run;

  %global pm_blocks_out pm_blocks_data;
  %let pm_blocks_out = %upcase(&out);
  %let pm_blocks_data = %upcase(&data);

  %let nout = %pm_nobs(&out);
  %if &quiet = 0 %then %do;
    data _null_;
      length s $ 40;
      x = &sig_th;
      if x = 0 then s = '0';
      else do;
        e = floor(log10(abs(x)));
        if e < -4 or e >= 6 then do;
          m = round(x / 10**e, 1e-5);
          if abs(m) >= 10 then do; m = m / 10; e = e + 1; end;
          s = cats(put(m, best12.), 'e', ifc(e < 0, '-', '+'), put(abs(e), z2.));
        end;
        else s = strip(put(round(x, 10**(e - 5)), best12.));
      end;
      call symputx('pm__sg', s);
      call symputx('pm__wn', strip(put(&window, best15.)));
    run;
    %put NOTE: physmerge: &pm__nin SNPs -> &nout blocks (window=&pm__wn, sig_th=&pm__sg, reward=%lowcase(&reward)).;
  %end;

  proc datasets lib=work nolist nowarn; delete _pm_prep _pm_chk _pm_raw _pm_col; quit;
%mend physmerge;

%macro pm_annotate(blocks=pm_blocks, data=, out=pm_blocks_annot);

  %local addcols roles i cv exact;
  %let roles = 'POSITION' 'VALUE' 'CHROM' 'ID';
  %let cv = ;
  %if %length(%pm_vtype(&data, CHROM)) %then %let cv = CHROM;
  %if %symexist(pm_out) %then %do;
    %if %upcase(&data) = %superq(pm_out) %then %do;
      %let roles = 'POSITION' 'VALUE' "%upcase(&pm_chrom_var)" "%upcase(&pm_pos_var)"
        "%upcase(&pm_id_var)" "%upcase(&pm_value_var)";
      %let cv = &pm_chrom_var;
    %end;
  %end;
  %let exact = 0;
  %if %symexist(pm_blocks_out) %then %do;
    %if %upcase(&blocks) = %superq(pm_blocks_out) and %upcase(&data) = %superq(pm_blocks_data) %then
      %let exact = 1;
  %end;

  proc sql noprint;
    select strip(name) into :addcols separated by ' '
      from dictionary.columns
      where libname = 'WORK' and memname = %upcase("&data")
      and upcase(name) not in (&roles)
      order by varnum;
  quit;

  data _pm_dat;
    set &data;
    _seq = _n_;
  run;

  %if &exact %then %do;
    proc sql;
      create table &out as
        select b.*
        %if %length(&addcols) %then %do i = 1 %to %sysfunc(countw(&addcols, %str( )));
          , d.%scan(&addcols, &i, %str( ))
        %end;
        from &blocks as b
          left join _pm_rows as r on r.serial = b.serial
          left join _pm_dat as d on d._seq = r._seq
        order by b.serial;
    quit;
  %end;
  %else %do;
    %put WARNING: physmerge: &blocks is not the last output of physmerge on &data%str(;) matching blocks by position.;
    proc sort data=_pm_dat nodupkey; by %if %length(&cv) %then &cv; POSITION; run;
    proc sql;
      create table &out as
        select b.*
        %if %length(&addcols) %then %do i = 1 %to %sysfunc(countw(&addcols, %str( )));
          , d.%scan(&addcols, &i, %str( ))
        %end;
        from &blocks as b
          left join _pm_dat as d on d.POSITION = b.rps_BP
          %if %length(&cv) %then and cats(d.&cv) = b.CHROM;
        order by b.serial;
    quit;
  %end;
  proc datasets lib=work nolist nowarn; delete _pm_dat; quit;
%mend pm_annotate;

%macro pm_mkdir(dir);
  %local i part path rc;
  %let path = ;
  %if %qsubstr(&dir, 1, 1) = / %then %let path = /;
  %do i = 1 %to %sysfunc(countw(&dir, /\));
    %let part = %scan(&dir, &i, /\);
    %if %length(&path) and %qsubstr(&path, %length(&path), 1) ne / %then %let path = &path/;
    %let path = &path&part;
    %if %sysfunc(fileexist(&path)) = 0 %then %do;
      %if %length(&path) = %length(&part) %then %let rc = %sysfunc(dcreate(&part, .));
      %else %let rc = %sysfunc(dcreate(&part, %qsubstr(&path, 1, %eval(%length(&path) - %length(&part)))));
    %end;
  %end;
%mend pm_mkdir;

%macro pm_export(blocks=pm_blocks, file=, dir=, quiet=0);
  %local v t n nc ns;
  %if (%length(&file) = 0) = (%length(&dir) = 0) %then %do;
    %put ERROR: give either file= or dir=, not both.;
    %return;
  %end;
  %if %length(%pm_vtype(&blocks, rps_ID)) %then %let v = rps_ID;
  %else %let v = rps_BP;
  %let t = %pm_vtype(&blocks, &v);
  %let n = %pm_nobs(&blocks);
  %if &n = 0 %then %put WARNING: physmerge: no blocks to export.;

  %if %length(&file) %then %do;
    data _null_;
      file "&file" lrecl=32767;
      do while (not _eof);
        set &blocks end=_eof;
        length _pmline $ 32767;
        %if &t = C %then %do; _pmline = &v; %end;
        %else %do; _pmline = strip(put(&v, best32.)); %end;
        put _pmline;
      end;
      stop;
    run;
    %if &quiet = 0 %then %put NOTE: physmerge: wrote &n IDs to &file..;
  %end;
  %else %do;
    %if %length(%pm_vtype(&blocks, CHROM)) = 0 %then %do;
      %put ERROR: a per-chromosome list needs a chromosome column.; %return;
    %end;
    %pm_mkdir(&dir);
    %if &n > 0 %then %do;
      proc sql noprint;
        select count(distinct CHROM), count(distinct translate(cats(CHROM), '___', '/\.'))
          into :nc trimmed, :ns trimmed from &blocks;
      quit;
      %if &nc ne &ns %then %do;
        %put ERROR: chromosome names collide once made safe for a file name.; %return;
      %end;
      data _null_;
        set &blocks;
        length _path $ 1024 _safe $ 256 _pmline $ 32767;
        _safe = translate(cats(CHROM), '___', '/\.');
        _path = cats("&dir", "/snp_ch", _safe, ".txt");
        %if &t = C %then %do; _pmline = &v; %end;
        %else %do; _pmline = strip(put(&v, best32.)); %end;
        file _pmout filevar=_path lrecl=32767;
        put _pmline;
      run;
    %end;
    %else %let nc = 0;
    %if &quiet = 0 %then %put NOTE: physmerge: wrote &nc chromosome file(s) to &dir..;
  %end;
%mend pm_export;
