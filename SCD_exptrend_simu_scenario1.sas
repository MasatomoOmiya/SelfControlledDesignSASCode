/* Note: The following SAS program is used to generate data for each scenario 
in order to investigate the impact of exposure time trends in a self-controlled case series design through a simulation study, 
and to execute the macro program "%sim_exptimetrend". 
Before running this program, please download "sccs.sas", "element.sas", and "poisreg.sas" from the Open University website
in the UK at http://statistics.open.ac.uk/sccs/index.htm, save them on your computer, 
and specify the location where the programs are saved in &macrofile.*/

/*Scenario 1*/
%macro sim_exptimetrend_sc1(replis=, replin=, Nday=, Nsample=);
proc datasets;
	delete result_sccs result_cco result_sbi result_ctc result_cctc;
quit;
%do repli= &replis %to &replin;
%let seed1 = (&repli*12345);
%let seed2 = (&repli*54321);
%let seed3 = (&repli*4321);

data simd01;
call streaminit(&seed1);
do ID = 1 to &Nsample;
  retain cov1 cov2;
  cov1 = rand("Bernoulli", 0.05);
  cov2 = rand("Bernoulli", 0.5);
  do day = 1 to &Nday;
  output;
  end;
end;
run;

data simd02 exp01(keep=id day exp);
call streaminit(&seed1);
  set simd01;
  %if &trend=U %then %do;
  exp = rand ("Bernoulli", exp(&theta0 + &theta1*cov1 + &theta2*cov2)*(1 +((&thetat.-1)/364)*(day-1)));
  output simd02;
  if exp=1 then output exp01;
  %end;
  %else %if &trend=C %then %do;
  pai=constant("pi");
  exp = rand ("Bernoulli", exp(&theta0 + &theta1*cov1 + &theta2*cov2)*(1+0.4*cos(2*pai*(day/365))));
  output simd02;
  if exp=1 then output exp01;
  %end;
  %else %if &trend=M %then %do;
  if  0<=day<=90 then do;
  exp = rand ("Bernoulli", 0);end;
  else if  90<day<=180 then do;
  exp = rand ("Bernoulli", (exp(&thetam.)/(180-90))*day + exp(&thetam.) - (((exp(&thetam.))*180)/(180-90))); end;
  else if  180<day<=270 then do;
  exp = rand ("Bernoulli", (exp(&thetam.)/(180-270))*day-(exp(&thetam.)*270/(180-270)));end;
  else if  270<day<=365 then do;
  exp = rand ("Bernoulli", 0);end;
  output simd02;
  if exp=1 then output exp01;
  %end;
run;

data simd03;
  set simd02;
  retain exped;
  by ID;
  if first.id then exped=0;
  if exp=1 then exped=min(day+(&exprisk. - 1), &Nday);
run;

data simd04(drop=exped);
  set simd03;
  if day =< exped then exp=1;
run;

data event01(keep=id day event);
call streaminit(&seed2);
  set simd04;
  event = rand ("Bernoulli", exp(&beta0 + &beta1*exp + &beta2*cov1 + &beta3*cov2));
  if event=1 then output;
run;

data exp02;
  set exp01;
  retain numexp;
  by id;
  if first.id then numexp=0;
  numexp=numexp+1;
run;

proc sql;
  select max(numexp) into :maxexp2
  from exp02;
quit;

proc transpose data=exp02 out=exp03(drop=_NAME_) prefix=exp;
  var day;
  by id;
run;

data cens01;
 call streaminit(&seed2);
 set event01;
 censor = rand("Bernoulli", &pcens);
run;

data cens01x(keep=id cday);
 set cens01;
 if censor=1;
 cday=day;
run;

proc sort data=cens01x;
 by id cday;
run;

data cens01xx;
  set cens01x;
  by id;
  if first.id;
run;

data event02;
  set event01;
  retain numevt;
  by id;
  if first.id then numevt=0;
  numevt=numevt+1;
run;

proc sql;
  select max(numevt) into :maxevent
  from event02;
quit;

proc transpose data=event02 out=event03(drop=_NAME_) prefix=event;
  var day;
  by id;
run;

data ctc_contid(keep=id);
  set simd04;
run;

proc sort data=ctc_contid nodupkey;
  by id;
run;

data noevent(keep=id day event);
  merge ctc_contid event01;
  by id;
  if event=. then event=0;
  if event=1 then delete;
run;

proc transpose data=noevent out=noevent01(drop=_NAME_) prefix=event;
  var day;
  by id;
run;

data noevent02;
  set noevent01;
  if event1=. then event1=0;
run;

data ctc_cont01;
  merge exp03 noevent02;
  by id;
  if event1=. then delete;
run;

data simd06;
	merge exp03 event03;
	by id;
run;

proc sql;
create table simd06x as
select a.*, b.cday
from simd06 a left join cens01xx b on a.ID=b.ID;
quit;

data simd07;
  set simd06x;
  startobs=1;
  endobs= min(&Nday, cday);
  dob=0;
  if event1=. then delete;
run;

%let agerange=1 &Nday;
%let risk=0 (&riskp-1);
%let age=;
%let season=;
%let semi=;

%sccs(data=simd07, dob_raw=dob, pid=id, events= %do i=1 %to &maxevent; event&i  %end;, vacc=%do i=1 %to &maxexp2; exp&i  %end;, startst=startobs, endst=endobs, outdata=sccs01);
%poisreg(data=sccs01,y=nevt,covar=int riskr1,offset=offset,elim=ID,prntyn=N,title="SCD_exptrend");

data out2(keep= estimate sterr ll ul expest expll expul);
  set out;
  where params like '%RISK%';
run;

proc append base=result_sccs data=out2;
run;

/*****************************************************************************************************************************
CCO_1:2
*****************************************************************************************************************************/

data simd38_C(keep=id %do i=1 %to &maxexp2; exp&i  %end; event1 contst conted contst2 conted2 casest caseed end);
  set simd07;
  if event1 < (3*&riskp+&lag) or endobs < event1 then delete;
  end=min(endobs, event1);
  contst=event1 - ((3*&riskp + &lag)-1);
  conted=event1 - (2*&riskp + &lag);
  contst2=event1 - ((2*&riskp+&lag)-1);
  conted2=event1 - (&riskp + &lag);
  casest=event1 - (&riskp-1);
  caseed=event1;
run;

data simd39_C;
  set simd38_C;
  %do i=1 %to &maxexp2;
  cont&i =(contst <= exp&i and exp&i <= conted);
  %end;
  %do i=1 %to &maxexp2;
  cont2&i =(contst2 <= exp&i and exp&i <= conted2);
  %end;
  %do j=1 %to &maxexp2;
  case&j = (casest <= exp&j and exp&j <= caseed);
  %end;
  cont=0
  %do i=1 %to &maxexp2;
  + cont&i
  %end;;
  cont_2=0
  %do i=1 %to &maxexp2;
  + cont2&i
  %end;;
  case=0
  %do j=1 %to &maxexp2;
  + case&j
  %end;;
run;

data simd40_C(keep=id exp);
  set simd39_C;
  exp=1;
  if case => 1 then output;
run;
proc sort data=simd40_C nodupkey;
  by ID;
run;
proc sql;
  create table simd41_C as
  select a.id, b.exp
  from simd38_C a left join simd40_C b on a.ID=b.ID;
quit; 

data simd42_C;
  set simd41_C;
  period=1;
  if exp=. then exp=0;
run;

data simd43_C(keep=id exp);
  set simd39_C;
  exp=1;
  if cont => 1 then output;
run;
proc sort data=simd43_C nodupkey;
  by ID;
run;

proc sql;
  create table simd44_C as
  select a.id, b.exp 
  from simd38_C a left join simd43_C b on a.ID=b.ID;
quit;  

data simd45_C;
  set simd44_C;
  period=0;
  if exp=. then exp=0;
run;

data simd46_C(keep=id exp);
  set simd39_C;
  exp=1;
  if cont_2 => 1 then output;
run;

proc sort data=simd46_C nodupkey;
  by ID;
run;

proc sql;
  create table simd47_C as
  select a.id, b.exp 
  from simd38_C a left join simd46_C b on a.ID=b.ID;
quit;  

data simd48_C;
  set simd47_C;
  period=0;
  if exp=. then exp=0;
run;

proc sort data=simd48_C nodupkey;
  by ID;
run;

data cco01phreg_C;
  set simd42_C simd45_C simd48_C ;
  time=2-period;
run; 

proc sort data=cco01phreg_C;
  by id time;
run;

proc phreg data=cco01phreg_C outest=outcco_C covout ;
  model time*period(0)=exp/ ties=Breslow;
  strata id;
run;

proc transpose data=outcco_C out=outcco2_Cpara (drop=_NAME_) prefix=param;
  var exp;
run;

proc transpose data=outcco_C out=outcco2_Cconv (drop=_NAME_ _LABEL_) prefix=conv;
  var _status_;
run;

data outcco2_C;
  merge outcco2_Cpara outcco2_Cconv;
run;

data outcco3_C (keep=estimate sterr ll ul expest expll expul conv1 conv2);
  set outcco2_C;
  estimate=param1;
  sterr=sqrt(param2);
  ll=param1-1.96*sqrt(param2);
  ul=param1+1.96*sqrt(param2);
  expest=exp(estimate);
  expll=exp(ll);
  expul=exp(ul);
run;

proc append base=result_cco data=outcco3_C;
run;

/*****************************************************************************************************************************
SBI
*****************************************************************************************************************************/

data simd38_BC(keep=id %do i=1 %to &maxexp2; exp&i  %end; event1 contst conted contst2 conted2 casest caseed end);
  set simd07;
  if event1-(2*&riskp+&lag) < 0 or endobs < (event1 + &lag + &riskp) then delete;
  end=min(endobs, event1 + &riskp + &lag);
  contst=event1 - ((2*&riskp+&lag)-1);
  conted=event1 - (&riskp + &lag);
  contst2=event1 + &lag+1;
  conted2=event1 + (&riskp + &lag);
  casest=event1 - (&riskp-1);
  caseed=event1;
run;

data simd39_BC;
  set simd38_BC;
  %do i=1 %to &maxexp2;
  cont&i =(contst <= exp&i and exp&i <= conted);
  %end;
  %do i=1 %to &maxexp2;
  cont2&i =(contst2 <= exp&i and exp&i <= conted2);
  %end;
  %do j=1 %to &maxexp2;
  case&j = (casest <= exp&j and exp&j <= caseed);
  %end;
  cont=0
  %do i=1 %to &maxexp2;
  + cont&i
  %end;;
  cont_2=0
  %do i=1 %to &maxexp2;
  + cont2&i
  %end;;
  case=0
  %do j=1 %to &maxexp2;
  + case&j
  %end;;
run;

data simd40_BC(keep=id exp);
  set simd39_BC;
  exp=1;
  if case => 1 then output;
run;

proc sort data=simd40_BC nodupkey;
  by ID;
run;

proc sql;
  create table simd41_BC as
  select a.id, b.exp
  from simd38_BC a left join simd40_BC b on a.ID=b.ID;
quit; 

data simd42_BC;
  set simd41_BC;
  period=1;
  if exp=. then exp=0;
run;

data simd43_BC(keep=id exp);
  set simd39_BC;
  exp=1;
  if cont => 1 then output;
run;

proc sort data=simd43_BC nodupkey;
  by ID;
run;

proc sql;
  create table simd44_BC as
  select a.id, b.exp 
  from simd38_BC a left join simd43_BC b on a.ID=b.ID;
quit;  

data simd45_BC;
  set simd44_BC;
  period=0;
  if exp=. then exp=0;
run;

data simd46_BC(keep=id exp);
  set simd39_BC;
  exp=1;
  if cont_2 => 1 then output;
run;

proc sort data=simd46_BC nodupkey;
  by ID;
run;

proc sql;
  create table simd47_BC as
  select a.id, b.exp 
  from simd38_BC a left join simd46_BC b on a.ID=b.ID;
quit;  

data simd48_BC;
  set simd47_BC;
  period=0;
  if exp=. then exp=0;
run;

proc sort data=simd48_BC nodupkey;
  by ID;
run;

data cco01phreg_BC;
  set simd42_BC simd45_BC simd48_BC ;
  time=2-period;
run; 

proc sort data=cco01phreg_BC ;
  by id time;
run;

proc phreg data=cco01phreg_BC outest=outcco_BC covout ;
  model time*period(0)=exp/ ties=Breslow;
  strata id;
run;

proc transpose data=outcco_BC out=outcco2_BCpara (drop=_NAME_) prefix=param;
  var exp;
run;

proc transpose data=outcco_BC out=outcco2_BCconv (drop=_NAME_ _LABEL_) prefix=conv;
  var _status_;
run;

data outcco2_BC;
 merge outcco2_BCpara outcco2_BCconv;
run;

data outcco3_BC (keep=estimate sterr ll ul expest expll expul conv1 conv2);
  set outcco2_BC;
  estimate=param1;
  sterr=sqrt(param2);
  ll=param1-1.96*sqrt(param2);
  ul=param1+1.96*sqrt(param2);
  expest=exp(estimate);
  expll=exp(ll);
  expul=exp(ul);
run;

proc append base=result_SBI data=outcco3_BC;
run;

/*****************************************************************************************************************************
CTC
*****************************************************************************************************************************/

data ctc_control;
  set ctc_cont01;
  startobs=1;
  endobs= min(&Nday, cday);
  if event1=0 then group=0;
  else group=1;
run;

data ctc_case;
  set simd07;
  group=1;
run;

proc sql;
  create table ctc_data01
  as select b.group, a.id , b.id as cont_id, %do i=1 %to &maxexp2; b.exp&i, %end; b.startobs, b.endobs, a.event1, ranuni(&seed3.) as rand_num
  from ctc_case a, ctc_control b
  having 0 <= rand_num <= 0.01; 
;
quit;

proc sort data=ctc_data01;
  by id rand_num;
run;

data ctc_data02;
  set ctc_data01;
  by id rand_num;
  retain mcount;
  if first.id then mcount=0;
  mcount+1;
  if mcount le 2;
run;

proc sql;
  create table ctc_data_all as
  select group, id ,%do i=1 %to &maxexp2; exp&i, %end; startobs, endobs, event1 from ctc_case
  union all
  select group, id ,%do i=1 %to &maxexp2; exp&i, %end; startobs, endobs, event1 from ctc_data02
  order by id, group desc;
quit;

data ctc_data_all2;
  set ctc_data_all;
  by id;
  retain no;
  if first.id then no=0;
  no+1;
  output ctc_data_all2;
run;

data not_ctc;
  set ctc_data_all2;
  by id;
  if last.id then do;
  if no le 2 then output not_ctc; 
  end;
run;

proc sort data=ctc_data_all2; by id;run;
proc sort data=not_ctc;by id;run;

data simd37_ctc;
  merge ctc_data_all2 not_ctc(in=b_);
  by id;
  if b_ then delete;
run;

proc sort data=simd37_ctc;
  by id descending group;
run;

data simd38_ctc(keep=id %do i=1 %to &maxexp2; exp&i  %end; event1 contst conted contst2 conted2 casest caseed end group);
  set simd37_ctc;
  if event1 < (3*&riskp+&lag) or endobs < event1 then delete;
  end=min(endobs, event1);
  contst=event1 - ((3*&riskp + &lag)-1);
  conted=event1 - (2*&riskp + &lag);
  contst2=event1 - ((2*&riskp+&lag)-1);
  conted2=event1 - (&riskp + &lag);
  casest=event1 - (&riskp-1);
  caseed=event1;
run;

data simd38_ctc;
  set simd38_ctc;
  by id;
  if first.id then seq=0;
  seq+1;
run;

data simd39_ctc;
  set simd38_ctc;
  %do i=1 %to &maxexp2;
  cont&i =(contst <= exp&i and exp&i <= conted);
  %end;
  %do i=1 %to &maxexp2;
  cont2&i =(contst2 <= exp&i and exp&i <= conted2);
  %end;
  %do j=1 %to &maxexp2;
  case&j = (casest <= exp&j and exp&j <= caseed);
  %end;
  cont=0
  %do i=1 %to &maxexp2;
  + cont&i
  %end;;
  cont_2=0
  %do i=1 %to &maxexp2;
  + cont2&i
  %end;;
  case=0
  %do j=1 %to &maxexp2;
  + case&j
  %end;;
run;

data simd40_ctc(keep=id exp group seq);
  set simd39_ctc;
  exp=1;
  if case => 1 then output;
run;

proc sort data=simd40_ctc;
  by ID seq;
run;

proc sql;
  create table simd41_ctc as
  select a.id, b.exp, a.group, a.seq
  from simd38_ctc a left join simd40_ctc b on a.ID=b.ID and a.seq=b.seq;
quit; 

data simd42_ctc;
  set simd41_ctc;
  period=1;
  ctcid=1;
  if exp=. then exp=0;
run;

data simd43_ctc(keep=id exp group seq);
  set simd39_ctc;
  exp=1;
  if cont => 1 then output;
run;

proc sort data=simd43_ctc;
	by ID seq;
run;

proc sql;
  create table simd44_ctc as
  select a.id, b.exp, a.group, a.seq
  from simd38_ctc a left join simd43_ctc b on a.ID=b.ID and a.seq=b.seq;
quit;  

data simd45_ctc;
  set simd44_ctc;
  period=0;
  ctcid=2;
  if exp=. then exp=0;
run;

data simd46_ctc(keep=id exp group seq);
  set simd39_ctc;
  exp=1;
  if cont_2 => 1 then output;
run;

proc sort data=simd46_ctc;
  by ID seq;
run;

proc sql;
  create table simd47_ctc as
  select a.id, b.exp, a.group, a.seq
  from simd38_ctc a left join simd46_ctc b on a.ID=b.ID and a.seq=b.seq;
quit;  
data simd48_ctc;
  set simd47_ctc;
  period=0;
  ctcid=3;
  if exp=. then exp=0;
run;

proc sort data=simd48_ctc ;
  by ID seq;
run;

data cco01logistic_ctc;
  set simd42_ctc simd45_ctc simd48_ctc ;
  time=2-period;
run; 

proc sort data=cco01logistic_ctc ;
  by id seq ctcid;
run;

data cco01logistic_ctc;
  set cco01logistic_ctc;
  count=_n_;
  ctc_id=ceil(count/3);
run;

proc logistic data=cco01logistic_ctc outest=outcco_ctc covout descending;
  strata ctc_id;
  model period=exp exp*group;
run;

proc transpose data=outcco_ctc out=outcco2_ctc0 (drop=_NAME_) prefix=param;
  var exp expgroup;
run;

proc transpose data=outcco_ctc out=outcco2_ctcconv (drop=_NAME_ _LABEL_) prefix=conv;
  var _status_;
run;

data outcco2_ctc0;
  set outcco2_ctc0;
  id=_N_;
run;

data outcco2_ctc1(drop=_label_) outcco2_ctc2(drop=_label_);
  set outcco2_ctc0;
  if id=1 then output outcco2_ctc1;
  if id=2 then output outcco2_ctc2;
run;

data outcco2_ctc1(drop=id param3);
  set outcco2_ctc1;
  kekka='timetrend';
run;

data outcco2_ctc2(drop=id param2);
  set outcco2_ctc2;
  kekka='ctc';
  rename param3=param2;
run;

data outcco2_ctcpara;
  set outcco2_ctc1 outcco2_ctc2;
run;

proc sql;
  create table outcco2_ctc
  as select *
  from outcco2_ctcpara, outcco2_ctcconv;
quit;

data outcco3_ctc (keep=kekka estimate sterr ll ul expest expll expul conv1 conv2 conv3);
  set outcco2_ctc;
  estimate=param1;
  sterr=sqrt(param2);
  ll=param1-1.96*sqrt(param2);
  ul=param1+1.96*sqrt(param2);
  expest=exp(estimate);
  expll=exp(ll);
  expul=exp(ul);
run;

proc append base=result_ctc data=outcco3_ctc;
run;

/*****************************************************************************************************************************
CCTC
*****************************************************************************************************************************/
data cactc_case cactc_cont;
  set simd07;
run;

data cactc_caseid(keep=id event1 group) cactc_case;
  set cactc_case;
  group=1;
run;

data cactc_cont;
  set cactc_cont;
  group=0;
run;

proc sql;
  create table cactc_data01
  as select b.group, a.id , b.id as cont_id, %do i=1 %to &maxexp2; exp&i, %end; b.startobs, b.endobs, b.event1 as refevent1, a.event1, 
  ranuni(&seed3.) as rand_num 
  from cactc_caseid a, cactc_cont b
  where (a.event1 < b.event1 - &cactcday.)
  having 0 <= rand_num <= 0.1
;
quit;

proc sort data=cactc_data01;
  by id rand_num;
run;

data cactc_data02;
  set cactc_data01;
  by id rand_num;
  retain mcount;
  if first.id then mcount=0;
  mcount+1;
  if mcount le 2;
run;

proc sql;
  create table cactc_data_all as
  select group, id ,%do i=1 %to &maxexp2; exp&i, %end; startobs, endobs, event1 from cactc_case
  union all
  select group, id ,%do i=1 %to &maxexp2; exp&i, %end; startobs, endobs, event1 from cactc_data02
  order by id, group desc;
quit;

data cactc_data_all2;
  set cactc_data_all;
  by id;
  retain no;
  if first.id then no=0;
  no+1;
  output cactc_data_all2;
run;

data not_cactc;
  set cactc_data_all2;
  by id;
  if last.id then do;
  if no le 2 then output not_cactc;
  end;
run;

proc sort data=cactc_data_all2;
  by id;
run;

proc sort data=not_cactc;
  by id;
run;

data simd37_cactc;
  merge cactc_data_all2 not_cactc(in=b_);
  by id;
  if b_ then delete;
run;

data simd37_cactc;
  set simd37_cactc;
  retain endobs2;
  if group=1 then endobs2=endobs;
  by id;
run;

proc sort data=simd37_cactc;
  by id descending group;
run;

data simd38_cactc(keep=id %do i=1 %to &maxexp2; exp&i  %end; event1 contst conted contst2 conted2 casest caseed end group);
  set simd37_cactc;
  if (event1 < (3*&riskp+&lag) or endobs2 < event1) then delete;
  end=min(endobs2, event1);
  contst=event1 - ((3*&riskp + &lag)-1);
  conted=event1 - (2*&riskp + &lag);
  contst2=event1 - ((2*&riskp+&lag)-1);
  conted2=event1 - (&riskp + &lag);
  casest=event1 - (&riskp-1);
  caseed=event1;
run;

data simd38_cactc;
  set simd38_cactc;
  by id;
  if first.id then seq=0;
  seq+1;
run;

data simd39_cactc;
  set simd38_cactc;
  %do i=1 %to &maxexp2;
  cont&i =(contst <= exp&i and exp&i <= conted);
  %end;
  %do i=1 %to &maxexp2;
  cont2&i =(contst2 <= exp&i and exp&i <= conted2);
  %end;
  %do j=1 %to &maxexp2;
  case&j = (casest <= exp&j and exp&j <= caseed);
  %end;
  cont=0
  %do i=1 %to &maxexp2;
  + cont&i
  %end;;
  cont_2=0
  %do i=1 %to &maxexp2;
  + cont2&i
  %end;;
  case=0
  %do j=1 %to &maxexp2;
  + case&j
  %end;;
run;

data simd40_cactc(keep=id exp group seq);
  set simd39_cactc;
  exp=1;
  if case => 1 then output;
run;

proc sort data=simd40_cactc;
  by ID seq;
run;

proc sql;
  create table simd41_cactc as
  select a.id, b.exp, a.group, a.seq
  from simd38_cactc a left join simd40_cactc b on a.ID=b.ID and a.seq=b.seq;
quit; 

data simd42_cactc;
  set simd41_cactc;
  period=1;
  cactcid=1;
  if exp=. then exp=0;
run;

data simd43_cactc(keep=id exp group seq);
  set simd39_cactc;
  exp=1;
  if cont => 1 then output;
run;

proc sort data=simd43_cactc;
  by ID seq;
run;

proc sql;
  create table simd44_cactc as
  select a.id, b.exp, a.group, a.seq
  from simd38_cactc a left join simd43_cactc b on a.ID=b.ID and a.seq=b.seq;
quit;  

data simd45_cactc;
  set simd44_cactc;
  period=0;
  cactcid=2;
  if exp=. then exp=0;
run;

data simd46_cactc(keep=id exp group seq);
  set simd39_cactc;
  exp=1;
  if cont_2 => 1 then output;
run;

proc sort data=simd46_cactc;
  by ID /*group*/ seq;
run;

proc sql;
  create table simd47_cactc as
  select a.id, b.exp, a.group, a.seq
  from simd38_cactc a left join simd46_cactc b on a.ID=b.ID and a.seq=b.seq;
quit;  

data simd48_cactc;
  set simd47_cactc;
  period=0;
  cactcid=3;/*Add20231223*/
  if exp=. then exp=0;
run;

proc sort data=simd48_cactc ;
  by ID /*group*/ seq;
run;

data cco01logistic_cactc;
  set simd42_cactc simd45_cactc simd48_cactc ;
  time=2-period;
run; 

proc sort data=cco01logistic_cactc ;
  by id seq cactcid;
run;

data cco01logistic_cactc;
  set cco01logistic_cactc;
  count=_n_;
  cactc_id=ceil(count/3);
run;

proc logistic data=cco01logistic_cactc outest=outcco_cactc covout descending;
  strata cactc_id;
  model period=exp exp*group;
run;

proc transpose data=outcco_cactc out=outcco2_cactc0 (drop=_NAME_) prefix=param;
  var exp expgroup;
run;

proc transpose data=outcco_cactc out=outcco2_cactcconv (drop=_NAME_ _LABEL_) prefix=conv;
  var _status_;
run;

data outcco2_cactc0;
  set outcco2_cactc0;
  id=_N_;
run;

data outcco2_cactc1(drop=_label_) outcco2_cactc2(drop=_label_);
  set outcco2_cactc0;
  if id=1 then output outcco2_cactc1;
  if id=2 then output outcco2_cactc2;
run;

data outcco2_cactc1(drop=id param3);
  set outcco2_cactc1;
  kekka='timetrend';
run;

data outcco2_cactc2(drop=id param2);
  set outcco2_cactc2;
  kekka='cactc';
  rename param3=param2;
run;

data outcco2_cactcpara;
  set outcco2_cactc1 outcco2_cactc2;
run;

proc sql;
  create table outcco2_cactc
  as select *
  from outcco2_cactcpara, outcco2_cactcconv;
quit;

data outcco3_cactc (keep=kekka estimate sterr ll ul expest expll expul conv1 conv2 conv3);
  set outcco2_cactc;
  estimate=param1;
  sterr=sqrt(param2);
  ll=param1-1.96*sqrt(param2);
  ul=param1+1.96*sqrt(param2);
  expest=exp(estimate);
  expll=exp(ll);
  expul=exp(ul);
run;

proc append base=result_cctc data=outcco3_cactc;
run;
%end;

%mend sim_exptimetrend_sc1;

%let macrofile =C:\Users\sank1\OneDrive\デスクトップ\File;
%include "&macrofile\sccs.sas";
%include "&macrofile\element.sas";
%include "&macrofile\poisreg.sas";

/*Exposure Parameters*/
%let theta0 = log(1.0*10**(-3)); /*Baseline incidence of exposure (fixed)*/
%let theta1 = log(2.0); /*Effect of Cov1 (fixed)*/
%let theta2 = log(1.5); /*Effect of Cov2 (fixed)*/
%let trend =U; /*Exposure time trend  U: No Trend or Linear, C: Cosine fuction-pattern, M: Mountai-pattern*/
%let thetat =1.00; /*Parameter value of αTR No Trend: 1.0, Linear Trend x1.1: 1.1, x1.2: 1.2, x1.5: 1.5, x2: 2.0, x3: 3.0, x5: 5.0*/
%let thetam = log(10**(-3));/*Parameter of Mountain-pattern temporal change*/
%let exprisk=15;/*Length of exposure risk period*/

/*Event Parameters*/
%let beta0 = log(0.0003); /*Baseline incidence of event (fixed)*/
%let beta1 = log(1.0); /*Effect of exposure: x1.0, x1.5, or x3.0*/
%let beta2 = log(2.0); /*Effect of Cov1 (fixed)*/
%let beta3 = log(3.0); /*Effect of Cov2 (fixed)*/
%let riskp = 15;/*Length of risk period for analysis*/
%let pcens = 0.0;/*Probability of censoring at event time*/
%let pcens3 = 0.0;/*In Scenario 3, specify the cessation of exposure after the event.*/
%let lag=30;/*Lag period between Focal window and Referrence windows*/
%let sampleunit=2;/*Sampling units per 1 case*/
%let cactcday=30;/*The Selection criteria for sampling control in CCTC: Setting an appropriate case time point*/

%sim_exptimetrend_sc1(replis=1, replin=1, Nday=365, Nsample=100000);
