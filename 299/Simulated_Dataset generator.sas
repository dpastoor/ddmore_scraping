dm 'clear log';
dm 'clear list';
libname output "C:\Users\am00402567\Desktop\Conc QT Modeling\Misspecification Example";
options ls=130;
proc datasets kill; run; quit;
title;
ods listing;

data set1;
	do doselabel = 1 to 6;
		if doselabel = 1 then dose = 0;
		if doselabel = 2 then dose = 100;
		if doselabel = 3 then dose = 500;
		if doselabel = 4 then dose = 1000;
		if doselabel = 5 then dose = 1500;
		if doselabel = 6 then dose = 1750;
		if dose NE 0 then trt = 1; else trt = 0;
		if dose = 0 then trtstr = "Placebo"; else trtstr = "Active";
		do subjectctr = 1 to 6;
			sid = (doselabel)*100 + subjectctr;
			ka = 0.7*exp(rannor(3453434)*0.1);
			cl = 12000*exp(rannor(3456)*0.1);
			v = 400000*exp(rannor(774)*0.2);

			baseline = round(400 + rannor(324323)*12, 0.1);
			intercept = -2.5 + rannor(345343)*4;
			trteffect = 0 + 0.2*rannor(454);
			chgbaseeffect = -0.1 + rannor(234234)*0.02;
			chgbase = (baseline - 400);
			emax = 20 + rannor(3656)*2;
			ec50 = 400*exp(rannor(345343)*0.05);

			timeEffect1 = 0 + 0.2*rannor(34534);
			timeEffect2 = 0.1 + 0.2*rannor(34534);
			timeEffect3 = -0.2 + 0.2*rannor(34534);
			timeEffect4 = 0.3 + 0.2*rannor(34534);
			timeEffect5 = 0.1 + 0.2*rannor(34534);
			timeEffect6 = -0.3 + 0.2*rannor(34534);
			timeEffect7 = 0.1 + 0.2*rannor(34534);
			timeEffect8 = 0.05 + 0.2*rannor(34534);
			timeEffect9 = -0.2 + 0.2*rannor(34534);

			do timectr = 1 to 9;
				if timectr = 1 then nomtime = 0.5;
				if timectr = 2 then nomtime = 1;
				if timectr = 3 then nomtime = 2;
				if timectr = 4 then nomtime = 3;
				if timectr = 5 then nomtime = 4;
				if timectr = 6 then nomtime = 6;
				if timectr = 7 then nomtime = 8;
				if timectr = 8 then nomtime = 12;
				if timectr = 9 then nomtime = 24;
				time1 = 0; time2 = 0; time3 = 0; time4 = 0; time5 = 0; time6 = 0; time7 = 0; time8 = 0; time9 = 0;
				if nomtime = 0.5 then time1 = 1;
				if nomtime = 1 then time2 = 1;
				if nomtime = 2 then time3 = 1;
				if nomtime = 3 then time4 = 1;
				if nomtime = 4 then time5 = 1;
				if nomtime = 6 then time6 = 1;
				if nomtime = 8 then time7 = 1;
				if nomtime = 12 then time8 = 1;
				if nomtime = 24 then time9 = 1;
				timestr = put(nomtime, 5.1);
				if dose = 0 then trueconc = 0; 
					else trueconc = dose*1E6/v * (ka/(ka-(cl/v))) * (exp(-(cl/v)*nomtime) - exp(-ka*nomtime));
				conc = trueconc * exp(rannor(24454)*0.01);
			 	drugeffect = emax*trueconc/(ec50 + trueconc);
				dqtcf =  intercept + trteffect + timeEffect1*time1 - timeEffect2*time2 + timeEffect3*time3 
					+ timeEffect4*time4 - timeEffect5*time5 + timeEffect6*time6 + timeEffect7*time7 + timeEffect8*time8 + timeEffect9*time9
					+ chgbaseeffect*chgbase + drugeffect + rannor(2432)*4;
				conc = round(conc, 0.01);
				dqtcf = round(dqtcf, 0.01);
				total = baseline + dqtcf;
				output;
			end;
		end;
	end;
	drop cl v ka emax ec50 subjectctr timectr intercept trteffect;
	label trt = "Treatment";
	label timestr = "Nominal Time (h)";
	label trtstr = "Treatment";
	label dqtcf = "Observed Change From Baseline QTcF (msec)";		
	label conc = "Drug Concentration";
	label sid = "Subject Identifier";
	label baseline = "Baseline QTcF Interval (msec)";
	label nomtime = "Nominal Time (h)";
	label dose = "Dose";
proc sort; by sid; run; quit;
data output.concqt;	set set1;run; quit;


proc sort; by sid dose; run; quit;
proc means max noprint; var conc; by sid dose; output out=cmax max=cmax; run; quit;
proc sort data=cmax; by dose; run; quit;
proc means; var cmax; by dose; run; quit;
proc freq; tables nomtime*dose; run;

proc sort data=set1; by dose; run; quit;
proc means data=set1; var dqtcf baseline total; by dose; run; quit;


data unique;
	set set1;
	by sid;
	if first.sid then output;
proc sort; by dose; run; quit;
proc freq; table dose; run; quit;
proc means; var baseline; run; quit;
proc corr; var baseline chgbaseeffect; run;


proc nlmixed data=set1;
	parms theta1=-3.5, theta2=12, theta3=0, theta4=0, theta5=0, theta6=0, theta7=0, theta8=0, theta9=0, theta10=0, theta11=0.1, theta12=20, theta13=2, std1=2.2, std2=0.001, sigma=5;
	bounds std1>0, std2>0, sigma > 0;
	intercept = theta1 + eta1;
	emax = theta12 + eta2;
	ec50 = theta13;
	dv = conc;
	if conc = 0 then drugeffect = 0; 
		else drugeffect = emax*dv/(ec50 + dv);
	pred = intercept + theta2*trt + theta3*time1 + theta4*time2 + theta5*time3 + theta6*time4
		+ theta7*time5 + theta8*time6 + theta9*time7 + theta10*time8 + theta11*chgbase + drugeffect;
	residuals = pred - dqtcf;
	model dqtcf ~ normal(pred, sigma*sigma);
	estimate 'Emax' theta12;
	estimate 'EC50' theta13;
	estimate 'Variance(ETA1)' std1*std1;
	estimate 'Variance(ETA2)' std2*std2;
	estimate 'Residual Variance' sigma*sigma;
	random eta1 eta2  ~ normal([0,0], [std1*std1, 0, std2*std2]) subject=sid;
run; quit;


ods graphics on/ reset width=6 in border=off reset=index;
ods html style=statistical image_dpi =200;
	proc sgplot data=set1;
		title;
		refline 0 / axis=y	lineattrs=(color=gray thickness=1 pattern=2);
		scatter x=conc y=dqtcf / group=dose;
		loess x=conc y=dqtcf / lineattrs=(color=black thickness=2) legendlabel="LOESS Fit" nomarkers;
		reg x=conc y=dqtcf / nomarkers;
		keylegend / across=1 position=se location=inside;
		xaxis label = "Drug Concentration (ng/mL)";
		yaxis label = "Observed Change From Baseline QTcF Interval (msec)";
	run; quit;

	proc sgpanel data=set1  noautolegend;
		panelby dose / columns=3 rows=2;
		series x=nomtime y=dqtcf / group=sid markers;
	run; quit;
	proc sgpanel data=set1 noautolegend;
		panelby dose / columns=3 rows=2;
		series x=nomtime y=conc / group=sid markers;
	run; quit;
ods html close;
