%macro I_Squared(	dataset=,
			study_id=,
			exposure=,
			paraest_mat=,
			g_mat=
		)
		/ minoperator;

/*****************/
/*  GNU licenses */
/*****************/

/* 
	Copyright © 2017 Zhaoxue Tong
	Zhaoxue Tong, Renmin University of China, October 2017
	
	This SAS macro is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.
	
	This SAS macro is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
	
	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*******************/
/* About the macro */
/*******************/

/* 
   PURPOSE:
	Computes estimated I-squared statistic for one-stage individual participant data 
	meta-analysis with binary outcomes
   REQUIRED ARGUMENTS:
	dataset		= the input dataset
	study_id	= the variable in the dataset identifying each study
	exposure	= the variable in the dataset indicating the exposure status
	paraest_mat	= the name given to the ParameterEstimates= data set in the ODS OUTPUT statement;
			  stores the estimated parameters beta0 and beta1
	g_mat		= the name given to the G= data set in the ODS OUTPUT statement;
			  stores the estimated covariance matrix of Mu_0 and Mu_1
   
   NOTES:
	Please do the following with your original PROC GLIMMIX step:
		- include 's' or 'solution' as an option in the MODEL statement so that the output data set 
		  containing the fixed effects solutions will be produced
		- include 'g' as an option in the RANDOM statement so that the output data set containing the 
		  random effects covariance matrix will be produced
		- include a statement of the form 'ods output ParameterEstimates=<ds1> G=<ds2>;' to produce the 
		  input data sets for this macro
	Example:
		proc glimmix 
			data = try;
  			model label = x1 / solution;
  			random int x1 / sub = subj_id type = un g;
  			ods output ParameterEstimates = para_est G = g_matrix;
		run;
   OUTPUT:
	Estimated I-squared statistic
*/

/* Generate macro variable for output */
data _null_;
set &dataset;
call symput('exposure_name',vlabel(&exposure));
run; 

/* Check whether all the required macro parameters are specified */
%macro I_squared_checkpara(paravalue, paraname);
	%if %length(&paravalue) = 0 %then %do;
		%put ERROR: Value for macro parameter &paraname is missing;
		%abort;
		%end;
%mend I_squared_checkpara;

%I_squared_checkpara(&dataset, dataset);
%I_squared_checkpara(&study_id, study_id);
%I_squared_checkpara(&exposure, exposure);
%I_squared_checkpara(&paraest_mat, paraest_mat);
%I_squared_checkpara(&g_mat, g_mat);

/* Check whether the data set exists */
%macro I_squared_checkdataset(dset);
	%if %sysfunc(exist(&dset)) = 0 %then %do;
	%put ERROR: data set &dset does not exist;
	%abort;
	%end;
%mend I_squared_checkdataset;

%I_squared_checkdataset(&dataset);
%I_squared_checkdataset(&paraest_mat);
%I_squared_checkdataset(&g_mat);
	
%local dsid anobs whstmt n_obs n_aver;
%let n_obs = 0;

/* Make sure the data set can be opened*/
%let dsid = %sysfunc(open(&dataset, IS));
%if &dsid = 0 %then %do;
	%put ERROR: data set &dataset cannot be opened;
	stop;
	%end;

/* Get the number of studies in the input dataset */
%local n_study;
proc sql noprint;
	select count(distinct(&study_id)) into:n_study
		from &dataset;

/* Get the number of observations in the input dataset */

/* If SAS knows how many observations are in the data set, */
/* and if there’s no WHERE clause, */
/* you can get the number of observations directly without counting */
%let anobs = %sysfunc(attrn(&dsid, ANOBS));
%let whstmt = %sysfunc(attrn(&dsid, WHSTMT));
%if &anobs = 1 & &whstmt = 0 %then %do;
	%let n_obs = %sysfunc(attrn(&dsid, NLOBS));
	%end;

/* If SAS doesn't know the number of observations, */
/* or if you’re using a WHERE clause, */
/* you can obtain the answer by iterating through the data set */
%else
	%do %while (%sysfunc(fetch(&DSID)) = 0);
	%let n_obs = %eval(&counted. + 1);
	%end;

/* Get the average number of observations among all of the studies */
%local n_aver;
%let n_aver = %sysevalf(&n_obs / &n_study);

proc iml;
	/* specify estimated parameters */
	use &paraest_mat;
		read all var {Estimate} into Est;
	beta0_est = Est[{1}, ];
	beta1_est = Est[{2}, ];

    /* specify population mean and covariance */
	Mean = {0, 0};
	use &g_mat;
		read all var {Col1 Col2} into Cov;

	/* check whether Cov is positive definite */
	G = root(Cov, "NoError");
	if G = . then do;
		print
"The covariance matrix should be positive definite";
		stop;
	end;
	eigval = eigval(Cov);
	if any(eigval<=0) then do;
		print
"The covariance matrix should be positive definite";
		stop;
	end;	
	
	/* simulation*/
	call randseed(4321);
	Mu = RandNormal(&n_obs, Mean, Cov);
	Mu_0 = Mu[ ,{1}];
	Mu_1 = Mu[ ,{2}];

	v_1 = Mu_1 + beta1_est;
	v1 = var(v_1);

	use &dataset;
	read all var {&exposure} into x;
	y_est = (beta0_est + Mu_0) + (beta1_est + Mu_1)#x;
	p_est = 1 / (1+exp(-y_est));
	v_2 = 1 / (&n_aver*p_est#(1-p_est));
	v2 = mean(v_2);

	title "I Squared";
	i_squared_est = v1 / (v1 + v2);
	print "&exposure_name" i_squared_est v1 v2;

	title;

%let rc = %sysfunc(close(&dsid));

quit;

%mend I_Squared;
