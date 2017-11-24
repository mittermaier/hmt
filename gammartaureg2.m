function gammartaureg2=gammartaureg2(tau, pic)
% Second derivative in tau of residual part of fundamental equation for region 2
% Stefan Petersen
% 10/10/02

Stoffparameter_H2O;

% pic=ones(1,43)*pic;
% tau=ones(1,43)*tau-0.5;
% gammartaureg2=sum(nreg2.*pic.^ireg2.*jreg2.*tau.^(jreg2-1));

tauvec(1:length(tau),1:length(nreg2))=repmat(tau-0.5,1,length(nreg2));
picvec(1:length(pic),1:length(nreg2))=repmat(pic,1,length(nreg2));

nreg2=repmat(nreg2,length(tau),1);
jreg2=repmat(jreg2,length(tau),1);
ireg2=repmat(ireg2,length(tau),1);

gammartaureg2_tmp=nreg2.*picvec.^ireg2.*jreg2.*tauvec.^(jreg2-1);

gammartaureg2=sum(gammartaureg2_tmp')';

