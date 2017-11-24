function gammataureg1=gammataureg1(tau, pic)
% gammataureg1(tau, pic)
%
% tau    reduced temperature     dimensionless
% pic    reduced pressure        dimensionless
%
% IAPWS water properties 1997

% Unterfunktion zu den Stoffwerten für Wasser
% Stefan Petersen 10/10/02

Stoffparameter_H2O;

tauvec(1:length(tau),1:length(nreg1))=repmat(tau-1.222,1,length(nreg1));
picvec(1:length(pic),1:length(nreg1))=repmat(7.1-pic,1,length(nreg1));

nreg1=repmat(nreg1,length(tau),1);
ireg1=repmat(ireg1,length(tau),1);
jreg1=repmat(jreg1,length(tau),1);

gammataureg1_tmp=nreg1.*picvec.^ireg1.*jreg1.*tauvec.^(jreg1-1);

gammataureg1=sum(gammataureg1_tmp')';
