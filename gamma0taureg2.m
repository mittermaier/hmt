function gamma0taureg2=gamma0taureg2(tau, pic)
% First derivative in tau of ideal-gas part of fundamental equation for region 2
% Stefan Petersen
% 10/10/02

Stoffparameter_H2O;

% tau=ones(1,9)*tau;
% gamma0taureg2=sum(n0reg2.*j0reg2.*tau.^(j0reg2-1));



tauvec(1:length(tau),1:length(n0reg2))=repmat(tau,1,length(n0reg2));

n0reg2=repmat(n0reg2,length(tau),1);
j0reg2=repmat(j0reg2,length(tau),1);


gamma0taureg2_tmp=n0reg2.*j0reg2.*tauvec.^(j0reg2-1);

gamma0taureg2=sum(gamma0taureg2_tmp')';
