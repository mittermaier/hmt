function enthalpySatLiqTH2OLiBr_PK=enthalpySatLiqTH2OLiBr_PK(temp,konz)

% This function calculates the specific enthalpy of H2O-LiBr-solution 
% temp                 Temperature             in K
%
% konz                 mass fraction           in kg_LiBr/kg_Solution , e.g. 0.4 ... 0.75
%
% enthalpySatH2OLiBr   specific enthalpy       in kJ/kg
%
% Quelle: J. Patek & J. Klomfar, A computationally effective formulation 
%         of the thermodynamic properties of LiBr-H2O solutions from 
%         273 to 500 K over full composition range. Published in IJR


Stoffparameter_PK;
T_0W=221; % K

if length(konz)==1
    konz(1:length(temp),1)=konz(1,1);
end 
if length(temp)==1
    temp(1:length(konz),1)=temp(1,1);
end 

% Conversion of mass fraction to mole fraction 
konz_mol = (konz./M_LiBr)./(konz./M_LiBr + (1-konz)./M_W);


tempvec(1:length(temp),1:length(aTab7))=repmat(temp,1,length(aTab7));
konzvec(1:length(konz_mol),1:length(aTab7))=repmat(konz_mol,1,length(aTab7));

aTab7=repmat(aTab7,length(temp),1);
mTab7=repmat(mTab7,length(temp),1);
nTab7=repmat(nTab7,length(temp),1);
tTab7=repmat(tTab7,length(temp),1);

tauvec=TCritW./(tempvec-T_0W);

mult_tmp=aTab7.*konzvec.^mTab7.*(0.4-konzvec).^nTab7.*(tauvec).^tTab7;

mult = sum(mult_tmp')';

enthalpySatLiqW = enthalpySatLiqTW(temp).*1000.*M_W;               % in J/mol

% enthalpySatLiqTH2OLiBr_PK = (1-konz_mol).*enthalpySatLiqW + enthalpyCritWmol.*mult;  % in J/mol;

c=(konz_mol.*M_LiBr + (1-konz_mol).*M_W);  
enthalpySatLiqTH2OLiBr_PK = ((1-konz_mol).*enthalpySatLiqW + enthalpyCritWmol.*mult)./1000./c; % in kJ/kg;


