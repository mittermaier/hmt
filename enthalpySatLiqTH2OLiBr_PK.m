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

%Copyright 2017 Martin Mittermaier
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%    Dieses Programm ist Freie Software: Sie können es unter den Bedingungen
%    der GNU General Public License, wie von der Free Software Foundation,
%    Version 3 der Lizenz oder (nach Ihrer Wahl) jeder neueren
%    veröffentlichten Version, weiterverbreiten und/oder modifizieren.
%
%    Dieses Programm wird in der Hoffnung, dass es nützlich sein wird, aber
%    OHNE JEDE GEWÄHRLEISTUNG, bereitgestellt; sogar ohne die implizite
%    Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.
%    Siehe die GNU General Public License für weitere Details.
%
%    Sie sollten eine Kopie der GNU General Public License zusammen mit diesem
%    Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.
