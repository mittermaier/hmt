function absorptionsenthalpie2 = absorptionsenthalpie2_PK(temp, konz)

%
% temp                     	temperature           	in K
%                                                	[273.15 ... 463.15]
% konz                     	mass fraction         	in kg LiBr/kg solution 
%                                                	[0.4 ... 0.75]
% absorptionsenthalpie2_PK    	heat of sorption  	in kJ/kg
%
% This function computes the heat of sorption when steam is absorbed by a 
% LiBr-solutiion at the temperature 'temp' and the LiBr mass fraction 'konz'.
%
% The calculation was found in:
% Å. Jernqvist, H. Kockum 
% Simulation of falling film absorbers and generators (1996)
% P. 313 Equation (12) and (13)
%

% Change log 
% 2011-05-24	programmed		P.Schulze
% 2017-11-23 	first release		M.Mittermaier

%
% enthalpy of solution:
h = enthalpySatLiqTH2OLiBr_PK(temp, konz);
%
h_grad = dh_dw_PK(temp,konz);
% enthalpy of saturated steam is computed with function "enthalpySatVapTW"

hv = enthalpySatVapTW(temp);
absorptionsenthalpie2 = hv - (h + konz.*h_grad);


end;

