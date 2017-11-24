function pSatW=pSatW(temp)

% temp                     temperature            	K (Kelvin) valid from 273 K to 647 K
%                                                	
% pSatW			   saturation pressure		in bar
%
% This function calculates the equilibrium pressure of water based on the temperature input

% Change log 
% 2002-08-10	programmed at TU Berlin 		S.Petersen

Stoffparameter_H2O;
del = temp + nreg4(9) ./ (temp - nreg4(10));
Aco = del .^ 2 + nreg4(1) .* del + nreg4(2);
Bco = nreg4(3) .* del .^ 2 + nreg4(4) .* del + nreg4(5);
cco = nreg4(6) .* del .^ 2 + nreg4(7) .* del + nreg4(8);
pSatW = (2 .* cco ./ (-Bco + (Bco .^ 2 - 4 .* Aco .* cco) .^ 0.5)) .^ 4 .* 10;
