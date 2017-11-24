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
