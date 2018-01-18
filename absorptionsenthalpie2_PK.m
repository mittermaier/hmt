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


end

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
