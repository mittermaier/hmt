function aus=nu_LiBr_vec(T,xi)
% T                        temperature            	C (Celcius)
%                                                	[ ... ]
% xi                       mass fraction         	kg LiBr/kg solution 
%                                                	[0.4 ... 0.7]
% nu			   kinematic viscosity     	in m^2/s
%
% This function calculates the kinematic viscosity of aqueous LiBr
% solution. Input parameters are: temperature and mass fraction of LiBr
%
% The calculation is based on data from the PHD thesis: 
% Löwer, Harald 
% Thermodynamische und physikalische Eigenschaften der wässrigen Lithiumbromid-Lösung Karlsruhe
% 1960

% The regression alaysis was conducted in 2009 during a PHD study by:
% Wohlfeil, Arnold
% Wärme- und Stoffübertragung bei der Absorption an Rieselfilmen in Absorptionskälteanlagen
% Technische Universität Berlin
% 2009

% Change log 
% 2017-11-23 	first release		M.Mittermaier

T=T+273.15; %in K!
a1=-771.2238243;
a2=1.313986647;
a3=6661.751115;
a4=-0.176584923;
a5=-0.00049714;
a6=-327.2586134;
a7=0.000136026;
a8=0.069697926;
a9=-0.000103042;
a10=283.1165871;
a11=-5.806649869;
a12=-5.412275708;
a13=-6824.256192;
a14=-2993.737847;
a15=-816.7742305;
a16=35.8541012;
aus=1e-6*exp(a1+a2*T+a3*xi+a4*T.*xi+a5*T.^2+a6*xi.^2+a7*T.^2.*xi+a8*T.*xi.^2+a9*T.^2.*xi.^2+a10*log(T)+a11*log(T).^2+a12*log(T).^3+a13*log(1+xi)+a14*log(1+xi).^2+a15*log(1+xi).^3+a16*log(T).*log(1+xi));


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
