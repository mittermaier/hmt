function ausgabe=rho_LiBr_vec(T,xi)
% T                        temperature            	C (Celcius)
%                                                	[ ... ]
% xi                       mass fraction         	kg LiBr/kg solution 
%                                                	[0.4 ... 0.7]
% nu			   density     			in kg/m^3
%
% This function calculates the density of aqueous LiBr
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
a0=-158.792071;
a1=18502.67064;
a2=9.79173601;
a3=-158.6456762;
a4=-73846.71379;
a5=-0.025910389;
a6=679.2667079;
a7=0.465762844;
a8=-2.043875777;
a9=86252.5194;
a10=2.04622E-05;
a11=-794.2826466;
a12=2.390046288;
a13=-0.000453566;
a14=0.002046451;
a15=-0.00239635;
a16=2093.332625;
ausgabe=a0+a1*xi+a2*T+a3*xi.*T+a4*xi.^2+a5*T.^2+a6*xi.^2.*T+a7*xi.*T.^2+a8*xi.^2.*T.^2+a9*xi.^3+a10*T.^3+a11*xi.^3.*T+a12*xi.^3.*T.^2+a13*xi.*T.^3+a14*xi.^2.*T.^3+a15*xi.^3.*T.^3+a16*xi.^4;

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
