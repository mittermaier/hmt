function aus=cP_LiBr_vec(T,x)
%aus=dHdT(T,x);
%
% T                        temperature            	C (Celcius)
%                                                	[ ... ]
% x                        mass fraction         	kg LiBr/kg solution 
%                                                	[0.4 ... 0.7]
% cP			   specific heat capacity  	in J/kg K
%
% This function calculates the specific heat capacity of aqueous LiBr
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
% 2017-11-23 	Made suitable for vectors, first release		M.Mittermaier

T=T+273.15;
a0=6.462731914;
a1=-68.15825241;
a2=-0.017426854;
a3=0.520285681;
a4=5.800384892;
a5=4.1611E-05;
a6=-0.055755167;
a7=-0.001404767;
a8=0.000168802;
a9=28.85672066;
a10=-2.95603E-08;
a11=-0.19710322;
a12=0.000474334;
a13=1.25375E-06;
a14=-1.81967E-07;
a15=-3.38265E-07;
aus=a0+a1*x+a2*T+a3*x.*T+a4*x.^2+a5*T.^2+a6*x.^2.*T+a7*x.*T.^2+a8*x.^2.*T.^2+a9*x.^3+a10*T.^3+a11*x.^3.*T+a12*x.^3.*T.^2+a13*x.*T.^3+a14*x.^2.*T.^3+a15*x.^3.*T.^3;
aus=aus*1000;

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
