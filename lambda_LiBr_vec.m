function aus=lambda_LiBr_vec(T,xi)

% T                        temperature            	C (Celcius)
%                                                	[ ... ]
% xi                       mass fraction         	kg LiBr/kg solution 
%                                                	[0.4 ... 0.7]
% nu			   thermal conductivity		in W/(m*K)
%
% This function calculates the thermal conductivity of aqueous LiBr
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
a1=-2.828546042;
a2=0.027095169;
a3=-7.3592E-05;
a4=6.84271E-08;
a5=0.168719769;
a6=0.371295842;
a7=0.171540919;
a8=-0.002034796;
a9=2.4103E-06;
a10=-0.003341468;
a11=4.93796E-06;
aus=a1+a2*T+a3*T.^2+a4*T.^3+a5*xi+a6*xi.^2+a7*xi.^3+a8*T.*xi+a9*T.^2.*xi+a10*T.*xi.^2+a11*T.^2.*xi.^2;
