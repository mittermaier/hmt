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


