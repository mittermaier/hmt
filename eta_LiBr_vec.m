function aus=eta_LiBr_vec(T,xi)
% T                        temperature            	C (Celcius)
%                                                	[ ... ]
% xi                       mass fraction         	kg LiBr/kg solution 
%                                                	[0.4 ... 0.7]
% nu			   dynamic viscosity     	in Pa * s
%
% This function calculates the dynamic viscosity of aqueous LiBr
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

aus=nu_LiBr_vec(T,xi).*rho_LiBr_vec(T,xi);
