function aus=T_Siede(var1xi, var2p)
% var1xi                   mass fraction         	kg LiBr/kg solution 
%                                                	[0.4 ... 0.7]
% var2p                    pressure	         	Pa 
%                                                	
% T_Siede		   saturation temperature	C (celcius)
%
% This function calculates the saturation temperature (thermodynamic equilibrium temperature) of aqueous LiBr
% solution. Input parameters are: mass fraction of LiBr and the pressure 
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

a0=-0.00470858;
a1=-0.001276757;
a2=0.000145597;
a3=0.000428261;
a4=0.000948526;
a5=3.47501E-06;
a6=-0.000495401;
a7=-5.44472E-05;
a8=0.000110477;
a9=0.004915398;
a10=-7.21234E-08;
a11=-0.00058121;
a12=-2.23738E-05;
a13=2.39788E-06;
a14=-6.64049E-06;
a15=4.26683E-06;

lnp=ln(var2p);
hilf=a0+a1*var1xi+a2*lnp+a3*var1xi*lnp+a4*var1xi^2+a5*lnp^2+a6*var1xi^2*lnp+a7*var1xi*lnp^2+a8*var1xi^2*lnp^2+a9*var1xi^3+a10*lnp^3+a11*var1xi^3*lnp+a12*var1xi^3*lnp^2+a13*var1xi*lnp^3+a14*var1xi^2*lnp^3+a15*var1xi^3*lnp^3;
aus=-1/hilf-273.15;

