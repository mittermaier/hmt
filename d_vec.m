function aus=d_vec(T,Xi)
 
% T                        temperature            	C (Celcius)
%                                                	[ ... ]
% Xi                       mass fraction         	kg LiBr/kg solution 
%                                                	[0.4 ... 0.7]
% D			   mass diffusivity     	in m^2/s
%
% This function calculates the mass diffusivity of aqueous LiBr
% solution. Input parameters are: temperature and mass fraction of LiBr
%
% The calculation is based on the PHD thesis p.189: 

% Kim, Kwang J
% Heat and mass transfer enhancement in absorption cooling 
% Arizona State University
% 1992

% and personal communication between Kwang J. Kim and Felix Ziegler


% Change log 
% 2017-11-23 	Made suitable for vectors, first release		M.Mittermaier

m=1000*Xi/86.845./(1-Xi);
D25=(1.3528+0.19881*m-0.036382*m.^2+0.0020299*m.^3-0.000039375*m.^4)*1e-9;
D=D25.*(T+273.15)./298.15.*eta_LiBr_vec(25,Xi)./eta_LiBr_vec(T,Xi);
aus=D;
