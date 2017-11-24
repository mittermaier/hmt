function aus = Filmdicke_hilf(delta99,nu,u_0,g,x)
% delta99		given boundary layer thickness            		in m
%                                                				
% nu 			kinematic viskosity 					in m^2/s
%
% u0 			initial velocity 					in m/s
%										[0.1...0.5]
% g 			acceleration of gravity 				in m/s^2
%										[9.81]
% x 			sought longitudinal coordinate				in m
%
%This function is an auxiliary function to find the longitudinal coordinate
%for a given boundary layer thickness. 

% The calculation is published in the book "Boundary Layer Theory" by Schlichting and
%Gersten (2006)

%The function was programmed 2012 by P.Schulze during his bachelor thesis 
%under supervision of M.Mittermaier 

% Change log 
% 2017-11-23 	first release		M.Mittermaier

aus=delta99-5*(nu*x/(2*g*x+u_0^2)^0.5)^0.5;

