function aus=longitudgridpnt_GS_2(z_vektor,n,nu,g,u0)

% z_vektor		vector containing transversal grid points           	in m
%                                                				[ ... ]
% n			index of sought longitudinal grid point
%                                                				[ ... ]
% nu 			kinematic viskosity 					in m^2/s
%
% g 			acceleration of gravity 				in m/s^2
%										[9.81]
%u0 			initial velocity 					in m/s
%										[0.1...0.5]
%
%This function calculates longitudinale grid points. For the index n 
%the longitudinal grid point is adjusted in a way
%that the boundary layer thickness fits to the transversal grid point of the 
%given vector z_vektor at index n.
%Therefore the root of function 'Filmdicke_hilf' is sought.  

%x is the sought root. It is expected to be between 10^-19 and 1.
%10^-19 is chosen close to zero, but 0 would lead to an error in matlab. 

%The function was programmed 2012 by P.Schulze during his bachelor thesis 
%under supervision of M.Mittermaier 

% Change log 
% 2017-11-23 	first release		M.Mittermaier

aus=fzero(@(x)Filmdicke_hilf(z_vektor(1,n),nu,u0,g,x),[10^-19 1]);




