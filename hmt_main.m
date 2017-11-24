%
% This is the main program to calculate temperatures, velocities and mass
% fractions during absorption of steam in a laminar lithium bromide film flowing over an isothermal wall.  

%The program is based on the publication: 
%'The impact of viscosity on the combined heat, mass and momentum transfer in laminar liquid falling films'
%from Martin Mittermaier and Felix Ziegler in 2017
%DOI 10.1007/s00231-017-2219-9

%The program is made for open usage. Please cite the publication as soon as your work is based on this program.
%Contact: mittermaier.martin@gmail.com

%Assumptions:
% The velocity at the inlet is considered to be uniform and develops to a parabolic velocity profile typical for laminar falling films
% Initial temperatures and mass fractions are set below

%Overall the following assuptions have been made
%(a) incompressible liquid
%(b) negligible lateral velocity (z-coordinate)
%(c) no slip at the smooth wall
%(d) impermeable and fully wetted, isothermal wall
%(e) no shear stresses at the liquid vapour interface
%(f) constant pressure throughout the entire film
%(g) negligible change in temperature caused by dissipation throughout the film
%(h) negligible diffusion and heat conduction in flow direction as compared to convection
%(i) negligible vapour pressure of the solvent
%(j) steady state problem

%%
% Change log:

%2017-11-20 	first release		M.Mittermaier

%%
clear all

%x - longitudinale coordinate (in flow direction)
%y - transversale coordinate (across flow direction)

%Inputs
y_number=600;  %number of transversal nodes across the film the wall is denoted by index 1, the interface is denoted by index y_number+1

L=0.1; %overall flow length in [m]
u0=0.1; % initial uniform velocitiy at the inlet in longitudinal direction in [m/s]
gamma(1,1)=0.036; %initial mass flow of solution in longitudinal direction (as referred to a non-modelled z-coordinate) in [kg/(ms)]

p=1500; %pressure in [Pa] 

%T0=27+273.15; %choose initial temperature in [K] or choose initial mass
%fraction and use saturation conditions
w0=0.5; %initial mass fraction in [kg H2O/kg solution] used later to determine the saturation temperature 
%w0 = 1-Xi_Siede(T0-273.15,p); %for use with above set temperature to compute the correspondig mass fraction. 
% XiSiede expects input temperatur in [°C] and pressure in [Pa] and returns the mass fraction of LiBr-salt in [kg LiBr /kg solution].
% Since the water mass fraction is used throughout the entire program the
% value Xi_Siede is substracted from 1.


T0=T_Siede(1-w0,p)+273.15 + 3; %'T_Siede' calls the function for saturation temperature based on mass fraction and pressure, here 3 K are added.

Tw=T_Siede(1-w0,p)+273.15 - 6; %wall temperature in [K] (constant over flow length), here 6 K below saturation temperature 'T_Siede'

alpha = 90; %angle of the plate, 90° corresponds to vertical plate
g = 9.80665*sin((alpha/90)*pi/2); % resulting acceleration of gravity in [m^2/s]

%LONGITUDINAL GRID:
%In vicinity to the inlet the longitudinal grid spacing is very small. 
%The computation of longitudinal spacing is done with the help of an analytical equation for the approx. thickness of the boundary layer. 
%(please see Boundary Layer Theory by Schlichting and Gersten (page 160 in 10th edition 2006))
%By doing so, the increase in longitudinal spacing is adjusted in a way that with each consecutive longitudinal node
%the boundary layer thickness grows  only one transversal node further towards the free surface.

%With increasing flow length the boundary layer meets the free surface. 
%Subsequently the longitudinal grid spacing is exponentally increased until 100 times of initial spacing is reached.
%PLEASE NOTE THAT THIS PROCEDURE MIGHT NOT RESULT IN A NUMERICALLY INDEPENDANT SOLUTION!
%THE GROWTH OF LONGITUDINAL SPACING IS ALSO DEPENDANT ON PROPERTY DATA ETC. 
%MAYBE A LINEAR GROWTH IS MORE APPROPRIATE TO THE DATA YOU ARE USING.
%Once the longitudinal grid spacing is increased the steps become equally spaced. 

Faktor=1.1;  %determines the base of exponential growth of longitudinal grid
%=> the higher Faktor, the faster is the growth of the spacing between longitudinal grid points
%and the shorter the distance of transition between the grids)%HANDLE WITH CARE


%STRUCTURE OF LOOPS:
%There are 3 different kinds of loops nested in the program:
%The first loop runs over the entire flow length processing grid point by grid point
%until the flow length L is reachhed. Thus the caculations are conducted sequential.

%The second loop processes the coupling among velocities with the temperatures and mass fractions
%and is called multiple times for each single longitudinal grid point. 

%Within the coupling loop for every single longitudinal grid point one loop for solving the
%velocity profile and one loop for solving the temperature- and mass fraction distribution is executed.
%These are using Newton's method to solve the algebraic finite difference approximation of the PDE system. 

%RESIDUAL DEVIANCE
%The following vectors define the accepted residual deviance for a certain loop. 
%The conditions are formed as row vectors in order to fit the structure of the quantity that it is compared with.

Genauigkeit_U=repmat(10^-12,1,y_number+1);%Accepted residual deviance 
%of longitudinal velocity u (in flow direction) used in Newton loop for
%velocitiy distribution
Genauigkeit_V=repmat(10^-10,1,y_number+1);%Accepted residual deviance
%of transversal velocity v (normal to flow direction) used in Newton loop for
%velocitiy distribution
Genauigkeit_t_sub=repmat(10^-8,1,y_number+1);%Accepted residual deviance  
%of temperature used in Newton loop for the temperature and mass fraction distribution
Genauigkeit_W_sub=repmat(10^-14,1,y_number+1);%Accepted residual deviance  
%of mass fraction used in Newton loop for the temperature and mass fraction distribution
Genauigkeit_t=repmat(10^-7,1,y_number+1);%Accepted residual deviance  
%of the temperature in the coupling loop
Genauigkeit_W=repmat(10^-12,1,y_number+1);%Accepted residual deviance  
%of the mass fraction in the coupling loop
Genauigkeit_Delta=repmat(10^-14,1,y_number+1);%Accepted residual deviance  
%of the film thickness in the coupling loop

%MATRICES
%T - matrix of temperature (contains temperture values for each grid point)
%u - matrix of longitudinal velocities (contains velocities at each grid point in flow direction)
%v - matrix of transversal velocities 
%w - matrix of mass fraction values at each grid point

%INLET CONDITIONS:
%The first index indicates the longitudinale coordinate of the grid e.g. Inlet=1
%The second index indicates the transversale coordinate
T(1,1:y_number+1)=repmat(T0,1,y_number+1);%constante temperature T0 (except of the wall):
T(1,1)=Tw;%at the wall

v(1,1:y_number+1)=zeros(1,y_number+1);%v equals zero at the inlet, since no deceleration
%of longitudinal velocity u has occurred yet. 
%It is not meant as boundary condition but rather a initial value
%to start the Newton's method.
w(1,1:y_number)=repmat(w0,1,y_number);%constant mass fraction w0 at the inlet (except of the interface)
w(1,y_number+1)=1-Xi_Siede(T0-273.15,p);%at the interface

%saturation (corresponding to the inlet temperature) is assumed at the interface 
%Xi_Siede is a function based on data collected by Harald Löwer (1960), an a regression analysis conducted by Arnold Wohlfeil (2009)
%Inputs are temperature in °C and pressure in Pa 
%Output is the saturation mass fraction of LiBr in kg SALT/kg SOLUTION 

%The properties are computed for each grid point
%my - dynamic viscosity in Pas
%rho - density in kg/m^3
%cp - specific heat capacity  in J/(kgK)
%lambda - thermal conductivity in W/(m*K)
%ha - partial massspecific enthalpy of water in J/kg
%hs - partial massspecific enthalpy of LiBr in J/kg
%Diff - diffusion coefficient in m^2/s
%habs - mass specific heat of sorption in J/kg

%properties at the inlet (first "row"):
my(1,1:y_number+1)=eta_LiBr_vec(T(1,:)-273.15,1-w(1,:));
rho(1,1:y_number+1)=rho_LiBr_vec(T(1,:)-273.15,1-w(1,:));
cp(1,1:y_number+1)=cP_LiBr_vec(T(1,:)-273.15,1-w(1,:));
lambda(1,1:y_number+1)=lambda_LiBr_vec(T(1,:)-273.15,1-w(1,:));
ha(1,1:y_number+1)=1000*(enthalpySatLiqTH2OLiBr_PK(T(1,:)',1-w(1,:)')'+(1-w(1,:)).*dh_dw_PK(T(1,:),1-w(1,:)));%ha=h+(1-w)dh_dw 
%TEMPERATURE is in Kelvin here!
hs(1,1:y_number+1)=1000*(enthalpySatLiqTH2OLiBr_PK(T(1,:)',1-w(1,:)')'-w(1,:).*dh_dw_PK(T(1,:),1-w(1,:)));%hs=h-w*dh_dw
Diff(1,1:y_number+1)=d_vec(T(1,:)-273.15,1-w(1,:));
habs(1,1)=1000*absorptionsenthalpie2_PK(T(1,y_number+1),1-w(1,y_number+1));

%The functions eta_LiBr_vec, rho_LiBr_vec, cp_LiBr_vec, lambda_LiBr_vec
%are based on data collected by Harald Löwer (1960), an a regression analysis conducted by Arnold Wohlfeil (2009)
%Inputs ar temperature in °C LiBr mass fraction in kg SALT/kg SOLUTION

%The Function enthalpySatLiqTH2OLiBr_PK is a property function for the enthalpy of the aqueous solution based on 
%Patek und Klomfar (2006)
%Inputs are temperature in K und LIBr mass fraction
%Output is the mass specific enthalpy in kJ/kg 
%dh_dw_PK ist the partial derivative of the specific enthalpy with respect to water mass fraction

%d_vec computes the diffusion coefficient based on data from  Kwang J. Kim (1992) 

%habs computes the mass specific heat of sorption based on the wquation h_abs=h_vapor-ha (partial massspecific enthalpy of water) 
%The data for steam are based on IAPWS (1997) ha is calculated by:(ha=h+(1-w)dh_dw).

rho_quer=sum(rho(1,:),2)/(y_number+1);%mean value at the inlet for calculation of the initial film thickness
my_quer=sum(my(1,:),2)/(y_number+1);%mean value at the inlet for calculation of the initial film thickness
%sum(...,2) (dimension input for rows)

%%
%ANALYTIC NUSSELT FILM THICKNESS and corresponding inlet velocity
%delta_NU=nthroot((3*gamma*my_quer/(g*(rho_quer)^2)),3);
%delta(1,1)=delta_NU;
%u0=gamma(1,1)/(rho_quer*delta_NU);
%%

%INLET VELOCITY
u(1,1)=0;%at the wall: no slip condition =0
u(1,2:y_number+1)=repmat(u0,1,y_number);%constant longitudinal velocity (except of the velocity at the wall)

u_quer=sum(u(1,:),2)/(y_number);%mean value of longitudinal velocity at the inlet
delta(1,1)=gamma(1,1)/(rho_quer*u_quer);%calculation of the film thicknes at inlet (with homogeneous velocity at the inlet)
%%

dy(1,1:y_number)=transversalgridpnt(y_number,delta(1,1));%transversal grid points at the inlet 

%Definition of some abbreviations used in the energy balance later:
psi=ha-hs;
%psi(n,1:y_number+1)=0; %in case the partial heat of mixing should be neglected
phi=rho.*Diff;

%Teta is a matrix containing averaged transversal grid spaces
Teta(1,1:y_number-1)=((dy(1,1:y_number-1)+dy(1,2:y_number))*0.5).^-1;

%z becomes the matrix that saves the transversal coordinates of each row 

z_hilf=tril(repmat(dy(1,1:y_number),y_number+1,1),-1);
%z_hilf contains a (y_number+1) by (y_number)-matrix of the form:
%0  0  0  ...  0
%dy 0  0  ...  0
%dy dy 0  ...  0
% :  : :   :   :
%dy dy dy ...  dy

%Summation of the rows of the matrix leads to a vector z which contains all transversal coordinates 
%e.g. 0,dy,2dy,...,y_number*dy(=delta):
z(1,1:y_number+1)=sum(z_hilf,2);

%A_3 is the vector which contains all longitudinale coordinates of the grid 
A_3(1,1)=0;%at the inlet =0

%check_2 is a kind of switch. It becomes 1 when the boundary layer thickness meets with the interface
check_2=0;%check_2 is initially zero


%n is the row index of the film flow
n=2;%becomes 2 (First grid point after the known inlet)

%%START of the loop over the entire flow length 
while(A_3(n-1,1)<L)% break onecs the flow length L is covered

    i=1;%Counter of coupling loop (starts at 1)
    %t,W and Delta are intermediate substitutes for the 
    %final results of T,w und delta (temperature, mass fraction, film thickness)
    %t,W and Delta are also the objectives within the coupling loop 
    %therefore the values are also used in the break condition of the coupling loop 
    %The three variables only save the values for the current position in the film
    %For comparision of two subsequent iterations two values are saved as intermediat values 
    %Thus t,W have the dimension 2 by (y_number+1) and Delta 2 by 1
    t=zeros(2,y_number+1);%Initialising of t
    W=zeros(2,y_number+1);%Initialising of W
    Delta=zeros(2,1);%Initialising of Delta
    
    %As initial value the substitutes get the final values of the last solution (n-1):
    t(1:2,1:y_number+1)=repmat(T(n-1,:),2,1);
    W(1:2,1:y_number+1)=repmat(w(n-1,:),2,1);
    Delta(1:2,1)=repmat(delta(n-1,1),2,1);
    %The velocities u and v get also the values of row n-1:
    u(n,:)=u(n-1,:);
    v(n,:)=v(n-1,:);
    
    %Start of the coupling loop:
    while (0==0)%Break condition follows further below
        dy(n,1:y_number)=transversalgridpnt(y_number,Delta(mod(i,2)+1));%Computation of transversal grid spacing in row n
        
        z_hilf=tril(repmat(dy(n,1:y_number),y_number+1,1),-1);
        z(n,1:y_number+1)=sum(z_hilf,2);%transversal coordinates row  n
        z_strich=z(n,:);%the n-th row will be cut out of z 
        
        
        %Calculation of properties (with substitute values for T und w)
        my(n,1:y_number+1)=eta_LiBr_vec(t(mod(i,2)+1,:)-273.15,1-W(mod(i,2)+1,:));
        rho(n,1:y_number+1)=rho_LiBr_vec(t(mod(i,2)+1,:)-273.15,1-W(mod(i,2)+1,:));
        cp(n,1:y_number+1)=cP_LiBr_vec(t(mod(i,2)+1,:)-273.15,1-W(mod(i,2)+1,:));
        lambda(n,1:y_number+1)=lambda_LiBr_vec(t(mod(i,2)+1,:)-273.15,1-W(mod(i,2)+1,:));
        ha(n,1:y_number+1)=1000*(enthalpySatLiqTH2OLiBr_PK(t(mod(i,2)+1,:)',1-W(mod(i,2)+1,:)')'+(1-W(mod(i,2)+1,:)).*dh_dw_PK(t(mod(i,2)+1,:),1-W(mod(i,2)+1,:)));
	%TEMPERATURE is in Kelvin here!
        hs(n,1:y_number+1)=1000*(enthalpySatLiqTH2OLiBr_PK(t(mod(i,2)+1,:)',1-W(mod(i,2)+1,:)')'-W(mod(i,2)+1,:).*dh_dw_PK(t(mod(i,2)+1,:),1-W(mod(i,2)+1,:)));
        Diff(n,1:y_number+1)=d_vec(t(mod(i,2)+1,:)-273.15,1-W(mod(i,2)+1,:));
        habs(n,1)=1000*absorptionsenthalpie2_PK(t(mod(i,2)+1,y_number+1),1-W(mod(i,2)+1,y_number+1));


        %% Some abbreviations used for solving the energy balance
        psi=ha-hs;
        phi=rho.*Diff;
        Teta(n,1:y_number-1)=((dy(n,1:y_number-1)+dy(n,2:y_number))*0.5).^-1; %#ok<SAGROW>
        %%
        rho_quer=sum(rho(n,:),2)/(y_number+1);%mean density in row n
        my_quer=sum(my(n,:),2)/(y_number+1);%mean dynamic viscosity in row n
        nu_quer=my_quer/rho_quer; %mean kinematic viscosity in 
        %row n (for calculation of the approx. boundary layer thickness)
        if(check_2~=1) %check, whether the thickness of the boundary layer has already 
            %reached the film thickness
            %longitudinal grid points within the developing velocity region:
            A_3(n,1)=longitudgridpnt_GS_2(z_strich,n,nu_quer,g,u0);
            %longitudgridpnt_GS_2 returns for given transversal points meaning given
            %matrix (z_strich), index (n), kinematic viscosity
            %(nu_quer), gravitational acceleration g and initial velocity u0
            %the respective longitudnale coordinate.
	    %The end of the hydrodynamic boundary layer is then at z_strich(n)
            
            %Approximation of the boundary layer thickness to provide a creterion for the variable check_2
            
            u_unendlich(n-1,1)=(2*g*A_3(n,1)+u0^2)^0.5; %velocity in the undisturbed bulk flow (potential flow, free fall)
            delta_99(n,1)=5*(nu_quer*A_3(n,1)/u_unendlich(n-1,1))^0.5;%Approximation of the boundary layer thickness, 
	    %equation taken from Schlichting and Gersten: "Boundary Layer Theorie")
            delta_99_index=find(min(abs(z_strich(1,:)-delta_99(n,1)))==abs(z_strich(1,:)-delta_99(n,1)));
	    %returns the transversal index in vicinity of the boundary layer thickness 
	    %(should be n since the longitudinal grid has been adjusted to be there)
        else%exponential growth from the inlet reagion or other kind of grid function:
            %When the variable check_2 is set to 1, in the grid function N_2 is set to the value of the last row n, 
            
             A_3(n,1)=A_3(N_2-1,1)+longitudgridpnt_GS_exp_lin(n+2-N_2,A_3(N_2,1)-A_3(N_2-1,1),Faktor);

            %longitudgridpnt_GS_exp_lin calculates the longitudinal grid points based on the last spacing 
            %on the last spacing of the former grid function (A_3(N_2,1)-A_3(N_2-1,1)) 
            %and the basis (Faktor)
            %The spacing in longitudinal direction increases until the spacing is 100 times the inital value
	    %Further on the spacing remains constant. 
	    %PLEASE NOTE:
	    %The increase might be too pronounced for the problem you are investigating. Please check the numerical independance. 
	    %Change the grid function if necessary 
	    %
            %longitudgridpnt_GS_exp_lin initially responds 0, therefore
            %the last result of the function longitudgridpnt_GS_2 has to be added.
            %
            %The first calculations are:
            %A_3(N_2+1)=A_3(N_2-1)+(A_3(N_2)-A_3(N_2-1))*Faktor
            %A_3(N_2+2)=A_3(N_2+1)+(A_3(N_2)-A_3(N_2-1))*Faktor^2
            %A_3(N_2+3)=A_3(N_2+2)+(A_3(N_2)-A_3(N_2-1))*Faktor^3
	    %...
        end
        dx(n-1,1)=A_3(n,1)-A_3(n-1,1);%Calculation of the step size in longitudinal direction (x-direction) in row n
        %
%% 
        o=1;%Counter of the Newton loop for the velocity field
        %U and V are the substitutes of u und v and the objective values of the
        %iterations within the following loop
        U=zeros(2,y_number+1);%Initialising of U as substitute of u
        V=zeros(2,y_number+1);%Initialising of V as substitute of v
        
        U(1:2,1:y_number+1)=repmat(u(n,:),2,1);%Transfer of values from velocities u to U 
        %
        V(1:2,1:y_number+1)=repmat(v(n,:),2,1);%Transfer of values from velocities v to V 
        
        %Start of the Newton loop to solve for velocities 
            while (0==0)%
            %The discretized equations of both the continuity and the boundary layer equation 
            %are arranged from index 2 (1 corresponds to the wall) 
            %up to the index value y_number (the interface corresponds to y_number +1)
            %As solution ones gets both velocity components: longitudinal velocity  u from 2 to
            %y_number and transversal velocity component v from index  3 to y_number+1. 

	    %These velocities are saved
            %in the vector 'uv_Vektor' (first entries are values for v followed by entries for u)
            %The system of equations is written in a matrix notation the way H*uv_vektor=I
            %I represents the inhomogeneity; H represents the
            %matrix of coefficients of the velocities and is devided in 
            %4 quadrants. Each one is made by a 
            %(y_number-1) by (y_number-1) matrix
	    %Since the values for transversal velocities v_i
            %have an offset of one index number as compared to the velocities u_i,
            %the entries of v_i are not placed in the principal diagonal 
            %but in the lower secondary diagonal
            
            %First quadrant of matrix of coefficients (H) (coefficients of v in equation of continuity):
            %principal diagonal(coefficients of v_(i+1)):
            diag_1_2(1:y_number-1)=rho(n,2:y_number)./dy(n,2:y_number);
            %lower secondary diagonal (coefficients of v_i):
            diag_1_3(1:y_number-2)=(rho(n,4:y_number+1)-(2.*rho(n,3:y_number)))./dy(n,3:y_number);
            %the vectors above become the diagonal in the matrix: 
            H_1_2=diag(diag_1_2); %diag_1_2 as principal diagonal
            H_1_3=diag(diag_1_3,-1);%the second input -1 puts the values into the lower secondary diagonal
            %The first quadrant of H (called H_1) is completed by summation:
            H_1=H_1_2+H_1_3;
            
            %Second quadrant of matrix of coefficients (H) (coefficients of u in equation of continuity):
            %principal diagonal (coefficients of u_i):
            diag_2_2(1:y_number-1)=((2.*rho(n,2:y_number))-rho(n-1,2:y_number))./dx(n-1,1);
            %Second quadrant of H (contains principal diagonal only):
            H_2=diag(diag_2_2);
            
            %Third quadrant matrix of coefficients (H) (coefficients of v in 
            %boundary layer equation):
            diag_3_3=zeros(1,y_number-2);%initialising of vector
            %lower secondary diagonal  (coefficients of v_i):
            %equations 3 to y_number-1:
            diag_3_3(1:y_number-3)=(U(mod(o,2)+1,4:y_number)-U(mod(o,2)+1,3:y_number-1)).*rho(n,3:y_number-1)./dy(n,3:y_number-1);
            %equation at y_number (with boundary for u: no shear stress at the free surface):
            diag_3_3(y_number-2)=0;
            H_3=diag(diag_3_3,-1);%third quadrant of H (called H_3)
            
            %Fourth quadrant of matrix of coefficients (H) (coefficients of u in  
            %boundary layer equation):
            diag_4_2=zeros(1,y_number-1);%initialising of vector to speed up matlab a bit
            %principal diagonal (coefficients of u_i):
            %Equations from 2 to y_number-1:
            diag_4_2(1:y_number-2)=(rho(n,2:y_number-1).*(U(mod(o,2)+1,2:y_number-1)-u(n-1,2:y_number-1))./dx(n-1,1))+((my(n,3:y_number)-my(n,2:y_number-1))./(dy(n,2:y_number-1).^2))+(Teta(n,1:y_number-2).*my(n,2:y_number-1)./(dy(n,2:y_number-1)))+(Teta(n,1:y_number-2).*my(n,2:y_number-1)./(dy(n,1:y_number-2)));
            %equation at grid point 'y_number' (with boundary for u: no shear stress at the free surface):
            diag_4_2(y_number-1)=rho(n,y_number)*(U(mod(o,2)+1,y_number)-u(n-1,y_number))/dx(n-1,1)+((Teta(n,y_number-1)*my(n,y_number))/(dy(n,y_number-1)));
            %upper secondary diagonal (coefficients of u_i+1):
            diag_4_1(1:y_number-2)=-(((Teta(n,1:y_number-2).*my(n,2:y_number-1))./dy(n,2:y_number-1))+((my(n,3:y_number)-my(n,2:y_number-1))./(dy(n,2:y_number-1).^2)));
            %lower secondary diagonal (coefficients of u_i-1):
            diag_4_3(1:y_number-2)=-(Teta(n,2:y_number-1).*my(n,3:y_number))./(dy(n,2:y_number-1));
            %building of matrices from the vectors:
            H_4_2=diag(diag_4_2);
            H_4_1=diag(diag_4_1,1);%the second input 1 puts the values into the upper secondary diagonal
            H_4_3=diag(diag_4_3,-1);
            H_4=H_4_1+H_4_2+H_4_3;%the forth quadrant of H (called H_4)
            
            %The full matrix of coefficients contains the 4 quadrants:
            H=[H_1 H_2; H_3 H_4];
            
            %uv_vektor contains the results of longitudinal and transversal velocities 
            %Since the intermediate results have to be saved the vector contains  
            %the values V_i from i=3 to i=y_number+1 and the values of U_i from 
            %i=2 to i=y_number:
            uv_vektor=[V(mod(o,2)+1,3:y_number+1)'; U(mod(o,2)+1,2:y_number)'];
            
            %I is a column vector which contains the inhomgeniteity of both
            %equations, continuity and boundary layer (from index i=2
            %to i=y_number):
            I=[rho(n,2:y_number)'.*u(n-1,2:y_number)'/dx(n-1,1); rho(n,2:y_number)'*g];
            
            %To use Newton's method of equation solving, the Jacobi matrix has to be formed (variabel J)  
            %Matrix J is of the same size and organisation as the matrix H  
            %
             
            %First quadrant of J (partial derivatives of v in  
            %equation of continuity):
            %principal diagonal(partial derivatives with respect to v_(i+1)):
            diag_1_2(1:y_number-1)=rho(n,2:y_number)./dy(n,2:y_number);
            %lower secondary diagonal (partial derivatives with respect to v_i):
            diag_1_3(1:y_number-2)=(rho(n,4:y_number+1)-(2.*rho(n,3:y_number)))./dy(n,3:y_number);
            %Creating matrices out of the diagonals:
            J_1_2=diag(diag_1_2);
            J_1_3=diag(diag_1_3,-1);
            J_1=J_1_2+J_1_3;%First quadrant of Jacobi matrix
            
            %Second quadrant of J (partial derivatives of u in 
            %equation of continuity):
            %principal diagonal(partial derivatives with respect to u_i):
            diag_2_2(1:y_number-1)=((2.*rho(n,2:y_number))-rho(n-1,2:y_number))./dx(n-1,1);
            J_2=diag(diag_2_2);%Second quadrant of Jacobi matrix
            
            %Third quadrant of J (partial derivatives of v in 
            %boundary layer equation):
            %lower secondary diagonal (partial derivatives with respect to v_i):
            %equations from 3 to y_number-1:
            diag_3_3(1:y_number-3)=(U(mod(o,2)+1,4:y_number)-U(mod(o,2)+1,3:y_number-1)).*rho(n,3:y_number-1)./dy(n,3:y_number-1);
            %equation at the grid point 'y_number' (with boundary for u: no shear stress at the free surface):
            diag_3_3(y_number-2)=0;
            J_3=diag(diag_3_3,-1);%Third quadrant of J
            
            %Fourth quadrant of J (partial derivatives of u in 
            %boundary layer equation):
            %principal diagonal (partial derivatives with respect to u_i):
            %equation at the grid point '2' (with boundary: v_2=0):
            diag_4_2(1)=(rho(n,2).*((2*U(mod(o,2)+1,2)-u(n-1,2))/dx(n-1,1)))+((my(n,3)-my(n,2))/(dy(n,2)^2))+((Teta(n,1)*my(n,2))/(dy(n,2)))+((Teta(n,1)*my(n,2))/(dy(n,1)));
            %equations from grid points '3 to y_number':
            diag_4_2(2:y_number-2)=(rho(n,3:y_number-1).*((2*U(mod(o,2)+1,3:y_number-1)-u(n-1,3:y_number-1))./dx(n-1,1)))+((my(n,4:y_number)-my(n,3:y_number-1))./(dy(n,3:y_number-1).^2))+((Teta(n,2:y_number-2).*my(n,3:y_number-1))./dy(n,3:y_number-1))+((Teta(n,2:y_number-2).*my(n,3:y_number-1))./dy(n,2:y_number-2))-((rho(n,3:y_number-1).*V(mod(o,2)+1,3:y_number-1))./(dy(n,3:y_number-1)));
            %equation at the grid point 'y_number' (with boundary for u: no shear stress at the free surface):
            diag_4_2(y_number-1)=(rho(n,y_number).*((2*U(mod(o,2)+1,y_number)-u(n-1,y_number))/dx(n-1,1)))+((Teta(n,y_number-1)*my(n,y_number))/dy(n,y_number-1));
            %upper diagonal (partial derivatives with respect to u_i+1):
            %equations at the grid point '2' (with boundary: v_2=0):
            diag_4_1(1)=-((Teta(n,1)*my(n,2)/dy(n,2))+((my(n,3)-my(n,2))/(dy(n,2)^2)));
            %equations at the grid points '2 bis y_number-1':         
            diag_4_1(2:y_number-2)=-((Teta(n,2:y_number-2).*my(n,3:y_number-1)./dy(n,3:y_number-1))+((my(n,4:y_number)-my(n,3:y_number-1))./(dy(n,3:y_number-1).^2))-((rho(n,3:y_number-1).*V(mod(o,2)+1,3:y_number-1))./dy(n,3:y_number-1)));   
            %lower secondary diagonal (partial derivatives with respect to u_i-1):
            diag_4_3(1:y_number-2)=-Teta(n,2:y_number-1).*my(n,3:y_number)./(dy(n,2:y_number-1));
            %Creating matrices out of diagonals:
            J_4_2=diag(diag_4_2);
            J_4_1=diag(diag_4_1,1);
            J_4_3=diag(diag_4_3,-1);
            J_4=J_4_1+J_4_2+J_4_3;%Fourth quadrant of J
            
            %The full Jacobi matrix of coefficients contains the 4 quadrants:
            J=[J_1 J_2; J_3 J_4];
            
            %the residua 'delta_y' are the differences between 2 iterations
            %of Newton's method:
            delta_y=sparse(J)\(I-H*uv_vektor);%A\b solves for x of the equation system  Ax=b
            
            %Creating V and U values for the following iterations:
            if(mod(o,2)==0) %query to see if 'o' is even or odd (depending on the result
                %the content of U/V(2,:) or U/V(1,:) will be overwritten)

                %residua of velocities are added:
                %V(2,:)=V(1,:)+1.half of vector of residua delta_y:
                V(mod(o,2)+2,3:y_number+1)=V(mod(o,2)+1,3:y_number+1)+delta_y(1:y_number-1,1)';
                %U(2,:)=U(1,:)+2.half of vector of residua delta_y:
                U(mod(o,2)+2,2:y_number)=U(mod(o,2)+1,2:y_number)+delta_y(y_number:2*y_number-2,1)';
                %Boundary condition at the interface:
                U(mod(o,2)+2,y_number+1)=U(mod(o,2)+2,y_number);
            else
                %residua of velocities are added:
                %V(1,:)=V(2,:)+1.half of vector of residua delta_y:
                V(mod(o,2),3:y_number+1)=V(mod(o,2)+1,3:y_number+1)+delta_y(1:y_number-1,1)';
                %U(1,:)=U(2,:)+2.half of vector of residua delta_y:
                U(mod(o,2),2:y_number)=U(mod(o,2)+1,2:y_number)+delta_y(y_number:2*y_number-2,1)';
                %Boundary condition at the interface:
                U(mod(o,2),y_number+1)=U(mod(o,2),y_number);
            end
            o=o+1;%counter of iteration 
            %break condition of Newton loop for velocities 
            %(residua must be smaller than the values defined in variables 'Genauigkeit...' 
            %break occurs also if loop is repeated more than 1000 times 
            %(=>forced break to avoid infinite loop)):
            if((abs(U(1,:)-U(2,:))<=Genauigkeit_U)&(abs(V(1,:)-V(2,:))<=Genauigkeit_V)|(o>10^3))
                %a warning is creatied if this is the case:
                if(o>10^3)
                    warning('velocity loop did not converge within 1000 repetitions')
                end
                break
            end
        end%end Newton loop for velocity distribution 
        
  %%      
        %the final solutions of the substitute variables are transferred to the matrix of results:
        u(n,:)=U(mod(o,2)+1,:);
        v(n,:)=V(mod(o,2)+1,:);

q=1;%Counter of temperature and mass fraction loop based on Newton's method (starts at 1)
        %t_sub und W_sub are the substitutes of t und W and also 
        %the objectives of the following iterations in the Newton loop
        %
        t_sub=zeros(2,y_number+1);%initialising of t_sub
        W_sub=zeros(2,y_number+1);%initialising of W_sub
        t_sub(1:2,1:y_number+1)=repmat(t(mod(i,2)+1,:),2,1);%Substitute 
        %t_sub gets the resulting values t of the finished iteration before
        W_sub(1:2,1:y_number+1)=repmat(W(mod(i,2)+1,:),2,1);%Substitute 
        %W_sub gets the resulting values W of the finished iteration before
           
        %Start of the Newton loop to solve for the temperature and mass fraction distribution
        %
        while (0==0)%Abbruchbedingung in Z.722
            %The following code has a very similar structure to the loop for  
            %solving the velocity field. The solution of t_sub and 
            %W_sub are saved in a vector called 'Tw_vektor'
	    %As before first the values of W_sub then values for t_sub
            %The system of equations will be saved in matrix notation  
            %in the form F*Tw_vektor=G 
            %G is the inhomogeneity (as I before) and F is the 
            %matrix of coefficients for temperatures and mass fractions (as H before)

	    %The matrix of coefficients comprised of 6 arrays. 
            %However, not all of the same size. The structure 'F' is designed as
            %[Array1 Array2; Array3 Array4; Array5 Array6]. The
            %array 1 and 3 are matrices of the dimension (y_number-1)by (y_number)

            %array 2 and 4 are matrices of the dimension (y_number-1) by (y_number-1)
            %Array 5 is a row vector of the size y_number and
            %array 6 is a row vector of the size  y_number-1 
            %The reason for the different 'sizes' or dimensions is:
            %due to the different boundary conditions, mass fractions are 
	    %calculated from i=2 to i=y_number+1, but temperatures are just 
            %calculated from i=2 to i=y_number. The temperature at the interface is 
            %determined by the function for saturation conditions (assumtion in the model) 
            
           
            %Array 1 in F contains the coefficients of w of the equation 
            %of mass transfer or component balance):
            %Array1 is not a square matrix. Nevertheless, the pricipal diagonal
            %(Diag_1_2) is the one with  
            %concurrent indices for rows and columns. 
            %The secondary diagonals are defined in a way that  
            %upper secondary diagonal (Diag_1_1) is of  the same size as
            %the pricipal diagonal. However, the lower secondary 
            %diagonal (Diag_1_3) is shorter by one entry
            Diag_1_3=zeros(1,y_number-2);%definition of the size of Diag_1_3
            Diag_1_2=zeros(1,y_number-1);%Definition of the size of Diag_1_2
            Diag_1_1=zeros(1,y_number-1);%Definition of the size of Diag_1_1
            %principal diagonal (coefficients of w_i):
            %equation at the grid point 2 (with boundary: impermeabel wall -> no gradient of w):
            Diag_1_2(1)=rho(n,2)*u(n,2)/dx(n-1,1)+Teta(n,1)*phi(n,2)/(dy(n,2));        
            %equations of domain 3 to y_number-1:
            Diag_1_2(2:y_number-2)=(rho(n,3:y_number-1).*((u(n,3:y_number-1)./dx(n-1,1))+(v(n,3:y_number-1)./dy(n,2:y_number-2)))-((phi(n,4:y_number)-phi(n,3:y_number-1))./(dy(n,3:y_number-1).*dy(n,2:y_number-2)))+(phi(n,3:y_number-1).*Teta(n,2:y_number-2)./dy(n,3:y_number-1))+(phi(n,3:y_number-1).*Teta(n,2:y_number-2)./dy(n,2:y_number-2)));
            
            %equation at grid point y_number:
            Diag_1_2(y_number-1)=(rho(n,y_number)*((u(n,y_number)/dx(n-1,1))+(v(n,y_number)/dy(n,y_number-1)))-((phi(n,y_number+1)-(phi(n,y_number)))/(dy(n,y_number)*dy(n,y_number-1)))+(phi(n,y_number)*Teta(n,y_number-1)/dy(n,y_number))+(phi(n,y_number)*Teta(n,y_number-1)/dy(n,y_number-1)));
            
            %upper secondary diagonal (coefficients of w_i+1):
            %equation of domain  2 to y_number-1:
            Diag_1_1(1:y_number-2)=-phi(n,2:y_number-1).*Teta(n,1:y_number-2)./(dy(n,2:y_number-1));
            %equation from at grid point 'y_number':
            Diag_1_1(y_number-1)=-phi(n,y_number).*Teta(n,y_number-1)./(dy(n,y_number));
            %lower secondary diagonal (coefficients of w_i-1):
            Diag_1_3(1:y_number-3)=-((rho(n,3:y_number-1).*v(n,3:y_number-1)./dy(n,2:y_number-2))-((phi(n,4:y_number)-phi(n,3:y_number-1))./(dy(n,3:y_number-1).*dy(n,2:y_number-2)))+(phi(n,3:y_number-1).*Teta(n,2:y_number-2)./dy(n,2:y_number-2)));
            
            %equation at the grid point 'y_number':
            Diag_1_3(y_number-2)=-(rho(n,y_number)*(v(n,y_number)/dy(n,y_number-1))-((phi(n,y_number+1)-(phi(n,y_number)))/(dy(n,y_number)*dy(n,y_number-1)))+(phi(n,y_number)*Teta(n,y_number-1)/dy(n,y_number-1)));
            
            %Subsequently the matrices are created. This time the dimensions  
            %have to be adjusted to attain a square matrix.
            %Thus a column vector is added filled with zeros. 
            F_1_2=[diag(Diag_1_2) zeros(y_number-1,1)];
            F_1_1=[zeros(y_number-1,1) diag(Diag_1_1)];
            F_1_3=[diag(Diag_1_3,-1) zeros(y_number-1,1)]; 
            F_1=F_1_1+F_1_2+F_1_3;%creating array1 of F
            
            %Array 2 of F (coefficients of T in component balance): 
            %since the tempeatur is not part of the component balance
            %the array is filled with zeros
            F_2=zeros(y_number-1);%creating array 2 of F
          
            %Array 3 of matrix F (coefficients of w in energy balance):
            %array 3 is like array 1 not of squared shape.
            Diag_3_1=zeros(1,y_number-1);%definition of the size of Diag_3_1
            Diag_3_2=zeros(1,y_number-1);%definition of the size of Diag_3_2            
            Diag_3_3=zeros(1,y_number-2);%definition of the size of Diag_3_3
            %principal diagonal (coefficients of w_i):
            %equation at grid point '2' (with boundary: impermeabel wall -> no gradient of w):
            Diag_3_2(1)=(rho(n,2)*u(n,2)*psi(n,2)/dx(n-1,1))+(phi(n,2)*psi(n,2)*Teta(n,1)/dy(n,2));
            %equations of domain  3 to y_number-1:
            Diag_3_2(2:y_number-2)=(rho(n,3:y_number-1).*psi(n,3:y_number-1).*((u(n,3:y_number-1)./dx(n-1,1))+(v(n,3:y_number-1)./dy(n,2:y_number-2)))-((phi(n,4:y_number).*psi(n,4:y_number)-(phi(n,3:y_number-1).*psi(n,3:y_number-1)))./(dy(n,3:y_number-1).*dy(n,2:y_number-2)))+(phi(n,3:y_number-1).*psi(n,3:y_number-1).*Teta(n,2:y_number-2)./dy(n,3:y_number-1))+(phi(n,3:y_number-1).*psi(n,3:y_number-1).*Teta(n,2:y_number-2)./dy(n,2:y_number-2)));
            %equation at the grid point 'y_number' :
            Diag_3_2(y_number-1)=(rho(n,y_number)*psi(n,y_number)*((u(n,y_number)/dx(n-1,1))+(v(n,y_number)/dy(n,y_number-1)))-(((phi(n,y_number+1)*psi(n,y_number+1))-(phi(n,y_number)*psi(n,y_number)))/(dy(n,y_number)*dy(n,y_number-1)))+(phi(n,y_number)*psi(n,y_number)*Teta(n,y_number-1)/dy(n,y_number))+(phi(n,y_number)*psi(n,y_number)*Teta(n,y_number-1)/dy(n,y_number-1)));
                   
            %upper secondary diagonal (coefficients of w_i+1):
            %equations from 2 to y_number-1:
            Diag_3_1(1:y_number-2)=-phi(n,2:y_number-1).*psi(n,2:y_number-1).*Teta(n,1:y_number-2)./(dy(n,2:y_number-1));
            %equation at grid point 'y_number' (with saturation condition):
            Diag_3_1(y_number-1)=-(((Teta(n,y_number-1)*lambda(n,y_number)*(T_Siede(1-W_sub(mod(q,2)+1,y_number+1),p)+273.15))/(W_sub(mod(q,2)+1,y_number+1)*(dy(n,y_number))))+(Teta(n,y_number-1)*phi(n,y_number)*psi(n,y_number)/dy(n,y_number)));
            
            %lower secondary diagonal (coefficients of w_i-1):
            %equations of domain  3 to y-Anzahl-1
            Diag_3_3(1:y_number-3)=-((rho(n,3:y_number-1).*v(n,3:y_number-1).*psi(n,3:y_number-1)./dy(n,2:y_number-2))-((phi(n,4:y_number).*psi(n,4:y_number)-(phi(n,3:y_number-1).*psi(n,3:y_number-1)))./(dy(n,3:y_number-1).*dy(n,2:y_number-2)))+(phi(n,3:y_number-1).*psi(n,3:y_number-1).*Teta(n,2:y_number-2)./dy(n,2:y_number-2)));
            %equation at the grid point 'y_number':
            Diag_3_3(y_number-2)=-((rho(n,y_number)*v(n,y_number)*psi(n,y_number)/dy(n,y_number-1))-(((phi(n,y_number+1)*psi(n,y_number+1))-(phi(n,y_number)*psi(n,y_number)))/(dy(n,y_number)*dy(n,y_number-1)))+(phi(n,y_number)*psi(n,y_number)*Teta(n,y_number-1)/dy(n,y_number-1)));

            %Generating the (non-quadratic) matrices from formerly defined diagonals:
            F_3_2=[diag(Diag_3_2) zeros(y_number-1,1)];
            F_3_1=[zeros(y_number-1,1) diag(Diag_3_1)];
            F_3_3=[diag(Diag_3_3,-1) zeros(y_number-1,1)];
            F_3=F_3_1+F_3_2+F_3_3;%creating array 3 of F
            
            
            %Array 4 of F (coefficients of T in the energy balance):
            %Array 4 is a quadratic matrix
            %principal diagonal (coefficients of T_i):
            %equation at grid point '2' (with boundary: v_2=0):
            Diag_4_2(1)=rho(n,2)*cp(n,2)*u(n,2)/dx(n-1,1)-((lambda(n,3)-lambda(n,2))/(dy(n,2)*dy(n,1)))+(Teta(n,1)*lambda(n,2))/(dy(n,2))+(Teta(n,1)*lambda(n,2))/(dy(n,1));
            %equations of domain 3 to y_number:
            Diag_4_2(2:y_number-1)=(rho(n,3:y_number).*cp(n,3:y_number).*((u(n,3:y_number)./dx(n-1,1))+(v(n,3:y_number)./dy(n,2:y_number-1)))-((lambda(n,4:y_number+1)-(lambda(n,3:y_number)))./(dy(n,3:y_number).*dy(n,2:y_number-1)))+(Teta(n,2:y_number-1).*lambda(n,3:y_number)./dy(n,3:y_number))+(Teta(n,2:y_number-1).*lambda(n,3:y_number)./dy(n,2:y_number-1)));
            %upper secondary diagonal (coefficients of T_i+1):
            %equations of domain 2 to y_number-1
            Diag_4_1(1:y_number-2)=-Teta(n,1:y_number-2).*lambda(n,2:y_number-1)./(dy(n,2:y_number-1));
            %untere Nebendiagonale (Koeffizienten von T_i-1):
            %equations of domain 3 to y_number
            Diag_4_3(1:y_number-2)=-((rho(n,3:y_number).*v(n,3:y_number).*cp(n,3:y_number)./dy(n,2:y_number-1))-((lambda(n,4:y_number+1)-(lambda(n,3:y_number)))./(dy(n,3:y_number).*dy(n,2:y_number-1)))+(Teta(n,2:y_number-1).*lambda(n,3:y_number)./dy(n,2:y_number-1)));
            
            %Generating the matrix from the diagonals:
            F_4_2=diag(Diag_4_2);
            F_4_1=diag(Diag_4_1,1);
            F_4_3=diag(Diag_4_3,-1);
            F_4=F_4_1+F_4_2+F_4_3;%generating array 4 of F
            
            %Array 5 of F (coefficients of w in the coupling condition between heat and mass transfer):
            F_5=zeros(1,y_number);%definition of the size of array5
            %coefficient of w_y_number:
            % (with boundary of uni-directional diffusion)
            F_5(y_number-1)=rho(n,y_number+1)*Diff(n,y_number+1)*habs(n,1)/((1-W_sub(mod(q,2)+1,y_number+1))*dy(n,y_number));
            %coefficient of w_y_number+1:
            % % (with boundary of uni-directional diffusion)
            F_5(y_number)=(T_Siede(1-W_sub(mod(q,2)+1,y_number+1),p)+273.15)*lambda(n,y_number+1)/(W_sub(mod(q,2)+1,y_number+1)*dy(n,y_number))-rho(n,y_number+1)*Diff(n,y_number+1)*habs(n,1)/((1-W_sub(mod(q,2)+1,y_number+1))*dy(n,y_number));
            
            %Array 6 of F (coefficients of T in the coupling condition between heat and mass transfer):
            F_6=zeros(1,y_number-1);%definition of the size of array 6
            %coefficient at y_number:
            F_6(y_number-1)=-lambda(n,y_number+1)/dy(n,y_number);
            
            %Generating the matrix of coefficients from the 6 arrays:
            F=[F_1 F_2; F_3 F_4; F_5 F_6];
            
            %the variable 'Tw_vektor' is the vector for saving the intermediate  
            %results of the iteration runs 
            %it contains the values of W_sub_i from i=2 to i=y_number+1 and the 
            %values of t_sub_i from i=2 to i=y_number:
            Tw_vektor=[W_sub(mod(q,2)+1,2:y_number+1)'; t_sub(mod(q,2)+1,2:y_number)'];
            
            %G is the column vector for the inhomogeneity of both the
            %equation for mass transfer and the energy balance from index i=2
            %to i=y_number. Additionally, G contains the inhomogeneity  
            %(=0) of the coupling condition between heat and mass transfer.
            %%
            G=[rho(n,2:y_number)'.*u(n,2:y_number)'.*w(n-1,2:y_number)'./dx(n-1,1); rho(n,2)*u(n,2)*cp(n,2)*T(n-1,2)/dx(n-1,1)+ rho(n,2)*u(n,2)*psi(n,2)*w(n-1,2)/dx(n-1,1)-(((lambda(n,3)-(lambda(n,2)))/(dy(n,2)*dy(n,1)))-(Teta(n,1)*lambda(n,2)/dy(n,1)))*Tw; rho(n,3:y_number)'.*u(n,3:y_number)'.*cp(n,3:y_number)'.*T(n-1,3:y_number)'./dx(n-1,1)+ rho(n,3:y_number)'.*u(n,3:y_number)'.*psi(n,3:y_number)'.*w(n-1,3:y_number)'./dx(n-1,1); 0];
            
            %%          
            %Again the Jacobi matrix is necessary for applying Newton's method on 
            %the energy and component balance
            %The Jacobi matrix is now called R...
            %The structure resembles the matrix for the coefficients
            
            %Array 1 of R (partial derivatives of w in  
            %component balance):
            %principal diagonal contains (partial derivatives with respect to w_i):
            %equation at grid point '2' (with boundary of impermeable wall):
            Diag_1_2(1)=rho(n,2)*u(n,2)/dx(n-1,1)+Teta(n,1)*phi(n,2)/(dy(n,2));
            %equations in the domain '3 to y_number-1':
            Diag_1_2(2:y_number-2)=(rho(n,3:y_number-1).*((u(n,3:y_number-1)./dx(n-1,1))+(v(n,3:y_number-1)./dy(n,2:y_number-2)))-((phi(n,4:y_number)-phi(n,3:y_number-1))./(dy(n,3:y_number-1).*dy(n,2:y_number-2)))+(phi(n,3:y_number-1).*Teta(n,2:y_number-2)./dy(n,3:y_number-1))+(phi(n,3:y_number-1).*Teta(n,2:y_number-2)./dy(n,2:y_number-2)));
            %equation at grid point 'y_number':
            Diag_1_2(y_number-1)=(rho(n,y_number)*((u(n,y_number)/dx(n-1,1))+(v(n,y_number)/dy(n,y_number-1)))-((phi(n,y_number+1)-(phi(n,y_number)))/(dy(n,y_number)*dy(n,y_number-1)))+(phi(n,y_number)*Teta(n,y_number-1)/dy(n,y_number))+(phi(n,y_number)*Teta(n,y_number-1)/dy(n,y_number-1)));
            
            %upper secondary diagonal (partial derivatives with respect to w_i+1):
            %equations in the domain '2 to y_number-1':
            Diag_1_1(1:y_number-2)=-phi(n,2:y_number-1).*Teta(n,1:y_number-2)./(dy(n,2:y_number-1));
            %equation at the grid point 'y_number':
            Diag_1_1(y_number-1)=-phi(n,y_number).*Teta(n,y_number-1)./(dy(n,y_number));
            
            %lower secondary diagonal (partial derivatives with respect to w_i-1):
            Diag_1_3(1:y_number-3)=-((rho(n,3:y_number-1).*v(n,3:y_number-1)./dy(n,2:y_number-2))-((phi(n,4:y_number)-phi(n,3:y_number-1))./(dy(n,3:y_number-1).*dy(n,2:y_number-2)))+(phi(n,3:y_number-1).*Teta(n,2:y_number-2)./dy(n,2:y_number-2)));
            %equation at grid point 'y_number':
            Diag_1_3(y_number-2)=-(rho(n,y_number)*(v(n,y_number)/dy(n,y_number-1))-((phi(n,y_number+1)-(phi(n,y_number)))/(dy(n,y_number)*dy(n,y_number-1)))+(phi(n,y_number)*Teta(n,y_number-1)/dy(n,y_number-1)));
            
            %Generating the matrices from the diagonals:
            R_1_2=[diag(Diag_1_2) zeros(y_number-1,1)];
            R_1_1=[zeros(y_number-1,1) diag(Diag_1_1)];
            R_1_3=[diag(Diag_1_3,-1) zeros(y_number-1,1)];
            R_1=R_1_1+R_1_2+R_1_3;%Creating Array 1 of the Jacobi matrix
                      
            %Array 2 of R (partial derivatives of T in the
            %component balance):
            %Again, it is full of zeros, since the temperature is not present in component balance:
            R_2=zeros(y_number-1);%Creating Array 2 of the Jacobi matrix
            
         
            %Array 3 of R (partial derivatives of w in the 
            %energy balance):
            %Array 3 is like Array 1 a non-square array
            %principal diagonal (partial derivatives with respect to w_i):
            %equations at the grid point '2' (with boundary of impermeable wall):
            Diag_3_2(1)=(rho(n,2)*u(n,2)*psi(n,2)/dx(n-1,1))+(phi(n,2)*psi(n,2)*Teta(n,1)/dy(n,2));
            %Equations in domain 3 to y_number-1:
            Diag_3_2(2:y_number-2)=(rho(n,3:y_number-1).*psi(n,3:y_number-1).*((u(n,3:y_number-1)./dx(n-1,1))+(v(n,3:y_number-1)./dy(n,2:y_number-2)))-((phi(n,4:y_number).*psi(n,4:y_number)-(phi(n,3:y_number-1).*psi(n,3:y_number-1)))./(dy(n,3:y_number-1).*dy(n,2:y_number-2)))+(phi(n,3:y_number-1).*psi(n,3:y_number-1).*Teta(n,2:y_number-2)./dy(n,3:y_number-1))+(phi(n,3:y_number-1).*psi(n,3:y_number-1).*Teta(n,2:y_number-2)./dy(n,2:y_number-2)));
            
            %Equation at the grid point 'y_number':
            Diag_3_2(y_number-1)=(rho(n,y_number)*psi(n,y_number)*((u(n,y_number)/dx(n-1,1))+(v(n,y_number)/dy(n,y_number-1)))-(((phi(n,y_number+1)*psi(n,y_number+1))-(phi(n,y_number)*psi(n,y_number)))/(dy(n,y_number)*dy(n,y_number-1)))+(phi(n,y_number)*psi(n,y_number)*Teta(n,y_number-1)/dy(n,y_number))+(phi(n,y_number)*psi(n,y_number)*Teta(n,y_number-1)/dy(n,y_number-1)));

            %upper secondary diagonal (partial derivatives with respect to w_i+1):
            %Equations in domain '2 to y_number-1':
            Diag_3_1(1:y_number-2)=-phi(n,2:y_number-1).*psi(n,2:y_number-1).*Teta(n,1:y_number-2)./(dy(n,2:y_number-1));
     
            %equation at grid point 'y_number' (with derivative of saturation condition):
            Diag_3_1(y_number-1)=(-(Teta(n,y_number-1)*phi(n,y_number)*psi(n,y_number))-Teta(n,y_number-1)*lambda(n,y_number)*dTSiede_dw_Loewer(1-W_sub(mod(q,2)+1,y_number+1),p))/(dy(n,y_number))-0;

            %lower secondary diagonal (partial derivatives with respect to w_i-1):
            Diag_3_3(1:y_number-3)=-((rho(n,3:y_number-1).*v(n,3:y_number-1).*psi(n,3:y_number-1)./dy(n,2:y_number-2))-((phi(n,4:y_number).*psi(n,4:y_number)-(phi(n,3:y_number-1).*psi(n,3:y_number-1)))./(dy(n,3:y_number-1).*dy(n,2:y_number-2)))+(phi(n,3:y_number-1).*psi(n,3:y_number-1).*Teta(n,2:y_number-2)./dy(n,2:y_number-2)));

            %equation at grid point 'y_number':             
            Diag_3_3(y_number-2)=-(rho(n,y_number)*(v(n,y_number)*psi(n,y_number)/dy(n,y_number-1))-(((phi(n,y_number+1)*psi(n,y_number+1))-(phi(n,y_number)*psi(n,y_number)))/(dy(n,y_number)*dy(n,y_number-1)))+(phi(n,y_number)*psi(n,y_number)*Teta(n,y_number-1)/dy(n,y_number-1)));
            
            %Generating the matrices from the diagonals:
            R_3_2=[diag(Diag_3_2) zeros(y_number-1,1)];
            R_3_1=[zeros(y_number-1,1) diag(Diag_3_1)];
            R_3_3=[diag(Diag_3_3,-1) zeros(y_number-1,1)];
            R_3=R_3_1+R_3_2+R_3_3;%Creating Array 3 of the Jacobi matrix

            %Array 4 of R (partial derivatives of T in the energy balance: 
            %Array 4  is a quadratic matrix
            %principal diagonal (partial derivatives with respect to T_i):
            %equation at the grid point '2' (with boundary: v_2=0):
            Diag_4_2(1)=rho(n,2)*cp(n,2)*u(n,2)/dx(n-1,1)-((lambda(n,3)-lambda(n,2))/(dy(n,2)*dy(n,1)))+(Teta(n,1)*lambda(n,2))/(dy(n,2))+(Teta(n,1)*lambda(n,2))/(dy(n,1));
            %equations in the domain 3 to y_number:
            Diag_4_2(2:y_number-1)=(rho(n,3:y_number).*cp(n,3:y_number).*((u(n,3:y_number)./dx(n-1,1))+(v(n,3:y_number)./dy(n,2:y_number-1)))-((lambda(n,4:y_number+1)-(lambda(n,3:y_number)))./(dy(n,3:y_number).*dy(n,2:y_number-1)))+(Teta(n,2:y_number-1).*lambda(n,3:y_number)./dy(n,3:y_number))+(Teta(n,2:y_number-1).*lambda(n,3:y_number)./dy(n,2:y_number-1)));
            %upper secondary diagonal (partial derivatives with respect to  T_i+1):          
            Diag_4_1(1:y_number-2)=-Teta(n,1:y_number-2).*lambda(n,2:y_number-1)./(dy(n,2:y_number-1));
            %lower secondary diagonal (partial derivatives with respect to  T_i-1):
            %equations in the domain 3 to y_number
            Diag_4_3(1:y_number-2)=-((rho(n,3:y_number).*v(n,3:y_number).*cp(n,3:y_number)./dy(n,2:y_number-1))-((lambda(n,4:y_number+1)-(lambda(n,3:y_number)))./(dy(n,3:y_number).*dy(n,2:y_number-1)))+(Teta(n,2:y_number-1).*lambda(n,3:y_number)./dy(n,2:y_number-1)));
            
            %Generating the matrices from the diagonals:
            R_4_2=diag(Diag_4_2);
            R_4_1=diag(Diag_4_1,1);
            R_4_3=diag(Diag_4_3,-1);
            R_4=R_4_1+R_4_2+R_4_3;%Creating Array 4 of the Jacobi matrix
            
            %Array 5 of R (partial derivatives of w in the  
            %coupling condition between heat and mass transfer:
            R_5=zeros(1,y_number);%definition of the size of array 5
            
            %partial derivatives with respect to w_y_number:
            R_5(y_number-1)=rho(n,y_number+1)*Diff(n,y_number+1)*habs(n,1)/((1-W_sub(mod(q,2)+1,y_number+1))*dy(n,y_number));
                        
            %partial derivatives with respect to w_y_number+1 (under consideration of uni-drectional diffusion):
            R_5(y_number)=(rho(n,y_number+1)*Diff(n,y_number+1)*habs(n,1)*(W_sub(mod(q,2)+1,y_number)-1)+lambda(n,y_number+1)*(W_sub(mod(q,2)+1,y_number)-1)^2*dTSiede_dw_Loewer(1-W_sub(mod(q,2)+1,y_number+1),p))/((W_sub(mod(q,2)+1,y_number)-1)^2*dy(n,y_number));

            %Array 6 of R (partial derivatives of T in 
            %coupling condition between heat and mass transfer:
            R_6=zeros(1,y_number-1);%definition of the size of array 6
            %partial derivatives with respect to T_y_number:
            R_6(y_number-1)=-lambda(n,y_number+1)/dy(n,y_number);
           
            %Generating the full Jacobi matrix R out of the 6 arrays:
            R=[R_1 R_2; R_3 R_4; R_5 R_6];
            
            %the residua 'delta_x' are the differences between 2 iterations
            %of Newton's method: 
            delta_x=sparse(R)\(G-F*Tw_vektor); %similar to delta_y in the velocity field above
            
            %Creating W_sub and t_sub values for the following iterations: 
            if(mod(q,2)==0)%query to see if 'q' is even or odd (depending on the result the first part is computed or the else
		%residua of velocities are added:
                %W_sub(2,:)=W_sub(1,:)+1.half of vector of residua delta_x:
                W_sub(mod(q,2)+2,2:y_number+1)=W_sub(mod(q,2)+1,2:y_number+1)+delta_x(1:y_number,1)';
                %t_sub(2,:)=t_sub(1,:)+2.half of vector of residua delta_x:
                t_sub(mod(q,2)+2,2:y_number)=t_sub(mod(q,2)+1,2:y_number)+delta_x(y_number+1:2*y_number-1,1)';
                %boundary condition at the wall:
                W_sub(mod(q,2)+2,1)=W_sub(mod(q,2)+2,2);
                %saturation condition at the interface
                %(to solve for the temperature there):
                t_sub(mod(q,2)+2,y_number+1)=T_Siede(1-W_sub(mod(q,2)+2,y_number+1),p)+273.15;
            else
                %W_sub(1,:)=W_sub(2,:)+1.half of vector of residua  delta_x:
                W_sub(mod(q,2),2:y_number+1)=W_sub(mod(q,2)+1,2:y_number+1)+delta_x(1:y_number,1)';
                %t_sub(1,:)=t_sub(2,:)+2.half of vector of residua delta_x:
                t_sub(mod(q,2),2:y_number)=t_sub(mod(q,2)+1,2:y_number)+delta_x(y_number+1:2*y_number-1,1)';
                %boundary condition at the wall:
                W_sub(mod(q,2),1)=W_sub(mod(q,2),2);
                %saturation condition at the interface
                %(to solve for the temperature there):
                t_sub(mod(q,2),y_number+1)=T_Siede(1-W_sub(mod(q,2),y_number+1),p)+273.15;
            end

            q=q+1;%counter of iteration 
            %break condition of Newton loop for velocities 
            %(residua must be smaller than the values defined in variables 'Genauigkeit...' 
            %break occurs also if loop is repeated more than 1000 times 
            %(=>forced break to avoid infinite loop)):
             if((abs(t_sub(1,:)-t_sub(2,:))<=Genauigkeit_t_sub)&(abs(W_sub(1,:)-W_sub(2,:))<=Genauigkeit_W_sub)|(q>10^3))
                %a warning is creatied if this is the case:
                if(q>10^3)
                    warning('energy or component balance did not converge within 1000 repetitions')
                end
                break
            end
        end%end Newton loop for temperature and mass fraction distribution        
        
        
        %%the final solutions of the substitute variables t_sub and W_sub are transferred to the matrix of results:
        if(mod(i,2)==0)
        %
        t(mod(i,2)+2,:)=t_sub(mod(q,2)+1,:);
        W(mod(i,2)+2,:)=W_sub(mod(q,2)+1,:);
                
%%
            %Calculation of the film thickness for the next iteration of the coupling loop:
            %
	    %Delta(mod(i,2)+2,1)=Delta(mod(i,2)+1,1); % use this eq. if a constant film
	    %thickness should be considered and comment out the following:
            Delta(mod(i,2)+2,1)=((0.5*(sum(dy(n-1,1:y_number).*((rho(n-1,2:y_number+1).*u(n-1,2:y_number+1))+(rho(n-1,1:y_number).*u(n-1,1:y_number))),2))+((rho(n,y_number+1)*Diff(n,y_number+1)*(W(mod(i,2)+2,y_number+1)-W(mod(i,2)+2,y_number))/((1-W(mod(i,2)+2,y_number+1))*dy(n,y_number))+rho(n-1,y_number+1)*Diff(n-1,y_number+1)*(w(n-1,y_number+1)-w(n-1,y_number))/((1-w(n-1,y_number+1))*dy(n-1,y_number)))/2)*dx(n-1,1))/(0.5*(sum(dy(n-1,1:y_number).*((rho(n,2:y_number+1).*u(n,2:y_number+1))+(rho(n,1:y_number).*u(n,1:y_number))),2))))*delta(n-1,1);


        else 
            t(mod(i,2),:)=t_sub(mod(q,2)+1,:);
            W(mod(i,2),:)=W_sub(mod(q,2)+1,:);
             
            %Calculation of the film thickness for the next iteration of the coupling loop:
	    %Delta(mod(i,2),1)=Delta(mod(i,2)+1,1);  %use this eq. if a constant film
	    %thickness should be considered and comment out the following:
            Delta(mod(i,2),1)=((0.5*(sum(dy(n-1,1:y_number).*((rho(n-1,2:y_number+1).*u(n-1,2:y_number+1))+(rho(n-1,1:y_number).*u(n-1,1:y_number))),2))+((rho(n,y_number+1)*Diff(n,y_number+1)*(W(mod(i,2),y_number+1)-W(mod(i,2),y_number))/((1-W(mod(i,2),y_number+1))*dy(n,y_number))+rho(n-1,y_number+1)*Diff(n-1,y_number+1)*(w(n-1,y_number+1)-w(n-1,y_number))/((1-w(n-1,y_number+1))*dy(n-1,y_number)))/2)*dx(n-1,1))/(0.5*(sum(dy(n-1,1:y_number).*((rho(n,2:y_number+1).*u(n,2:y_number+1))+(rho(n,1:y_number).*u(n,1:y_number))),2))))*delta(n-1,1);

             
        end
        i=i+1;%counter of the coupling loop (coupling between the velocity field and the temperature and mass fraction distribution


        %The break conditions are the differences between 
        %temperatures, mass fractions and film thickness values of two  
        %subsequent iterations of the coupling loop. 
        %If these are smaller than the previously defined values for the variable 'Genauigkeit_...'  
        %the loop stops
        if ((abs(W(1,:)-W(2,:)) <= Genauigkeit_W) & (abs(t(1,:)-t(2,:)) <= Genauigkeit_t) & (abs(Delta(1,1)-Delta(2,1)) <= Genauigkeit_Delta) | (i > 10^2))
            %a warning is creatied if this is the case:
            if (i > 10^3)
                warning('AcctError:NoClient','coupling loop did not converge')
            end
            break
        end
    end%End of coupling loop 
    
    %The final values for temperature and mass fraction distribution are saved in the matrix for results:
    T(n,:)=t(mod(i,2)+1,:);
    w(n,:)=W(mod(i,2)+1,:);
    %The determined film thickness is put into the vector 'delta'
    delta(n,1)=Delta(mod(i,2)+1,1);
    if(check_2==0)
        N_2=n;%N_2 indicates the longitudinal index where the 
        %hydrodynamic boundary layer meets the free surface
    end
    
    if (delta_99_index>y_number-1)
         check_2=1;%check_2 becomes 1 umgeschaltet, as soon the 
         %boundary layer thickness is one grid point away from the free surface
         
     end
    %Output of current position and progress in command window:
    display(['Step# ',num2str(n),'--- longitudinal step size dx:',num2str(dx(n-1,1)),'  ---- current x-coordinate: ',num2str(A_3(n))]); 
    n=n+1;%counter of loop over flow lenght L
end
