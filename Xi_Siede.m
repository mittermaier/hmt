function aus=Xi_Siede(T,p)
% T			temperature         		C (celcius)
%                                                	
% p                    	pressure	         	Pa 
%                                                	
% Xi_Siede		saturation mass fraction	kg LiBr/kg solution 
%
% This function calculates the saturation mass fraction (mass fraction at thermodynamic equilibrium) of aqueous LiBr
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

epsilon=1e-2;
Startwert=0.01; %initial guess value; 0.0 not defined!
x0=Startwert;
x1=Startwert;

if T_Siede(Startwert,p)-T > 0
   while x0==Startwert
       if T_Siede(x1+epsilon,p)-T > 0
           x1=x1+epsilon;
       end
       if T_Siede(x1+epsilon,p)-T < 0
           x0=x1+epsilon;
       end
   end
end

if T_Siede(Startwert,p)-T < 0
    while x1==Startwert
        if T_Siede(x0+epsilon,p)-T < 0
            x0=x0+epsilon;
        end
        if T_Siede(x0+epsilon,p)-T > 0
            x1=x0+epsilon;
        end
    end
end

xalt=x0;
xaltalt=x1;
xneu=x1; 
durchlauf=0;
while abs(xneu-xalt) > 1e-7
   durchlauf=durchlauf+1;
   if durchlauf > 1
       xaltalt=xalt;
       xalt=xneu;
   end
   if abs((T_Siede(xalt,p)-T)-(T_Siede(xaltalt,p)-T))>0
       xneu=xalt-(xalt-xaltalt)/((T_Siede(xalt,p)-T)-(T_Siede(xaltalt,p)-T))*(T_Siede(xalt,p)-T);
   else
       break;
   end
end
aus=xneu;
