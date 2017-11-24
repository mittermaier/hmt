function enthalpyreg2=enthalpyreg2(temp, press)
% '
% ' specific enthalpy in region 2
% ' enthalpyreg2 in kJ/kg
% ' temperature in K
% ' pressure in bar

tau = 540 ./ temp;
pic = 0.1 .* press;
rgas_water = 461.526;
enthalpyreg2 = 0.001 .* rgas_water .* temp .* tau .* (gamma0taureg2(tau, pic) + gammartaureg2(tau, pic));


