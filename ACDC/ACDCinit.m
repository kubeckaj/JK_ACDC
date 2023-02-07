function [] = ACDCinit(Temp)
%ACDCINIT Summary of this function goes here
%   Detailed explanation goes here

  command="perl ../ACDC/acdc_2020_07_20.pl";
  command=command+" --temperature "+num2str(Temp);
  command=command+" --e ../ACDC_input/HS298.15K.txt";
  command=command+" --i ../ACDC_input/input.inp";
  command=command+" --keep_useless_collisions";
  %command=command+" --disable_nonmonomers";
  command=command+" --dip ../ACDC_input/dip_pol_298.15K.txt";
  %command=command+" --no_generic_ions";
  command=command+" --use_cs --cs exp_loss --exp_loss_coefficient 0.001 --exp_loss_exponent -1.6 ";
  %command=command+" --use_cs --cs exp_loss --exp_loss_coefficient "+num2str(1e-3)+" --exp_loss_exponent "+num2str(-1.6);
  %commadn=command+" --use_wl --wl CLOUD4_simple --use_dilution"
  ok=system(command);
  if ok ~= 0,cd ..;error('Trouble with perl script!');end

end

%command="perl ../ACDC/acdc_2019_08_14.pl"; DOES NOT WORK YET FOR ME
%command="perl ../ACDC/acdc_2017_09_04.pl";
%command="perl ../ACDC/create_acdc_2016_09_30.pl";
