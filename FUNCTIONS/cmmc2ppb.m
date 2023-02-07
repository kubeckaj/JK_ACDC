function y = cmmc2ppb(c,temp) 
  % Avogadro's constant
  Na=6.02214179e23;
  % Boltzmann constant [J/K]
  kB=1.3806504e-23;
  R  = 8.314;
  
  p=101325;
  
  y = c*10^9/p*R*temp*10^6/Na
  
end