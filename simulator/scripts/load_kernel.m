addpath('~/simulator/mice/src/mice/')
addpath('~/simulator/mice/lib/')

currentDir = pwd;
cd ../kernel/

% Leapsecond kernel;
cspice_furnsh('naif0012.tls');

% Planet ephemeris kernel;
% cspice_furnsh('de430.bsp');
cspice_furnsh('de440.bsp');
cspice_furnsh('mar097.bsp');
cspice_furnsh('jup310.bsp');
cspice_furnsh('sat360.bsp');
cspice_furnsh('ura090_1.bsp');
cspice_furnsh('nep081.bsp');

% Gravity constant kernel;
cspice_furnsh('gm_de431.tpc');

% Planet constant kernel;
cspice_furnsh('pck00010.tpc');
cspice_furnsh('earth_fixed.tf');
% https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/
%cspice_furnsh('earth_latest_high_prec.bpc'); 

cd ../scripts/



