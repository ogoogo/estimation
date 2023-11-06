function path_setting(situation)

if situation == "mac"
addpath('~/issl/research/simulator/mice/src/mice/')
addpath('~/issl/research/simulator/mice/lib/')
else
addpath('~/issl/research/simulator/mice/src/mice/')
addpath('~/issl/research/simulator/mice/lib/')
end

if cspice_ktotal( 'ALL' ) >= 1
else
    disp('kernels settings')
    cspice_furnsh('../kernel/naif0012.tls');
    cspice_furnsh('../kernel/de440.bsp');
end

end