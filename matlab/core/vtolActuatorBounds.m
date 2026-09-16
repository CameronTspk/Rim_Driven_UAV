function [lowerBound, upperBound, rateBound, names] = vtolActuatorBounds(p, architecture)
%VTOLACTUATORBOUNDS Physical position and slew limits in actuator order.
%
% Twin: [T1 alpha1 beta1 T2 alpha2 beta2].'
% Tail: [T1 alpha1 beta1 T2 alpha2 beta2 T3 delta].'

key = vtolArchitectureKey(architecture);

lowerBound = [p.main.Tmin; -p.main.alphaMax; -p.main.betaMax; ...
              p.main.Tmin; -p.main.alphaMax; -p.main.betaMax];
upperBound = [p.main.Tmax;  p.main.alphaMax;  p.main.betaMax; ...
              p.main.Tmax;  p.main.alphaMax;  p.main.betaMax];
rateBound = [p.main.TdotMax; p.main.alphaDotMax; p.main.betaDotMax; ...
             p.main.TdotMax; p.main.alphaDotMax; p.main.betaDotMax];
names = {'T1','alpha1','beta1','T2','alpha2','beta2'};

if strcmp(key, 'tail')
    lowerBound = [lowerBound; p.tail.Tmin; -p.tail.deltaMax];
    upperBound = [upperBound; p.tail.Tmax;  p.tail.deltaMax];
    rateBound = [rateBound; p.tail.TdotMax; p.tail.deltaDotMax];
    names = [names, {'T3','delta'}];
end
end
