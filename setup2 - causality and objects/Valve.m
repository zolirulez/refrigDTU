classdef Valve < matlab.mixin.Copyable
   properties
      Kv
      tau
      Dm
      DDm
      record
      h
      pInlet
      pOutlet
      hInlet
      dInlet
   end
   methods
      function flow(valve,capacityRatio,dInlet,pInlet,pOutlet)
         % Kv value is in non-SI units!
         deltap = max([0,pInlet-pOutlet]);
         Dm = 8.7841e-06*capacityRatio*valve.Kv*sqrt(dInlet*(deltap)) ;
         valve.DDm = (Dm - valve.Dm)/valve.tau;
      end
      function enthalpy(valve,hInlet)
          valve.h = hInlet;
      end
      function DDm = process(valve,t,x,capacityRatio)
          % Inputs
%           dInlet = Inputs.dInlet;
%           pInlet = Inputs.pInlet;
%           pOutlet = Inputs.pOutlet;
%           hInlet = Inputs.hInlet;
%           capacityRatio = Inputs.capacityRatio;
          % State
          valve.Dm = x;
          % Process
          valve.flow(capacityRatio,valve.dInlet,valve.pInlet,valve.pOutlet);
          valve.enthalpy(valve.hInlet);
          % Derivatives
          DDm = valve.DDm;
      end
      function initialize(valve,Kv,Tau,Initial)
          valve.Kv = Kv;
          valve.tau = Tau;
          valve.Dm = Initial;
          % Record
          Record.t = [];
          Record.x = [];
          valve.record = Record;
      end
      function reinitialize(valve,Initial)
          valve.Dm = Initial;
          valve.record.t = [];
          valve.record.x = [];
      end
   end
end