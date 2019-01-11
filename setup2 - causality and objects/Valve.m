classdef Valve < matlab.mixin.Copyable
   properties
      Kv
      tau
      DmState
      DDm
      record
      hOutlet
   end
   methods
      function flow(valve,capacityRatio,dInlet,pInlet,pOutlet)
         % Kv value is in non-SI units!
         Dm = 8.7841e-06*capacityRatio*valve.Kv*sqrt(dInlet*(pInlet-pOutlet));
         valve.DDm = (Dm - valve.DmState)/valve.tau;
      end
      function enthalpy(valve,hInlet)
          valve.hOutlet = hInlet;
      end
      function DDm = process(valve,t,x,Inputs)
          % Inputs
          dInlet = Inputs.dInlet;
          pInlet = Inputs.pInlet;
          pOutlet = Inputs.pOutlet;
          hInlet = Inputs.hInlet;
          capacityRatio = Inputs.capacityRatio;
          % State
          valve.DmState = x;
          % Process
          valve.flow(capacityRatio,dInlet,pInlet,pOutlet);
          valve.enthalpy(hInlet);
          % Derivatives
          DDm = valve.DDm;
      end
      function initialize(valve,Kv,Tau,Initial)
          valve.Kv = Kv;
          valve.tau = Tau;
          valve.DmState = Initial;
          % Record
          Record.t = [];
          Record.x = [];
          valve.record = Record;
      end
      function reinitialize(valve,Initial)
          valve.DmState = Initial;
          valve.record.t = [];
          valve.record.x = [];
      end
   end
end