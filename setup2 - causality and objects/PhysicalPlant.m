classdef PhysicalPlant < matlab.mixin.Copyable
    properties
        parts
        ODEoptions
        process
        postProcess
        t
        x
    end
    methods
        function timestep(pp,t,inputs)
            % Function help:
            
            [T, X] = ode15s(@pp.process,[t(1) t(2)],pp.x,pp.ODEoptions,inputs,pp);
            pp.t = [pp.t; T];
            pp.x = X(end,:)';
            % Evaluation
            pp.postProcess(X,pp);
        end
        function initialize(pp,Parts,Process,PostProcess,Initial,ODEoptions)
            % Function help: Process and PostProcess are functions, while
            %   parts is a structure made up of objects. The objects are
            %   supposed to be already initialized.
            
            pp.process = Process;
            pp.parts = Parts;
            pp.postProcess = PostProcess;
            pp.ODEoptions = ODEoptions;
            pp.t = [];
            pp.x = Initial;
        end
    end
end