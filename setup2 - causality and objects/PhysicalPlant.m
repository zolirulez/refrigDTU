classdef PhysicalPlant < matlab.mixin.Copyable
    properties
        parts
        ODEoptions
        process
        preProcess
        postProcess
        t
        x
        Inputs
    end
    methods
        function timestep(pp,t)
            % Function help:
            
            [T, X] = ode15s(@pp.process,[t(1) t(2)],pp.x,pp.ODEoptions,pp);
            pp.t = [pp.t; T];
            pp.x = X(end,:)';
            % Evaluation
            pp.postProcess(X,pp);
        end
        function initialInputs(pp,Inputs)
            pp.Inputs = Inputs;
        end
        function initialize(pp,Parts,Process,PostProcess,Initial,ODEoptions)
            % Function help: Process and PostProcess are functions, while
            %   parts is a structure made up of objects. The objects are
            %   supposed to be already initialized.
            
            pp.process = Process;
            pp.parts = Parts;
            pp.postProcess = PostProcess;
            pp.ODEoptions = ODEoptions;
            %pp.ODEoptions.BDF = 'off';
            pp.t = [];
            pp.x = Initial;
        end
    end
end