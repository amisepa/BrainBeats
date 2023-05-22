classdef axisinfo
    %stores information specific to each axis.
    properties
        increaseOn = 'rightarrow'; 
        decreaseOn = 'leftarrow';
        defaultdata =[1 1 2 2 1];
        shiftSide = [1 0]; %left/top
        commandSide = [0 1]; %right/bottom
        referenceAxRepositionInd  = [0 1 0 -1];
    end
    methods
        function obj = axisinfo(increaseOnKey,...
                decreaseOnKey,...
                defaultData,...
                shiftsideInd,...
                commandsideInd,...
                referenceAxRepositionInd)
            obj.increaseOn=increaseOnKey;
            obj.decreaseOn = decreaseOnKey;
            obj.defaultdata = defaultData;
            obj.shiftSide = logical(shiftsideInd);
            obj.commandSide = logical(commandsideInd);
            obj.referenceAxRepositionInd=referenceAxRepositionInd;
            if strcmp(char(obj),'Z')
                warning('this has not been prepared for Z axis scrolling yet')
            end
        end %construction
        %% translating axes limits to patch data and vise versa
        function newData = axesLim_to_PatchData(obj,patchhandle,axeshandle)
            lims = getaxisLimits(obj,axeshandle);
            newData = lims(obj.defaultdata); %patch X/Y data
            setPatchData(obj,patchhandle,newData);
        end
        %% reference plot data determines tracker axes limits
        function trackerAxLim(obj,Lim,axeshandle)
            Lim=[nanmin(Lim),nanmax(Lim)];
            assert(~isempty(Lim),'plothandle data must contain at least one non-NaN element.');
            assert(any(~isnan(Lim)),'plothandle data must contain at least one non-NaN element.');
            errorIfBadLim=false;
            setaxisLimits(obj,axeshandle,Lim,errorIfBadLim)
        end
        
        function setPatchData(obj,patchhandle,newData)
            patchhandle.([char(obj),'Data'])=newData;
        end
        
        function setaxisLimits(obj,axeshandle,newLim,errorIfBadLim)
            try
                axeshandle.([char(obj),'Lim']) = newLim;
            catch ME
                ermsg=['Value must be a 1x2 vector of numeric type in which the',...
                    'second element is larger than the first and may be Inf'];
                if errorIfBadLim || isempty(strfind(ME.message,ermsg))
                    throw(ME)
                end
            end
            
        end
        
        function d = getPlotData(obj,gobjecthandle)
            propname=[char(obj),'Data'];
            assert(ismember(propname,fieldnames(gobjecthandle)),'gobjecthandle must have %s property',propname);
            d = gobjecthandle.(propname);
        end
        
        function l = getaxisLimits(obj,axeshandle)
            l = axeshandle.([char(obj),'Lim']);
        end
        

        %% adjust limits
        function wasmatched = moveObject(obj,axeshandle,KeyPressData,moveby)
            %determine whether keypress applies and detemine whether moveby
            %will increase or decrease axes limits
            
            switch KeyPressData.Key
                case obj.increaseOn;
                    coef = 1;
                case obj.decreaseOn;
                    coef = -1;
                otherwise
                    %key press does not apply
                    wasmatched=false;return;
            end
            wasmatched=true;
            % how much do you change the axes?
            
            moveby = coef*moveby; %apply coefficient.  
            clear movebyStruct coef
            
            %determine which side of the axes to adjust (determine if panning or zooming)
            modls={'shift','command'};
            modii = ismember(modls,KeyPressData.Modifier);
            if ~xor(modii(1),modii(2)) %perform panning
                ind = obj.shiftSide | obj.commandSide;
       %else performing zooming....
            elseif modii(1) %if shift 
                ind = obj.shiftSide;
            elseif modii(2)
                ind = obj.commandSide;
            end; clear modls modii
            
            %get new axes limits 
            
            Lim =  getaxisLimits(obj,axeshandle);
            Lim(ind) = Lim(ind)+moveby; clear moveby ind
            
            %decide to throw an error message if axes limits are bad 
            %regardless of whether or not you throw an error, nothing will be changed.
            
            errorIfBadLim = false;
            
            %don't set patch position. Set axes position (patch position will follow);
            setaxisLimits(obj,axeshandle,Lim,errorIfBadLim);
            
        end

    end
    
    enumeration
        %increaseOn    decreaseON  defaultdata shift side, command side
      X ('rightarrow','leftarrow', [1 1 2 2 1], [1 0] ,[0 1], [0 1 0 -1]) %puts scroller below reference ax
      Y ('uparrow',   'downarrow', [1 2 2 1 1], [0 1], [1 0], [0 0 -1 0]) %puts to left of reference ax
      Z ('to be defined','to be defined',[],false(1,2),false(1,2),[0 0 0 0]);
    end
    
    
end

