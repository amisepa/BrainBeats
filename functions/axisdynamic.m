classdef axisdynamic < handle & Handles
    %malleable properties specific to each axis.
    properties (SetObservable,AbortSet,SetAccess=private)
        private_referenceLim = [0,1];
        enable = true;
    end
    properties
        name = 'X';
        moveSpeed = 2 %percentage of total data (difference of reference limits) that you move by
        default_moveSpeed = 2; %user may wish to go back to this setting after modifying moveSpeed
        listen2=struct('enable',[],'referenceLim',[],'referenceAx',struct('Handle',[],'Lim',[],'Position',[]));
    end
    properties (Dependent)
        axisinfoObj
        referenceLim%limits of scrollbar
    end

    methods
        
        function obj = axisdynamic(axesName,DefaultaxesMoveSpeed,refAx)
            if nargin ~=0
                assert(iscellstr(axesName),'axes names must be cellstr');
                N = numel(axesName); assert(numel(DefaultaxesMoveSpeed)==N,'axesNames,axesMoveSpeed must be same length');
                obj(N,1) = axisdynamic;
                for n = 1:N
                    obj(n,1).name = axesName{n};
                    %if Y-axis then put on right side 
                    %(default position is below reference axes)
                    if strcmp(axesName{n},'Y')
                        obj(n,1).left_right_top_bottom = 'right'; 
                    end
                    obj(n,1).default_moveSpeed = DefaultaxesMoveSpeed(n);
                    obj(n,1).referenceAx= refAx;
                    % listeners
                    Listen2ObjProp(obj(n,1));
                    Listen2AxProp(obj(n,1))
                end
                [obj(:).moveSpeed]=deal(obj(:).default_moveSpeed);
                
            end
            
        end
        %% Listeners
        function Listen2ObjProp(obj)
            %these listeners listen to this object properties are made during construction and don't have to be updated.
            delete([...
                obj.listen2.enable,...
                obj.listen2.referenceLim,...
                obj.listen2.referenceAx.Handle]);
            %enable/disable listeners when to enable property of this object changes
            obj.listen2.enable = addlistener(obj,'enable','PostSet',@(src,ev)onEnable(ev.AffectedObject));
            
            %change X/Y limits of tracker (background of patch) to match reference limits
            obj.listen2.referenceLim = addlistener(obj,'private_referenceLim','PostSet',...
                @(src,ev)set(obj.Tracker,[obj.name,'Lim'],ev.AffectedObject.referenceLim));
            
            %Update listeners etc. when reference Axes handle changes
            obj.listen2.referenceAx.Handle = addlistener(obj,'referenceAx','PostSet', @(src,ev)new_referenceAx(ev.AffectedObject));
        end
        function Listen2AxProp(obj)
            %listen to reference axes properties (must update each time there is new axes)
            delete([...
                obj.listen2.referenceAx.Lim,...
                obj.listen2.referenceAx.Position])
            
            %change tracker position when reference axes moves.
            obj.listen2.referenceAx.Position = addlistener(...
                obj.referenceAx,'Position','PostSet',...
                @(src,ev)placeScrollbar(obj));
            
            %change patch X/YData to match reference Axes' X/Ylimits
            obj.listen2.referenceAx.Lim = addlistener(...
                obj.referenceAx,[obj.name,'Lim'],'PostSet',...
                @(src,ev)axesLim_to_PatchData(obj.axisinfoObj,obj.Patch,ev.AffectedObject));
        end
        
        %% listener callbacks
        function onEnable(obj) % triggered by obj.listen2.enable
            %decide whether listeners are on or off.
            listenls={obj.listen2.referenceAx.Lim, obj.listen2.referenceLim};
            goodii = ~cellfun(@isempty,listenls);
            if ~any(goodii); return;end
            listenls=listenls(goodii);
            set([listenls{:}],'enabled',obj.enable);
        end
        function new_referenceAx(obj) % triggered by obj.listen2.referenceAx.Handle 
            Listen2AxProp(obj);
        end
        
        %% set/get
        function set.name(obj,newval)
            obj.name=validatestring(newval,{'X','Y','Z'});
        end
        function aio = get.axisinfoObj(obj)
            aio = axisinfo.(obj.name);
        end
        %% Reference limits deteremine scrollbar limits
        function set.referenceLim(obj,VecOrHandle)
            if isnumeric(VecOrHandle)
                Vec = VecOrHandle;
            elseif ishandle(VecOrHandle)
                assert(isvalid(VecOrHandle),'handle must be valid');
                Vec = getPlotData(obj.axisinfoObj,VecOrHandle);
            else
                error('class %s invalid setting for referenceLim property',class(VecOrHandle))
            end
            assert(~isempty(Vec),'referenceLim cannot be empty');
            Lim=[nanmin(Vec),nanmax(Vec)];
            assert(~any(isnan(Lim)),'referenceLim cannot contain NaN')
            %set limit
            obj.private_referenceLim = Lim;
        end
        function rl = get.referenceLim(obj)
            rl = obj.private_referenceLim;
        end
        %% make scrollbar
        function buildScroller(obj)
            makeTracker(obj);
            makePatch(obj);
        end
        %% tracker
        function makeTracker(obj)
            makeTracker@Handles(obj)
            %set Tracker limits now
            obj.Tracker.([obj.name,'Lim']) = obj.referenceLim;
            %% ref ax pos
        end
        
        function putTrackerBehindReference(obj)
             if ~ishandle(obj.Tracker) || ~isvalid(obj.Tracker)
                 makeTracker(obj)
             end
            %transfer tracker to reference axes' parent
            obj.TrackerParent = obj.referenceAx.Parent; %tracker axes will be placed at very bottom.
            
            %shrink_referenceAx_to_fit_tracker(obj)
            placeScrollbar(obj)
            
            %send tracker to back of parent (so it is behind reference axes).
            uistack(obj.Tracker,'bottom');
            obj.listen2.referenceAx.Position.Enabled=true;
        end
        
        function shrink_referenceAx_to_fit_tracker(obj)
            orguni=obj.referenceAx.Units; obj.referenceAx.Units='pixel';
            set(obj.Tracker,'units','pixel',...
                'position',obj.referenceAx.Position);
            %index parts of the reference axes position vector
            %(i.e.[left,bottom,width,height]) that will change 
            %and the direction(increase/decrease) they will change in.
            refposind= obj.axisinfoObj.referenceAxRepositionInd;
            %thick and bufferspace are properties of Appearance class
            obj.referenceAx.Position = obj.referenceAx.Position + refposind*obj.thick + refposind*obj.bufferspace;
            set([obj.referenceAx,obj.Tracker],'units',orguni);
        end
        
        %% patch
        function makePatch(obj)
            %Parent is axes is the Tracker (slider background)
            %make slider limits depend on reference limits
            makePatch@Handles(obj);
            
            if ~obj.enable; return; end
            
            set(obj.referenceAx,[obj.name,'limmode'],'manual');
            
            stretchPatch_to_cover_opposite_ax(obj);
            
            %Make Patch limits listen to change in X/Y Lim of reference Ax
            
            %set patch limits to reference ax limits now.
            axesLim_to_PatchData(obj.axisinfoObj,obj.Patch,obj.referenceAx);
        end
        
        function stretchPatch_to_cover_opposite_ax(obj)
            switch obj.name
                case 'X'
                    otherax=axisinfo.Y;
                case 'Y'
                    otherax=axisinfo.X;
                case 'Z'
                    error('not prepared to deal with Z axis yet')
            end
            %if scoller is for X data, the Y limits of the tracker and patch will be [1 2]
            %if scoller is for Y data, the X limits of the tracker and patch will be [1 2]
            trackerAxLim(otherax,otherax.defaultdata,obj.Tracker);  
            setPatchData(otherax,obj.Patch,otherax.defaultdata);
        end
        %% move patch       
        function wasmatched = keypress(obj,KeyPressData)
            % determine if and how to move patch on arrow key press
            dataspan = max(obj.referenceLim)-min(obj.referenceLim);
            moveby =dataspan*obj.moveSpeed/100;
            assert(numel(moveby)==1,'moveby must be scalar');
            wasmatched = moveObject(obj.axisinfoObj,obj.referenceAx,KeyPressData,moveby);
        end
        
        %% move speed
        function adjust_moveSpeed(obj)
            % user can adjust move speed
            promptstr = sprintf('%s-axis speed (percentage of %s data the % s axis limits change by:)',char(obj.name),obj.name,obj.name);
            defval = {num2str(obj.moveSpeed)};
            Answer = inputdlg(promptstr,'Adjust pan/zoom speed.',1,defval);
            newval = str2double(Answer{1});
            if isnan(newval);
                uiwait(errordlg(sprintf('Move speed percentage must be numeric. %s axis speed did not change',char(obj(n,1).name)))); return;
            end
            obj.moveSpeed=newval;
        end
        
        function make_moveSpeed_adjustable(obj)
            arrayfun(@(x)...
                set(x.Patch,'buttondownfcn',@(src,ev)adjust_moveSpeed(x)),...
                obj);
        end

       
        
    end
    
end
