classdef Handles < handle & Appearance
    properties (SetObservable,AbortSet)
        referenceAx %axes whose limits will be adjusted
    end
    properties
        Tracker  = struct('Parent',[],'Color',[0.85 0.85 0.85]);%background of slider (initially a structure but will become axes handle)
        Patch  = struct('Parent',[],'FaceColor',[0.4 0.4 0.4])%bar on slider (initially a struct but will become patch handle)
    end
    properties (Dependent)
        TrackerParent
    end
    methods
        
        function propval = getOldprops(obj,TrackerOrPatch)
            %returns cell array of property-value pairs for tracker/patch
            TrackerOrPatch=validatestring(TrackerOrPatch,{'Tracker','Patch'});
            propsrc= obj.(TrackerOrPatch); clear TrackerOrPatch
            %make sure propsrc (which could be handle or struct) is struct
            propsrc=struct(propsrc);
            props = fieldnames(propsrc);
            
            val = struct2cell(propsrc);
            propval = cell(2*numel(props),1);
            propval(1:2:end-1) = props;
            propval(2:2:end) = val; 
        end
        
        function makeTracker(obj)
            %obj.Tracker may be an axes or a struct with axes properties.
            % determine tracker parent
            if isempty(obj.Tracker.Parent) || ~ishandle(obj.Tracker.Parent) ||  ~isvalid(obj.Tracker.Parent)
                %make new figure to put tracker on.
                trackerParent = figure;
            else
                trackerParent = obj.Tracker.Parent;
            end
            
            %retain preexisting properties
            oldProps = getOldprops(obj,'Tracker');
            
            if ishandle(obj.Tracker)
                delete(obj.Tracker);
            end
            
            %make tracker axes
            obj.Tracker=axes(oldProps{:},...%old property settings will be overriden if they are the same properties set below
                'parent',trackerParent,...
                'nextplot','add',...
                'xcolor','none','ycolor','none',...
                'xlimmode','manual','ylimmode','manual',...
                'xlim',[1 2],'ylim',[1 2],...
                'xtickmode','manual','ytickmode','manual',...
                'xtick',[],'ytick',[]);
        end
        
        function makePatch(obj)
            %obj.Patch may be a patch() or a struct with patch properties.
            %retain preexisting properties
            oldProps = getOldprops(obj,'Patch');
            %delete old patch
            if ishandle(obj.Patch)
                delete(obj.Patch);
            end
            if isempty(obj.Tracker) || ~isvalid(obj.Tracker)
                make_tracker(obj); %make patch parent if necessary
            end
            obj.Patch = patch(...
                oldProps{:},...  %old property settings will be overriden if they are the same properties set below
                'xdata',axisinfo.X.defaultdata,...
                'ydata',axisinfo.Y.defaultdata,...
                'parent',obj.Tracker,...
                'edgecolor','none');
        end
        
        function set.TrackerParent(obj,newParent)
            %delete old parent after moving tracker to new parent;
            oldparent = obj.Tracker.Parent;
            if ~isequal(oldparent,newParent)
                obj.Tracker.Parent = newParent;
                delete(oldparent); %important not to delete oldparent if it is the same as the newparent
            end
        end
        
        function p = get.TrackerParent(obj)
            p=obj.Tracker.Parent;
        end
        
        function set.Tracker(obj,newval)
            msg='Tracker may be an axes or a struct with axes properties';
            assert(ishandle(newval) | isstruct(newval),msg);
            obj.Tracker=newval;
            %put Tracker at very bottom of Parent;
            if ishandle(obj.Tracker)
                uistack(obj.Tracker,'bottom');
            end
        end
        
        function set.Patch(obj,newval)
            msg='Patch may be an patch or a struct with patch properties';
            assert(ishandle(newval) | isstruct(newval),msg);
            obj.Patch=newval;
        end
        
        
        
    end
end