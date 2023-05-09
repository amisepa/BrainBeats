classdef Appearance < handle 
    %position and color of scrollbar
    properties (SetObservable,AbortSet)
        thick = 20 %scrollbar will be 20 pixels thick 
        bufferspace = 20; %space between tracker and referenceAx
    end
    properties (Abstract)
        Tracker
        Patch
        referenceAx
        listen2
    end
    properties
        left_right_top_bottom = 'bottom'
    end
    

    properties (Dependent)
        patchColor = [0.4 0.4 0.4]; %dark gray
        trackerColor = [0.8 0.8 0.8]; %light gray
    end
    methods
        %% Color
        
        function c = get.patchColor(obj)
            c = obj.Patch.FaceColor;
        end
        
        function c = get.trackerColor(obj)
            c = obj.Tracker.Color;
        end
        
        function set.patchColor(obj,newval)
            obj.Patch.FaceColor = newval;
        end
        
        function set.trackerColor(obj,newval)
            obj.Tracker.Color = newval;
        end
        
        %example of changing color.
        function yellow(obj)
            obj.patchColor =   [0.85, 0.85, 0.4];
            obj.trackerColor = [1.0, 1.0, 0.7];
        end
        

        
        
        %% Position
        function set.left_right_top_bottom(obj,val)
            obj.left_right_top_bottom = validatestring(val,{'left','right','top','bottom'});
        end
        
        function placeScrollbar(obj)
            side=obj.left_right_top_bottom;
            orguni=obj.referenceAx.Units; 
            obj.referenceAx.Units='pixel';
            obj.Tracker.Units='pixel';
            refpos = obj.referenceAx.Position;
            trackerpos=refpos;
           
            switch side
                case 'right' %scroller is placed right of reference axes
                    trackerpos(1) = refpos(1) + refpos(3) + obj.bufferspace;
                    trackerpos(3) = obj.thick;
                    
                case 'left'
                    trackerpos(1) = refpos(1) - obj.thick - obj.bufferspace;
                    trackerpos(3) = obj.thick;
                case 'top'
                    trackerpos(2) = refpos(2)+refpos(4)+obj.bufferspace;
                    trackerpos(4) = obj.thick;
                case 'bottom'
                    trackerpos(2) = refpos(2)- obj.thick - obj.bufferspace;
                    trackerpos(4) = obj.thick;
            end
            
            obj.Tracker.Position = trackerpos;
            set(obj.referenceAx,'units',orguni); clear orguni
        end
        

    end

end