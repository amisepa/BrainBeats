classdef superscroll < handle
    %User can pan and zoom a plot (or other graphics object with
    %YData/XData) using arrow keys, shift and command. efzscroll works
    %by adjusting the limits of the axes parent the plot (or other
    %object).You can also adjust scroll/pan speed. 
    
    %This program will adapt automatically even if the YData and XData of the reference
    %plot change (providing that they do not change to [] or all NaN vectors).
    
    % This method should be compileable into stand alone apps.
    
    % To modify the position, color of scrollbars see Appearance.m
   
    % MADE BY ELLEN FRANCES ZAKRESKI 2016.
    
    %% KEY PRESS 
    % || Key       || Description
    % ||(+modifier)||
    %----------------------------------------------------------------
 %PANING-------------------------------------------------------------
    % || >         || pan right  (moves both left and right sides of patch)
    % || <         || pan left
    % || ^         || pan up     (moves both top and bottom sides of patch)
    % || v         || pan down
 %ZOOMING (horizontally)---------------------------------------------
    % || > shift   || zoom in horizontally (moves only left side of patch right)
    % || < shift   || zoom out horizontally (moves only left side left)
    % || > command || zoom out horizontally  (moves right side right)
    % || < command || zoom in horizontally (moves right side left) 
 %ZOOMING (vertically)-----------------------------------------------
    % || v shift   || zoom in vertically (moves only top side of patch down)
    % || ^ shift   || zoom out vertically (moves top side of patch up)
    % || v command || zoom out vertically  (move bottom side of patch down)
    % || ^ command || zoom in vertically (moves bottom side up) 
    
    properties
        referenceAx
        axdyn %axisdynamic({'X','Y'},[2,20],refax);
    end

    methods
        function obj=superscroll(referenceAx,axesList)
            %e.g. superscroll(referenceAx,{'X'}) or superscroll(referenceAx,{'X','Y'})
            if nargin ~= 0
                %determine movespeed
                axesList=cellstr(axesList);
                movespeed = NaN(size(axesList));
                movespeed(strcmp(axesList,'X'))=2; %movespeed for X-axis
                movespeed(strcmp(axesList,'Y'))=20; %movespeed for Y-axis
                
                obj.referenceAx = referenceAx;
                obj.axdyn = axisdynamic(axesList,movespeed,obj.referenceAx);
            end
        end
        
        function autoscrollbar(obj,plotHandle)
            %automatically adds all the functions
            %make an X axis and Y axis scroll bar
            %the scrollbar will be put behind axeslist.
            buildscrollerBehind(obj);
            
            %make scrollbar speed (movespeed) adjustable so that when the user
            %clicks on the patch, they can change the percentage of data moved
            %through.
            make_moveSpeed_adjustable(obj)
            
            %Make scrollbar depend on plot handle so if
            %the reference plot handle data limits change, the scrollbar limits
            %will automatically update.
            dependOnPlot(obj,plotHandle);
        end
        
        function set.referenceAx(obj,newval)
            assert(strcmp(newval.Type,'axes'),'referenceAx must be axes handle');
            set(newval,'nextplot','add','xlimmode','manual','ylimmode','manual');
            obj.referenceAx = newval;
        end
        
        function buildscrollerBehind(obj)
            arrayfun(@putTrackerBehindReference,obj.axdyn);
            arrayfun(@makePatch,obj.axdyn);
            prepFigure(obj,obj.axdyn(1).Tracker.Parent)
        end
        
        function prepFigure(obj,Parent) %add keypress and SizeChangedFcn
            %get the handle to the figure that must be in focus for arrow
            %keys to trigger pan/zoom (scrollbar moves)
            set(Parent,'sizechangedfcn',@(src,ev)SizeChangedFcn(obj));
            while ~strcmp(Parent.Type,'figure') %in case parent is uipanel
                Parent = Parent.Parent;
            end

            set(Parent,'keypressfcn',@(src,ev)keypress(obj,ev))
        end
        function keypress(obj,KeyPressData)
            wasmatched = false; n=0;
            while ~wasmatched && n<length(obj.axdyn)
                n=n+1;
                wasmatched = keypress(obj.axdyn(n,1),KeyPressData);
            end
        end
        function SizeChangedFcn(obj)
            arrayfun(@placeScrollbar,obj.axdyn) %see Appearance.m
        end
        %add ons
        function dependOnPlot(obj,ploth)
            % make scroller limits depend on data limits of plot handle or
            % other object with X/YData
            for n = 1:length(obj.axdyn)
                obj.axdyn(n,1).referenceLim = ploth;
            end
        end
        
        function make_moveSpeed_adjustable(obj)
            %User can adjust scroll speed by clicking on patch
            make_moveSpeed_adjustable(obj.axdyn)
        end
        
        
    end

    
end

