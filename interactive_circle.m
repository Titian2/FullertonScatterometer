function interactive_circle()
    %--- Load & show image -----------------------
    I = imread('peppers.png');      % ← change to your filename
    hFig = figure('Name','Interactive Circle',...
                  'NumberTitle','off',...
                  'KeyPressFcn',@onKeyPress);
    hAx = axes('Parent',hFig);
    imshow(I,'Parent',hAx); hold(hAx,'on');
    
    %--- Initial settings ------------------------
    defaultRadius = 20;               % starting radius (pixels)
    previewColor   = [0 1 0];         % green dashed preview
    finalColor     = [1 0 0];         % red solid final
    theta = linspace(0,2*pi,100);     % circle resolution
    
    %--- Preview circle (as a line) -------------
    hPreview = plot(hAx, nan, nan, '--', ...
                    'Color', previewColor, ...
                    'LineWidth', 1, ...
                    'Visible','off');
    
    %--- UI buttons for + / – -------------------
    uicontrol(hFig,'Style','pushbutton','String','+','FontSize',14,...
        'Units','normalized','Position',[0.02 0.02 0.06 0.06],...
        'Callback',@(s,e)changeRadius(+0.1));
    uicontrol(hFig,'Style','pushbutton','String','–','FontSize',14,...
        'Units','normalized','Position',[0.10 0.02 0.06 0.06],...
        'Callback',@(s,e)changeRadius(-0.1));
    
    %--- Mouse callbacks -------------------------
    set(hFig, ...
        'WindowButtonMotionFcn',@onMouseMove, ...
        'WindowButtonDownFcn',  @onMouseClick);
    
    %--- Nested callback functions --------------
    function onMouseMove(~,~)
        C = get(hAx,'CurrentPoint');
        cx = C(1,1);  cy = C(1,2);
        xlim = get(hAx,'XLim'); ylim = get(hAx,'YLim');
        if cx>=xlim(1) && cx<=xlim(2) && cy>=ylim(1) && cy<=ylim(2)
            xCircle = cx + defaultRadius*cos(theta);
            yCircle = cy + defaultRadius*sin(theta);
            set(hPreview, ...
                'XData', xCircle, ...
                'YData', yCircle, ...
                'Visible','on');
        else
            set(hPreview,'Visible','off');
        end
        drawnow limitrate
    end

    function onMouseClick(~,~)
        C = get(hAx,'CurrentPoint');
        cx = C(1,1);  cy = C(1,2);
        xCircle = cx + defaultRadius*cos(theta);
        yCircle = cy + defaultRadius*sin(theta);
        plot(hAx, xCircle, yCircle, ...
             'Color', finalColor, ...
             'LineWidth', 1);
    end

    function onKeyPress(~,evt)
        switch evt.Key
            case {'add','equal'}      % numpad + or =
                changeRadius(+0.1);
            case {'subtract','hyphen'} % numpad – or -
                changeRadius(-0.1);
        end
    end

    function changeRadius(delta)
        defaultRadius = max(0, defaultRadius + delta);
        title(hAx, sprintf('Radius = %.1f px', defaultRadius));
    end
end