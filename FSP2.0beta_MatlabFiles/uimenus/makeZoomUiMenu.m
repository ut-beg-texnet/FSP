% toggle Zoom
% By Rall
% Stanford

% make zoom uimenu
function makeZoomUiMenu(hDV)

TopZoomUimenu(1)=uimenu(hDV.hfig,'Label','Zoom');
% Zoom on/off menu
uiHandle1(1) = uimenu('parent',TopZoomUimenu,'Label','zoom',...
    'callback',{@toggleZoomCallback,hDV},'Checked','off','tag','zoom');
uiHandle1(2) = uimenu('parent',TopZoomUimenu,'Label','Data Cursor',...
    'callback',{@toggleZoomCallback,hDV},'Checked','off','tag','datacursor');

end

% callback to change figure zoom setting
function toggleZoomCallback(src,~,hDV)

siblings = get(get(src,'parent'),'children'); % uimenus at same level
siblings=siblings(siblings~=src); % at same level, and not this menu
% get(src,'tag')
cursormode = datacursormode(hDV.hfig);
switch get(src,'tag')
    case 'zoom'
        set(siblings,'checked','off')
        set(cursormode,'Enable','off')
        switch get(src,'checked')
            case 'off' %allow zooming
                zoom(hDV.hfig,'on');
                set(src,'checked','on')
            case 'on' %don't allow zooming
                zoom(hDV.hfig,'off');
                set(src,'checked','off')
        end
    case 'datacursor'
        zoom(hDV.hfig,'off');
        set(siblings,'checked','off')
        switch get(src,'checked')
            case 'off' %allow data cursor
                set(cursormode,'Enable','on')
                set(src,'checked','on')
            case 'on' %don't allow data cursor
                set(cursormode,'Enable','off')
                set(src,'checked','off')
        end
        
        
end
end