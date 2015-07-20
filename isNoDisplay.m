function v = isNoDisplay()

% Return true if matlab was called with the -nodisplay option, false
% otherwise

[~,d] = system('echo $DISPLAY');
if ~reMatch(d, ':\d+')
    v = true;
    return
end

% From http://www.mathworks.com/matlabcentral/newsreader/view_thread/136261
ss = get(0, 'Screensize');
v = all(ss(3:4) == [1 1]);
