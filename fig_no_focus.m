function fig_no_focus(f)

% Sets the current figure without raising its window.

try
  set(0, 'CurrentFigure', f);
catch
  figure(f);
end
