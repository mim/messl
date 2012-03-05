function replace_eps_figure_strings_with_tex(eps_fig_file, strings, tex_strings, preamble_files)
% replace_eps_figure_strings_with_tex(eps_fig_file, strings,
%                                     tex_strings, preamble_files)
%
% Use psfrag to replace strings in an eps file with TeX rendered
% strings (i.e. TeX rendered math).  eps_fig_file is the *full* path
% to an eps figure, strings is a cell array of strings in the figure
% to be replaces, tex_strings are the strings to replace them with.
% preamble_files is an optional cell array list of TeX headers
% containing any additional TeX preamble needed to properly render the
% given tex_strings (i.e. using any additional packages or defining
% macros Note that the amsmath package is included by default).

if nargin < 4;  preamble_files = {};  end
if ~iscell(preamble_files)
  preamble_files = {preamble_files};
end

if length(strings) ~= length(tex_strings)
  error('strings and tex_strings must be he same length');
end

tmpdir = tempname();
mkdir(tmpdir);

preamble_str = '';
for n = 1:length(preamble_files)
  preamble_str = sprintf('%s\\\\input{%s}\n', preamble_str, preamble_files{n});
end
psfrag_str = '';
for n = 1:length(strings)
  % usage: \psfrag{text}[posn][psposn][scale][rotate]{formula}
  psfrag_str = sprintf('%s\\\\psfrag{%s}{%s}\n', psfrag_str, ...
      strings{n}, strrep(tex_strings{n}, '\', '\\'));
end

[figpath figname ext] = fileparts(eps_fig_file);

tex = ['\\documentclass[10pt]{article}\n', ...
       '\\pagestyle{empty}\n', ...  %% don''t show the page number\n', ...
       '\\usepackage{psfrag}\n', ...
       '\\usepackage{graphicx}\n', ...
       preamble_str, ...
       '\\usepackage{amsmath}\n', ...
       '\\usepackage{bm}\n', ...
       '\\begin{document}\n', ...
       '\\begin{figure}\n', ...  
       psfrag_str, ...
       '\\includegraphics[scale=1.0]{' figname '}\n', ... 
       '\\end{figure}\n', ...
       '\\end{document}\n'];
fid = fopen(fullfile(tmpdir, 'figure.tex'), 'w');
fprintf(fid, tex);
fclose(fid); 
copyfile(eps_fig_file, fullfile(tmpdir, [figname ext]));

system(sprintf(['cd %s;' ...
                'latex figure > /dev/null;' ...
                'dvips -Ppdf figure;' ...
                'ps2epsi figure.ps figure.eps;'...
                'cp figure.eps %s;'...
                'rm -rf %s;'], tmpdir, eps_fig_file, tmpdir));
