function [verts,faces,volume] = msmsSurface0(xyzr) 
% This script computes the MSMS surface using the msms software. 
% xyzr: the coordinates of atoms and their radii. 

%eval(['save -ascii ' code '.xyzr xyzr']);
eval(['save -ascii temp.xyzr xyzr']);
%[status, cmdout] = system(['msms -if ' code '.xyzr -of temp']); 
[status, cmdout] = system(['msms -if temp.xyzr -of temp']);
%lines = splitlines(string(cmdout));
%tokens = split(strip(lines(end-2)));
%if status == 0 % success
where = find(cmdout==':', 1, 'last');
nextNewln = find(cmdout(where+1:end)==10);
%volume = str2double(tokens(3)); % total SES volume
volume = str2num(cmdout(where+1:where+nextNewln-1));
[verts] = readMSMSoutput('temp.vert');
[faces] = readMSMSoutput('temp.face');

if isempty(volume) || volume<0 volume = 0; end % MSMS failed

function [xyz] = readMSMSoutput(fname)
fid=fopen(fname,'r');
A = fread(fid);
fclose(fid);
newln = find(A==10); % newline characters
A(newln(4:end))=' ';
A = A(newln(3)+1:end)';
data = str2num(char(A));
xyz = data(:, 1:3);
%lines = splitlines(string(char(A')));
%lines = lines(4:end-1);
%tokens = split(strip(lines)); % split into columns
%N = size(lines,1);
%xyz = str2double(tokens(:, 1:3));