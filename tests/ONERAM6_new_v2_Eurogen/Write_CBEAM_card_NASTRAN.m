%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   NASTRAN CBEAM CARD WRITING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc;

OutputFile= 'CBEAM_file.bdf';

fid=fopen(OutputFile,'w');

% CBEAM file layout
% $       $       $       $       $       $       $       $       $
%CBEAM          1       1       1       2 1.00000 0.00000 0.00000
id = 20;
for i=1:20
    
    fprintf(fid,'$       $       $       $       $       $       $       $       $\n');
    fprintf(fid,'%-8s%8i%8i%8i%8i%8f%8f%8f\n','CBEAM',id,4,i,i+20,0,1,0);
    id = id+1;
    fprintf(fid,'%-8s%8i%8i%8i%8i%8f%8f%8f\n','CBEAM',id,4,i,i+40,0,1,0);
    id = id+1;
    fprintf(fid,'%-8s%8i%8i%8i%8i%8f%8f%8f\n','CBEAM',id,4,i,i+60,0,1,0);
    id = id+1;
    fprintf(fid,'%-8s%8i%8i%8i%8i%8f%8f%8f\n','CBEAM',id,4,i,i+80,0,1,0);
    fprintf(fid,'\n');
    id = id+1;
end

fclose(fid);