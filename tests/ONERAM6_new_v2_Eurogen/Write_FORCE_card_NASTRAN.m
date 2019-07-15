%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   NASTRAN FORCE CARD WRITING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc;

FileToRead='Output_convergence_iter5/Loads_Iter0.dat';

A = readmatrix(FileToRead);

OutputFile= 'FORCE_file.bdf';

fid=fopen(OutputFile,'w');

% Force file layout
% $       $       $       $       $       $       $       $       $
% FORCE          1      20       0      1.   0.         0.6.9807+3
%              set    node refSYS sc. fact   Lx       Ly      Lz

for i=1:length(A)
    
    fprintf(fid,'$       $       $       $       $       $       $       $       $\n');
    fprintf(fid,'%-8s%8i%8i%8i%8.5f%8.5e%8.5e%8.5e\n','FORCE', 1,A(i,1)+ 1, 0, 1,A(i,3),A(i,5),A(i,7));
    fprintf(fid,'\n');
    
end

fclose(fid);


