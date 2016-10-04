function ACdata = LoadAcFile(WF_path,filenumber,numCHR,numSFpfile)

% Load acoustic file and reshape according to the
% number of receivers
% WF_path: location of the files
% filenumber: which file is being loaded
% numCHR: number of receivers 
% numSFpfile: number of 'superframes' per file (Verasonics jargon)

ACfilename = [WF_path num2str(filenumber) '.ac'];
fid = fopen(ACfilename,'r');
ACdata = fread(fid,'int16');
fclose(fid);

% reshape to get one column per channel
ACdata = reshape(ACdata,[],numCHR,numSFpfile);   % 3D matrix with WF vs numCHR vs number of SF
ACdata = permute(ACdata,[1 3 2]);               % put numCHR as the last dimension before reshaping
ACdata = reshape(ACdata,[],numCHR,1);            % WF vs numCHRs

end
