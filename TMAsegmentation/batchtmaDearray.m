function batchtmaDearray(parentFolder,varargin)

ip = inputParser;
ip.addParamValue('buffer',1.2,@(x)(numel(x) == 1 & all(x > 0 )));  
ip.addParamValue('downSampleFactor',4,@(x)(numel(x) == 1 & all(x > 0 )));  
ip.addParamValue('writeTiff','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('writeMasks','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('outputFiles','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('outputCenters','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('useGrid','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('sample','TMA',@(x)(ismember(x,{'TMA','tissue'})));
ip.addParamValue('Docker',false,@islogical);
ip.addParamValue('DockerParams',0,@isstruct);
ip.addParamValue('outputChan',0,@(x)(isnumeric(x))); 
ip.addParamValue('cluster',false,@islogical);
ip.parse(varargin{:});          
p = ip.Results;  

% override default values because using Docker config file
if p.Docker==1 && isstruct(p.DockerParams)
    parentFolder = p.DockerParams.parentPath;
    modelPath = [p.DockerParams.modelPath filesep];
    outputPath = p.DockerParams.analysisPath;
else
    if nargin<1
        parentFolder =uigetdir;
    end
    modelPath ='';
    outputPath = parentFolder;
        
end

if isempty(parentFolder)
    error('You must select a folder to continue')
end

finalFolderList = [];
folders = dir(parentFolder);
for iFolder = 1:numel(folders)
    if isdir([parentFolder filesep folders(iFolder).name]) && ~isequal(folders(iFolder).name,'.') && ~isequal(folders(iFolder).name,'..')
        finalFolderList{end+1} = folders(iFolder).name;
    end
end
disp(['Found ' num2str(numel(finalFolderList)) ' folders(s)' ])

for iFolder = 1:numel(finalFolderList)
        pathName = finalFolderList(iFolder);
        finalFileList = [];
        fileListing = dir([parentFolder filesep char(pathName) filesep '*.ome.tif']);
        for iFile = 1:numel(fileListing)
            if ~isequal(fileListing(iFile).name,'.') && ~isequal(fileListing(iFile).name,'..')
                finalFileList{end+1} = fileListing(iFile).name;
            end
        end
        
        
        disp(['Found ' num2str(numel(finalFileList)) ' file(s) in folder ' num2str(iFolder) ' of ' int2str(numel(finalFolderList))])

        if ~isempty(finalFileList)
            for iFile = 1:numel(finalFileList)
                disp(['Processing filename ' finalFileList{iFile} ' - file ' num2str(iFile) ' of ' num2str(numel(finalFileList))])
               
                    tmaDearray([parentFolder filesep char(pathName) filesep finalFileList{iFile}],'buffer',p.buffer,...
                        'writeTiff',p.writeTiff,'writeMasks',p.writeMasks,'outputFiles',p.outputFiles,'modelPath', modelPath,...
                        'outputPath',[outputPath filesep 'dearray'],'outputChan',p.outputChan,'useGrid',p.useGrid,...
                        'cluster',p.cluster,'Docker',p.Docker,'sample',p.sample,'downSampleFactor',p.downSampleFactor);
                               
            end
        end
end

    disp('All Cores Complete!')