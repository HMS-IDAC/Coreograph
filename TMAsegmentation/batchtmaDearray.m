function batchtmaDearray(parentFolder,varargin)

ip = inputParser;
ip.addParamValue('buffer',1.5,@(x)(numel(x) == 1 & all(x > 0 )));  
ip.addParamValue('writeTiff',true,@islogical);
ip.addParamValue('writeMasks',true,@islogical);
ip.addParamValue('outputFiles',true,@islogical);
ip.addParamValue('Docker',false,@islogical);
ip.addParamValue('DockerParams',0,@isstruct);
ip.addParamValue('outputChan',1,@(x)(all(x > 0)));   
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
        fileListing = dir([parentFolder filesep char(pathName) filesep 'registration' filesep '*.ome.tif']);
        disp(['Found ' num2str(numel(fileListing)) ' file(s) in folder ' num2str(iFolder) ' of ' int2str(numel(finalFolderList))])

        if ~isempty(fileListing)
            for iFile = 1:numel(fileListing)
                disp(['Processing file ' num2str(iFile) ' of ' num2str(numel(fileListing))])
                tmaDearray([parentFolder filesep char(pathName) filesep 'registration' filesep fileListing(iFile).name],'buffer',p.buffer,'writeTiff',p.writeTiff,...
                    'writeMasks',p.writeMasks,'outputFiles',p.outputFiles,'modelPath', modelPath,'outputPath',outputPath,'outputChan',p.outputChan,'Docker',p.Docker);
            end
        end
end

    disp('All Cores Complete!')