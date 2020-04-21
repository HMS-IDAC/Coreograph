function batchtmaDearrayTest(parentFolder,varargin)

ip = inputParser;
ip.addParamValue('buffer',2,@(x)(numel(x) == 1 & all(x > 0 )));  
ip.addParamValue('downSampleFactor',5,@(x)(numel(x) == 1 & all(x > 0 )));  
ip.addParamValue('writeTiff','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('writeMasks','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('outputFiles','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('outputCenters','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('useGrid','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('sample','TMA',@(x)(ismember(x,{'TMA','tissue'})));
ip.addParamValue('Docker',false,@islogical);
ip.addParamValue('DockerParams',0,@isstruct);
ip.addParamValue('outputChan',1,@(x)(isnumeric(x))); 
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

fileListing = dir([parentFolder filesep '**' filesep '*.ome.tif']);   

for iFile = 1:numel(fileListing)
     
        disp(['Processing filename ' fileListing(iFile).name])
tic
            tmaDearrayTest([fileListing(iFile).folder filesep fileListing(iFile).name ],'buffer',p.buffer,...
                'writeTiff',p.writeTiff,'writeMasks',p.writeMasks,'outputFiles',p.outputFiles,'modelPath', modelPath,...
                'outputPath',[outputPath filesep 'dearray'],'outputChan',p.outputChan,'useGrid',p.useGrid,...
                'cluster',p.cluster,'Docker',p.Docker,'sample',p.sample,'downSampleFactor',p.downSampleFactor);
      toc
end

    disp('All Cores Complete!')