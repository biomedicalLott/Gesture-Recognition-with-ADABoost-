
function varargout = rLottADABoost(varargin)

%% This is an Interface for capturing digital rLottADABoost during sketching
%  wriiten for MAE 527 Intelligent CAD Interface course, Spring 2022
%  Students are supposed to use this platform as a base for thire class
%  project ad HWs. Furthere improvement of this interface should be done by
%  students.
%
%  ET Esfahani
%  Department of Mechanical and Aerospace Engineering,
%  University at Buffalo, State University of New York
% % % % % % % % % % ROBERT LOTT % % % % % % % % % %
% % % % % % % % % % 4/01/2022 % % % % % % % % % %
%% WHAT IT SUPPOSED TO DO
% Does simple template matching using euclidean distance between templates
% and a drawn shape 

% Edit the above text to modify the response to help rLottADABoost

% Last Modified by GUIDE v2.5 08-Jan-2023 16:46:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @rLottADABoost_OpeningFcn, ...
    'gui_OutputFcn',  @rLottADABoost_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before rLottADABoost is made visible.
function rLottADABoost_OpeningFcn(hObject, eventdata, handles, varargin)
warning('off','all')
disp('all warnings disabled')
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rLottADABoost (see VARARGIN)

% Choose default command line output for rLottADABoost
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes rLottADABoost wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global classifiers;
global ghandles;
global points;
global pID;
global drawing;
global sketchID;
global figures;
global figCount;
% createClassifierFileFromSampleMats();
loadClassifiers();

figures(1) = gca;

figCount = 1;
axes(figures(1));
% global plots;
ghandles = handles;
tic
axis([0 1 0 1])
hold on
set(gca,'XTick',[],'Ytick',[],'Box','on')
% set(ghandles.saveButton,'Enable','off')
set(ghandles.clearButton,'Enable','off')
set(ghandles.classifyButton,'Enable','off')


set(gcf,'WindowButtonDownFcn',@mouseDown,...
    'WindowButtonMotionFcn',@mouseMoving,...
    'WindowButtonUpFcn',@mouseUp)

points = zeros(1000,3);
pID = 0;
drawing = 0;
sketchID = 0;
% loadTemplates();
clc;

% helpdlg({'When loading in test or sample sets',...
%     'select the box, press load, and then classify',...
%     'to classify your own drawing',...
%     'please make sure all check marks are clicked off'},...
%         'Help')

% --- Executes on button press in classifyButton.
function classifyButton_Callback(hObject, eventdata, handles)
% hObject    handle to classifyButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ghandles;
global points;
global classifiers;
global pID;
drawing = adaBoostGesture(0, points(1:pID,:),0); 
set(ghandles.displayText,'string','Classifying...');

[H,L] = adaBoostGesture.getHypothesis(classifiers,drawing);
str1 = strcat("The classifier thinks it's ", string(L));
set(ghandles.displayText,'string', str1);

% 
% fprintf('%0.2f score for label %0.0f\n', [H;L]);


% classifiers = adaBoostGesture.classMeans(classifiers, 10000, 10);
% for i = 1:length(classifiers)
%    saveStruct.classifier(i) = classifiers(i).saveClassifier;
% end
% save('adaBoostClassifier.mat','saveStruct');



% 
%     for i = 1:numel(classifiers(:))
%        fprintf('Group %0.0f\n', i);
%         for j = 1:numel(classifiers(i).samples(:))
%            H = classifiers(i).samples(j).LikelyHypothesis;
%            L = classifiers(i).samples(j).LikelyLabels;
% %            fprintf('%0.2f score for label %0.0f\n', [H([1,end]);L([1,end])]);
%            
%        end
%        
%     end
% [bestTemplate, Scores] = adaBoostGesture.polarRecognize(drawing,classifiers);
% adaBoostGesture.compareToTemplates(drawing,classifiers)
% outputClass(bestTemplate,real(Scores))

% redraw(drawing);
% if((ghandles.TrainData.Value))
%     BeginTraining();
% elseif ghandles.samplePack.Value
%     classes = loadClasses();
%     beginClassifyingSampleGroup(classes);
% else
%     classes = loadClasses();
%     beginClassifyingUserInput(classes);    
% end

% --- Outputs from this function are returned to the command line.
function varargout = rLottADABoost_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in clearButton.
function clearButton_Callback(hObject, eventdata, handles)
% hObject    handle to clearButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ghandles;
global points;
global pID;
global sketchID;
global figures
global figCount;
global templateObjects;
% plots(1)
axes(figures(1));
cla
for i = 2: figCount
    clf(figures(i))
end
axes(figures(1));


points = zeros(1000,3);
pID = 0;
set(ghandles.clearButton,'Enable','off')
% set(ghandles.saveButton,'Enable','on')
set(ghandles.classifyButton,'Enable','off')
set(ghandles.displayText,'string', "Data Cleared");


function mouseDown(hObject, eventdata, handles)
global ghandles;
global points;
global pID;
global drawing;
% axes(ghandles.drawGraph);

cp = get(gca,'CurrentPoint');

if cp(1,1)<1 && cp(1,1)>0 && cp(1,2)<1 && cp(1,2)>0
    drawing = 1;
end


function mouseMoving(hObject, eventdata, handles)
global ghandles;
global points;
global pID;
global drawing;

if drawing
    t = toc;
    cp = get(gca,'CurrentPoint');
    pID = pID + 1;
    points(pID,1) = cp(1,1);
    points(pID,2) = cp(1,2);
    points(pID,3) = t;
    plot(cp(1,1),cp(1,2),'o')
    
end


function mouseUp(hObject, eventdata, handles)


global ghandles;
global points;
global pID;
global drawing;
global sketchID;

if drawing
    drawing = 0;
    sketchID = sketchID+1;
end
set(ghandles.clearButton,'Enable','on')
% set(ghandles.saveButton,'Enable','on')
set(ghandles.classifyButton,'Enable','on')



function redraw(drawing)
global figures    

    axes(figures(1));
    cla
    for i = 2: figCount
        clf(figures(i))
    end
    axes(figures(1));


    plot(drawing.X,drawing.Y,'rs')
    
function outputClass(bestTemplate,Score)
global ghandles;
    guessLen = length(bestTemplate);
    outputText = ['You drew ', num2str(bestTemplate(1)),...
        ' score: ',num2str(Score(1),3)];
    if(guessLen > 1)
        outputText = [outputText,...
            ' but it might be ', num2str(bestTemplate(2)),...
            ' score: ',num2str(Score(2),3)];
    end
        set(ghandles.displayText,'string', outputText)



% --- Executes on button press in saveButton.
% function saveButton_Callback(hObject, eventdata, handles)
% % hObject    handle to saveButton (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% global points;
% global pID;
% global sketchID;
% global ghandles;
% global classifiers;
% % fileNum = 1;
% if pID < 1000
%     points(pID+1:end,:)=[];
% end
% groundNum = 0
% fileName = ['Template_', num2str(groundNum),'.mat']
% while(isfile(fileName))
%     groundNum = groundNum +1;
%     fileName = ['Template_', num2str(groundNum),'.mat']
% end
% template = adaBoostGesture(0,points,string(groundNum), randi(100));
% saveStruct = template.saveobj()
% save(fileName,'saveStruct')

% groundNum = (get(ghandles.groundTruth,'string'));
% % fileNameTry = fileName ;
% fileNum = 0;
% fileName = ['data_', groundNum, num2str(fileNum), '.mat'];
% 
% while(isfile(fileName))
%     fileNum = fileNum +1;
%     fileName = ['data_', groundNum, num2str(fileNum), '.mat'];
% end
% save(fileName,'points')
%
% clearvars -except classifiers;
% txt = jsonencode(classifiers);
% [file] = uiputfile('digitMatrix.mat');
% save(file,'txt');
% uisave();
% disp("file Saved");

function saveStruct = saveObj(template)
saveStruct.groundTruth = template.groundTruth
saveStruct.ID = template.ID
saveStruct.X = template.X
saveStruct.Y = template.Y
saveStruct.center = template.center
saveStruct.theta = template.theta


% --- Executes on button press in loadButton.
function loadButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% LOAD FILE BY NAME: data_#
global points;
global pID;
global sketchID;
global figures;
global ghandles;
global classifiers;
%     loadTemplates()
%     createTemplateFiles()
% You can either open the training set, the sample set, or something else
% although it better be points because i don't have any logic setup for
% that yet!
%     loadFromDataPack('sampleSet.mat');
% 
% if(ghandles.TrainData.Value == 1)
%     loadFromDataPack('trainset.mat');
%     
% elseif(ghandles.samplePack.Value == 1)
%     loadFromDataPack('sampleSet.mat');
% else
%     uiopen('*.mat')
% end
% 
% disp("File ready");

set(ghandles.clearButton,'Enable','on')
% set(ghandles.saveButton,'Enable','on')
set(ghandles.classifyButton,'Enable','on')


% % % % % % % % USER MADE FUNCTIONS % % % % % % % %

function createClassifierFileFromSampleMats()
digit1 = [0:9];
digit2 = [1:10];
% classObjects = adaBoostGesture.empty(10,0);

% gets current directory
myDir = pwd;
% get Segmented Hand Data folder
source = strcat(myDir ,'\digits\');
myFiles = dir(fullfile(source,'*.mat'));
counter = 1;
for i = 1:length(digit1)
sampleObjects = adaBoostGesture.empty(10,0);
    for j = 1:length(digit2)
%         fileName = ['digit_',num2str(digit1(i)),'_',num2str(digit2(j))];
        load(myFiles(counter).name,'points');
        counter = counter+1;
        sampleObjects(j) = adaBoostGesture(0,points,digit1(i), digit2(j));
    end
    classObjects(i) = adaBoostGesture(1,[],digit1(i));
    classObjects(i).makeClassifier(sampleObjects);
end
for i = 1:length(classObjects)
   saveStruct.classifier(i) = classObjects(i).saveClassifier;
end
save('adaBoostClassifier.mat','saveStruct');


function loadClassifiers()
global classifiers;
classifiers = adaBoostGesture.empty(10,0);
    load('adaBoostClassifier.mat', 'saveStruct');
for i = 1:length(saveStruct.classifier(:))
%    newObj = adaBoostGesture(true);
%    newObj.loadobj(saveStruct(i));
% adaBoostGesture.loadClassifiers(saveStruct.classifier(i))
    classifiers(i) = adaBoostGesture.loadClassifiers(saveStruct.classifier(i));
end

    

function loadFromDataPack(fileName)
% Loads in json encoded data packs that i have already stored
global points;
global classifiers;
clear points;
%     

% uiopen(fileName)
% bigMat = jsondecode(txt);
% matSize = size(bigMat);
%          dataNums = [0:9];
%          fileNums = dataNums +1;
% %         txt =[];
% %         len = length(dataNums);
% templateHolder = adaBoostGesture.empty(10,0);
% for i = dataNums%1: matSize(1)
%     disp(num2str(i))
%     for j = 1%1: matSize(2)
%         points = [bigMat(i,j).x,bigMat(i,j).y,bigMat(i,j).t];
%     fileName = ['template_',num2str(i)];    

% fileName = ['digit_',num2str(i),'_',num2str(1),'.mat'];
%         load(fileName,'points');
%         fileStruct = load(fileName);
%         pointsVar  = fileStruct.points;
%         templateHolder(i+1,j) = adaBoostGesture(points,i);
%         template = adaBoostGesture(0,pointsVar,i);
%       templateHolder(i+1) = template;
%         clear pointsVar;
%     end


%     save(fileName, 'template')
% end

% classifiers = templateHolder;
function dataIDs = loadTemplateDict()
%             loads from csv, i did this because i kept forgetting what the
%             file was
dataIDs = readcell('templateDict.csv')
% dataIDs = 

function loadTemplates()
%             general load templates for app start
    global classifiers;
    dataIDs  = loadTemplateDict()
%     app.ClassesListBox.Items = string(dataIDs');
    templateCount = length(dataIDs);
%     dataNums = [1:templateCount];
    obj = [];
%             it will iterate through all template files based on the
%             template "dictionary" i already created. That way we don't
%             need to remember what files go where
    for i = 1:templateCount
        fileName = 'Template_' + string(dataIDs(i)) + '.mat'
        load(fileName, 'saveStruct');
%         points = [saveStruct.oldX,saveStruct.oldY,saveStruct.oldX];
%         template = adaBoostGesture(0,points,saveStruct.groundTruth,saveStruct.ID);
%         fileName = strjoin(['Template_',template.groundTruth,'.mat']);
%         fileName = strrep(fileName,' ','')
%         saveStruct = template.saveobj()
%         save(fileName, 'saveStruct')
        template = adaBoostGesture.loadobj(saveStruct);
        obj = [obj ;template];

    end
%         close all
    
    classifiers = obj;

function createTemplateFiles()
         dataNums = [0:9];
templateHolder = z.empty(10,0);
for i = dataNums%1: matSize(1)
    disp(num2str(i))
fileName = ['digit_',num2str(i),'_',num2str(1),'.mat'];
        fileStruct = load(fileName);
        pointsVar  = fileStruct.points;
        template = adaBoostGesture(0,pointsVar,i);
      templateHolder(i+1) = template;
      saveStruct = saveobj(template);
fileName = ['templates\template_',num2str(i),'.mat'];
  save(fileName, 'saveStruct')
end

function saveStruct = saveobj(template)
%             savedObj.
saveStruct.groundTruth = template.groundTruth     
saveStruct.ID = template.ID             
saveStruct.X = template.X             
saveStruct.Y = template.Y             
saveStruct.center = template.center          
saveStruct.theta = template.theta         

% 
% function combineTemplateFiles()
%    dataNums = [0:9];
%    
%     file = [];
%     for i = dataNums
%     fileName = ['template_',num2str(i),'.mat'];
%     file = [file;load(fileName.name)]
%         
%     end
%     save('Templates.mat',file);



% % % % % % % % UNUSED FUNCTIONS FROM MATLAB APP % % % % % % % %




% --- Executes during object creation, after setting all properties.
function groundTruth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to groundTruth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function speedThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speedThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function angleThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angleThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function curveThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to curveThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function polyText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to polyText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showOptions.
function showOptions_Callback(hObject, eventdata, handles)
% hObject    handle to showOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showOptions
global ghandles 
if(ghandles.showOptions.Value)
set(ghandles.loadButton,'Enable','on')
set(ghandles.loadButton,'Visible','on')
set(ghandles.saveButton,'Enable','on')
set(ghandles.saveButton,'Visible','on')
set(ghandles.showOptions,'Visible','on')
set(ghandles.groundTruth,'Visible','on')
return
end
set(ghandles.loadButton,'Enable','off')
set(ghandles.loadButton,'Visible','off')
set(ghandles.saveButton,'Enable','off')
set(ghandles.saveButton,'Visible','off')
set(ghandles.showOptions,'Visible','off')
set(ghandles.groundTruth,'Visible','off')
