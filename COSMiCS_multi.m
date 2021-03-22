    clear all;

clc;


disp('****************************************************************************************************');
disp('**  STRUCTURAL ANALYSIS OF MULTICOMPONENT AMYLOID SYSTEMS BY CHEMOMETRIC SAXS DATA DECOMPOSITION  **');
disp('****************************************************************************************************');
disp(' '); 

      

%% Selection of folders of multiple experiments

% Input the number of experiments

disp(   'Please enter the number of experiments performed'); 
numexp = input('Number of experiments[1] ');

if isempty(numexp)
    numexpv = 1;
end
numexpv = numexp;

% Selection of multiple folders and reading the files of all the
% experiments

for k=1:numexpv
    fprintf(  'Select the folder containing the data from experiment number %d (press any key to continue)\n', k);
    pause;

    folderInput = uigetdir('/home/');
    foldernames{k} = folderInput;
    disp([  'Directory: ', folderInput]);
    disp(' ');
    disp(' ');
    display(folderInput);
    disp('Files inside');
    disp(' ');
    listing = dir(folderInput);
    for i=3:5
        display(listing(i).name);
    end
    display('...');
    for i=(length(listing)-3):length(listing)
        display(listing(i).name);
    end
    
    
    disp(' ');
    str =   'What is the pattern of the files (e.g. curve*.dat ) ?? ';
    patron = input(str,'s');
    if isempty(patron)
        patron = 'curve*.dat';
    end

end

% Selection of output folder


disp(  'Where do you want to save the output? (press any key to continue)');
pause;

folderOutput = uigetdir('../');
folderOutput = [folderOutput, '/'];

OutputFolder = folderOutput;


disp([  'Directory: ', folderOutput]);
disp(' ');

%% Selection of the pattern of the files



%% Charge curves


% Initialization
OutFilesList = [];
foldernameslist = {};
NumberofCurvesPerExperimentlist = zeros(numexpv+1,1);
for a = 1:numexpv
    filePattern = fullfile(foldernames{a}, patron);
    OutFiles = dir(filePattern);
    NumberofCurvesPerExperiment = length(OutFiles);
    [filepath,name,ext] = fileparts(filePattern);
    foldernameslist = [foldernameslist;filepath];
    OutFilesList = [OutFilesList; OutFiles];
    NumberofCurvesPerExperimentlist(a+1,1) = NumberofCurvesPerExperiment;

end


ExperimentRanges(:,1) = cumsum(NumberofCurvesPerExperimentlist);


curvesInitial = cell(length(OutFilesList),1);
curves = cell(length(OutFilesList),1);
    

filePattern = fullfile(foldernames{a}, patron);
OutFiles = dir(filePattern);


basefilename = OutFilesList(1).name;
foldername = foldernames{1};

fullFileName = fullfile(foldername,basefilename);
textillo = ['dbtype ' fullFileName, ' 1:10'];
eval(textillo);

% Select number of header lines
% disp(' ');
% num = input('Number of header lines [3]: ');
num = 3;
if isempty(num)
    num = 3;
end
headerLines = num;

% Add SAXS curves

% disp(' ');
% num = input('Number of columns [3]: ');
num = 3;
if isempty(num)
    num = 3;
end
columns = num;
columnsvariable = ['indata = textscan(fid, '''];
for k = 1:columns
    columnsvariable = [columnsvariable, '%f '];
end
columnsvariable = [columnsvariable, ''', ''HeaderLines'',headerLines);'];

for k = 1:numexp
    for curvenumber = ExperimentRanges(k)+1:ExperimentRanges(k+1)
        basefilename = OutFilesList(curvenumber).name;
        foldername = foldernameslist{k};
        fullFileName = fullfile(foldername,basefilename);
            fid = fopen(fullFileName,'rt');
            eval(columnsvariable);
            fclose(fid);
            for i = 1:columns
                curvesInitial{curvenumber,1}(:,i) = indata{i};
            end
    end
    
end


% % Check length of the curves
% 
% % Sometimes the curves have different number of datapoints, we have to
% % put all with the same number of datapoints. We compare the length of all
% % curves with the first one. If someone is not equal long we remove the
% % final points so we get all curves with the same length. Check S values
% % after


for i=1:length(curvesInitial)
    long(i,:) = size(curvesInitial{i,1});
    esIgual(i,1) = isequal(long(1,1),long(i,1));
end

if all(esIgual(:)) == 0  %alguna de las curvas tiene distinto tamanyoo/Some of the curves have different s value
    puntosCurvasExp = min(long(:,1));
    for i=1:length(curvesInitial)
        compS(i,1) = isequal(curvesInitial{1,1}(1:puntosCurvasExp,1),curvesInitial{i,1}(1:puntosCurvasExp,1));
    end
    if all(compS(:,1)) == 1 %los valores S concuerdan/the S values agree
        for i=1:length(curvesInitial)
            curvesFinal{i,1}(1:puntosCurvasExp,:) = curvesInitial{i,1}(1:puntosCurvasExp,:);
        end
    else % los valores S no concuerdan/the S values don't agree
        clear idx;
        for n=1:length(curvesInitial)
            firstValoresSall(1,n) = curvesInitial{n,1}(1,1);
            [biggerS index] = max(firstValoresSall);
            idx(n,1) = find(curvesInitial{n}(:,1) == biggerS );
            curvesFinal{n,1} = curvesInitial{n,1}(idx(n,1):end,:);
        end
        for i=1:length(curves)
            long(i,:) = size(curves{i,1});
            esIgual(i,1) = isequal(long(1,1),long(i,1));
        end
        if all(esIgual(:)) == 0  %alguna de las curvas tiene distinto tamanyo
            puntosCurvasExp = min(long(:,1));
            for j=1:length(curves)
                curvesFinal{j,1} = curves{j,1}(1:puntosCurvasExp,:);
            end
        else
            for j=1:length(curves)
                curvesFinal{j,1} = curves{j,1};
            end
        end
    end
end



A = exist('curvesFinal','var');
if A==0
    curvesFinal = curvesInitial;
end
for i=1:length(curvesFinal)
    Intensities(:,i) = curvesFinal{i,1}(:,2);
    if columns >=3
        Errors(:,i) = curvesFinal{i,1}(:,3);
        
    end
end

%% Normalize Intensities
%Get the list of I0 values
I0values = Intensities(1,:);
I0values = transpose(I0values);
I0minperExperiment = min(I0values);
for i = 1:length(Intensities(1,:))
Intensities(:,i) = Intensities(:,i)*100/I0minperExperiment;
Errors(:,i) = Errors(:,i)*100/I0minperExperiment;
end
valoresS(:,1) = curvesFinal{1,1}(:,1);
%%
% Smoothen the Intensity Data


curvasTotales = length(curvesFinal);

fprintf('\nYou added %s files\n\n', num2str(curvasTotales));
% % Scale the two experiments
scalefactor = zeros(ExperimentRanges(2),1);
% get the lowest curve
[M,I] = min(valoresI(1,:));
for i = 1:length(valoresI(1,:))
    f1 = sum(valoresI(:,I).*valoresI(:,i)./valoresE(:,I).^2);
    f2 = sum(valoresI(:,i).^2./valoresE(:,I).^2);
%scalefactor(i,1) = valoresI(1,i)/min(valoresI(1,:));
    scalefactor(i,1) = f2/f1;
% scalefactor = ones(ExperimentRanges(2),1);
valoresI(:,i) = valoresI(:,i)/scalefactor(i,1);
valoresE(:,i) = valoresE(:,i)/scalefactor(i,1);
end
scalefilename = fullfile(folderOutput, 'scale.txt');
dlmwrite(scalefilename, scalefactor);

% Plot curves in semi-logarithmic scale

 disp(' ');
 str = 'Do you want to see the curves (semi-logarithmic scale)?? [Y]/N: ';
 respuesta = input(str,'s');
 if isempty(respuesta)
     respuesta = 'y';
 end
 if respuesta =='y'
     seleccion = 2;
     plots;
 end

seleccion = 2;
plots;
disp('Curves (semi-logarithmic scale)... Press any key to continue ');
%pause;


%% Check units



control = 0;
while control == 0
%     disp(' ');
%     str = 'Units of your curves: 1/Angstrom [A] or 1/nanometers (N): ';
%     unidades = input(str,'s');
    units = 'A';
    if isempty(units)
        units = 'A';
    end
    if units == 'a' || units == 'A'
        units = 'A';
        control = 1;
    elseif units == 'n' || units == 'N'
        units = 'N';
        control = 1;
    end
end
control = 0;
while control == 0
%     disp(' ');
%     str = 'Do you want to change the units?? Y/[N]: ';
%     respuesta = input(str,'s');
    respuesta = 'N';
    if isempty(respuesta)
        respuesta = 'n';
    end
    if respuesta =='y' || respuesta =='Y'
        control = 1;
        if units == 'A'
            units = 'N';
            valoresS = valoresS*10;
            disp('Now the units are nanometers');
        elseif units == 'N'
            units = 'A';
            valoresS = valoresS/10;
            disp('Now the units are Angstroms');
        end
    elseif respuesta == 'n' || respuesta == 'N'
            control = 1;
    end
end





%% Remove points at the beginning
    close all;
    seleccion = 2;
    plots;
% disp(' ');
% str = 'Do you want eliminate some data points at the beginning?? [Y]/N: ';
% respuesta = input(str,'s');
respuesta = 'N';
if isempty(respuesta)
    respuesta = 'y';
end
if respuesta=='y' || respuesta =='Y'
    valoresIcopy = Intensities;
    valoresScopy = valoresS;
    if columns >=3
    valoresEcopy = Errors;
    end
    likeIt = 'n';
    
    while likeIt~= 'y'
        seleccion = 2;
    plots;
        disp(' ');
        elimPoints = input('Number of data points [0]: ');
        close all;
        if isempty(elimPoints)
            elimPoints = 0;
        end
        if elimPoints ~=0
            Intensities(1:elimPoints,:) = [];
            valoresS(1:elimPoints,:) = [];
            if columns >=3
                Errors(1:elimPoints,:) = [];
            end
            puntosCurvasExp = length(Intensities);
        else
            elimPoints = 0;
        end
        seleccion = 2;
        plots;
        disp(' ');
        str = 'Do you like it?? Y/[N]: ';
        likeIt = input(str,'s');
        if isempty(likeIt)
            likeIt = 'n';
        end
        if likeIt == 'Y' 
            likeIt = 'y';
        else
            Intensities = valoresIcopy;
            valoresS = valoresScopy;
            if columns >=3
                Errors = valoresEcopy;
            end
        end
    end
else
    elimPoints = 0;
end


valoresIActual = Intensities;


% NOTE: From here I will use always the values in valoresS, valoresI and
% valoresE but I kept the variable curvasFibrilacion with the original data



%% Definition of number of species (PCA)

% Principal Component Analysis

% PCA of the absolute values is the best to decide the number of species

if length(OutFilesList)<10
    pcaLong = length(OutFiles);
else
    pcaLong = 10;
end

for i=1:pcaLong
    [uAbs,sAbs,vAbs,xAbs, SigmaAbs] = pcarep(Intensities,i);
    sigmaAbs(i,1) = SigmaAbs;
end

[uAbs,sAbs,vAbs,xAbs, SigmaAbs] = pcarep(Intensities,length(Intensities(1,:)));

    
for i=1:length(sAbs)
    eigenvalues(i,1) = sAbs(i,i);
%     eigen(i,1) = i;
%     eigen(i,2) = sAbs(i,i);
end
    
for i=1:length(uAbs)
    eigenvectors(i,1) = i;
end
   
% Save eigenvectors

for i=1:pcaLong
%     j = i+1;
    eigenvectors(:,i) = uAbs(:,i);
end

% Save text files

mkdir(folderOutput);
ruta = [folderOutput,'eigenvalues.txt'];
save(ruta, 'eigenvalues', '-ASCII');

% ruta = [folderOutput,'sigma.txt'];
% save(ruta, 'sigmaAbs', '-ASCII');

ruta = [folderOutput,'eigenvectors.txt'];
save(ruta, 'eigenvectors', '-ASCII');

% Elimine first column with the counter of first column so we keep the eigenvectors as variables

eigenvectors(:,1) = [];

% Plot PCA results (eigenvalues)

disp(' ');
str = 'Do you want see results of PCA?? Y/N [Y]: ';
respuesta = input(str,'s');
respuesta = 'N';
if isempty(respuesta)
    respuesta = 'y';
end

if respuesta=='y'
    % Plot eigenvalues
    figure();
    bar(eigenvalues(1:pcaLong,1));
    title('Eigenvalues');
    text(9,eigenvalues(1,1), num2str(eigenvalues(1:pcaLong,:)),'HorizontalAlignment','left','VerticalAlignment','baseline');
    
%     % Plot sigma
%     figure();
%     bar(sigmaAbs(1:pcaLong,1));
%     title('Sigma');
%     text(9,sigmaAbs(1,1), num2str(sigmaAbs(1:pcaLong,:)),'HorizontalAlignment','left','VerticalAlignment','baseline');
    
    % Plot eigenvectors
    figure();
    subplot(5,1,1);
    title('Eigenvectors');
    for i= 1:5
        subplot(5,1,i);
        plot(eigenvectors(:,i));
    end
    subplot(5,1,1);
    title('Eigenvectors');
       
%     disp('Press any key to continue');
%     pause;
    
end


% Select number of species

% disp(' ');
num = input('Number of species [3]: ');
%num = 3;
if isempty(num)
    numeroEspecies = 3;
else
    numeroEspecies = num;
end


    close all;
%% Cut matrices

pointsAbs = length(Intensities);
pointsHoltzer = pointsAbs;
pointsKratky = pointsAbs;
pointsPorod = pointsAbs;


seleccion = 1;
plots;
disp(' ');
str = 'Do you want to cut the matrices?? Y/N [Y]: ';
cortar = input(str,'s');
% cortar = 'Y';
if isempty(cortar)
    cortar = 'y';
end


if (cortar == 'y' || cortar =='Y')
    
    if units == 'N'
        corteAbs = 3.0;
        corteHoltzer = 1.6;
        corteKratky = 1.6;
        cortePorod = 1.6;
    elseif units == 'A'
        corteAbs = 0.50;
        corteHoltzer = 0.30;
        corteKratky = 0.30;
        cortePorod = 0.16;
    end
    pointsAbs = 0;
        pointsHoltzer = 0;
        pointsKratky = 0;
        pointsPorod = 0;
       
        
        while pointsAbs == 0
            if valoresS(i,1) >= corteAbs
                pointsAbs = i;
            end
            i= i+1;
        end
        i=1;
        while pointsHoltzer == 0
            if valoresS(i,1) >= corteHoltzer
                pointsHoltzer = i;
            end
            i= i+1;
        end
        i=1;
        while pointsKratky == 0
            if valoresS(i,1) >= corteKratky
                pointsKratky = i;
            end
            i= i+1;
        end
        i=1;
        while pointsPorod == 0
            if valoresS(i,1) >= cortePorod
                pointsPorod = i;
            end
            i= i+1;
        end
    likeIt = 'n';
    close all;
    seleccion = 1;
    plots;
    if units == 'N'
        disp(' ');
        disp('Default (1/nm):');
        disp('Cut for Absolute scale  = 3.0');
        disp('Cut for Holtzer scale  = 1.6');
        disp('Cut for Kratky scale = 1.6');
        disp('Cut for Porod scale  = 0.7');
    elseif units == 'A'
        disp(' ');
        disp('Default (1/???):');
        disp('Cut for Absolute scale  = 0.3');
        disp('Cut for Holtzer scale  = 0.16');
        disp('Cut for Kratky scale = 0.16');
        disp('Cut for Porod scale  = 0.07');
    end
    
%     str = 'Do you want default values?? [Y]/N: ';
%     likeIt = input(str,'s');
    likeIt = 'Y';

    if (isempty(likeIt) || likeIt == 'Y' || likeIt == 'y')
        likeIt = 'y';
        if units == 'N'
            corteAbs = 3.0;
            corteHoltzer = 1.6;
            corteKratky = 1.6;
            cortePorod = 0.7;
        elseif units == 'A'
            corteAbs = 0.5;
            corteHoltzer = 0.16;
            corteKratky = 0.16;
            cortePorod = 0.07;
        end
    else
        pointsAbs = length(Intensities);
        pointsHoltzer = pointsAbs;
        pointsKratky = pointsAbs;
        pointsPorod = pointsAbs;
    end
    

    while likeIt ~= 'y'
    
        corteAbs = input('Cut for Absolute scale: ');
        corteHoltzer = input('Cut for Holtzer scale: ');   
        corteKratky = input('Cut for Kratky scale: ');
        cortePorod = input('Cut for Porod scale: ');
        
        pointsAbs = 0;
        pointsHoltzer = 0;
        pointsKratky = 0;
        pointsPorod = 0; 
        while pointsAbs == 0
            if valoresS(i,1) >= corteAbs
                pointsAbs = i;
            end
            i= i+1;
        end
        i=1;
        while pointsHoltzer == 0
            if valoresS(i,1) >= corteHoltzer
                pointsHoltzer = i;
            end
            i= i+1;
        end
        i=1;
        while pointsKratky == 0
            if valoresS(i,1) >= corteKratky
                pointsKratky = i;
            end
            i= i+1;
        end
        i=1;
        while pointsPorod == 0
            if valoresS(i,1) >= cortePorod
                pointsPorod = i;
            end
            i= i+1;
        end
        
        close all;
        seleccion = 1;
        plots;
        str = 'Do you want these values?? [Y]/N: ';
        likeIt = input(str,'s');
        
        if (isempty(likeIt) || likeIt == 'Y' || likeIt == 'y')
            likeIt = 'y';
        end
        
        
        
    end

end
close all;


%% Definition of initial estimations

% Run pure for doing an estimation

close all;

disp(' ');disp(' ');
disp('Calculating initial estimations');
disp(' ');

for i=1:length(OutFiles)
    fprintf('Curve %s - %s\n',num2str(i), OutFiles(i).name);
end
disp(' ');
disp(' ');

[pureAbs, curvasEstInic(1,:)] = pure(valoresIActual,numeroEspecies,5);

legendCell = ' ';
for i=1:numeroEspecies
    semilogy(valoresS, valoresIActual(:,curvasEstInic(i)));
    legendCell=[legendCell, '''curve ',num2str(curvasEstInic(i))];
    legendCell=[legendCell, ''''];
    if i<numeroEspecies
        legendCell=[legendCell, ','];
    end
    hold on; 
end
texto = ['legend(', legendCell,')'];
eval(texto);
title('Initial estimations');
% pause;
    
%     
% str = ['\nDo you want use the initial estimations selected by pure?? [', num2str(curvasEstInic), '] [Y]/N: '];
% respuesta = input(str,'s');
respuesta = 'Y';
if isempty(respuesta)
    respuesta = 'y';
    pureOK = 1;
end


close all;

correct = 'n';


if respuesta=='n'
    pureOK = 0;
    while correct == 'n'
        for i=1:length(OutFiles)
            fprintf('Curve %s - %s\n',num2str(i), OutFiles(i).name);
        end
        disp(' ');
        disp('Introduce the curves that you want as initial estimations: ');
        i=1;
        while i<=numeroEspecies
            k = input('Initial estimation: ');
            if k<=curvasTotales
                curvasEstInic(1,i) = k;
                i=i+1;
            else
                fprintf(2,'This curve does not exist!!! Must be equal or smaller than %s\n', num2str(curvasTotales));
            end
        end
        figure();
        legendCell = ' ';
        semilogy(valoresS, valoresIActual(:,curvasEstInic(i)));
        legendCell=[legendCell, '''curve ',num2str(curvasEstInic(i))];
        legendCell=[legendCell, ''''];
        if i<numeroEspecies
            legendCell=[legendCell, ','];
        end
        hold on;
        texto = ['legend(', legendCell,')'];
        eval(texto);
        title('Initial estimations');
        str = ['\nDo you want use the initial estimations selected by pure?? [', num2str(curvasEstInic), '] [Y]/N: '];
        correct = input(str,'s');
        if isempty(correct)
            correct = 'y';
        end
        if correct == 'n' || correct == 'N'
            close all;
        end
        
        
    end
else
    pureOK = 1;
end

% str = ['Do you want sort the curves?? [Y]/N: '];
% respuesta = input(str,'s');
respuesta = 'Y';
if isempty(respuesta)
    respuesta = 'y';
end
if respuesta=='y' || respuesta == 'Y'
    curvasEstInic = sort(curvasEstInic);
%     InitEstKratky = sort(InitEstKratky);
%     InitEstHoltzer = sort(InitEstHoltzer);
%     InitEstPorod = sort(InitEstPorod);
end
% Initaverage = zeros(numeroEspecies,length(valoresIActual));
% for i = 1:numeroEspecies
%     Initaverage(i,:) = (valoresIActual(:,curvasEstInic(1,i))+valoresIActual(:,InitEstKratky(1,i))+valoresIActual(:,InitEstHoltzer(1,i))+valoresIActual(:,InitEstPorod(1,i)))/4;
% end

%% Select constraints

numMatrices = 2;
nsign = numeroEspecies;

tipoALS = 0;

while (tipoALS ~= 1) && (tipoALS ~= 2) && (tipoALS ~= 3)
    disp('');
    disp('');
    disp('Type of experiment');
    disp('(1) Fibrillation or folding (default)');
    disp('(2) Titration');
    disp('(3) SEC-SAXS experiment');
%     tipoALS = input('Select type of experiment: ');
    tipoALS = 1;
    if isempty(tipoALS)
        tipoALS = 1;
    end
end

if(tipoALS == 1)
    fprintf('\n');
    disp('Standard options for fibrillation or folding:');
    disp('- Non-negativity (fnnls) for concentrations and spectra');
    disp('- Closure for concentrations (equal to 1.0)');
    disp('Press any key to continue');
    %pause;
    
end

if(tipoALS == 2)
    fprintf('\n');
    disp('Standard options for titration:');
    disp('- Non-negativity (fnnls) for concentrations and spectra');
    disp('- Unimodality for concentrations');
    disp('- No closure. Normalization of spectra (same length)');
    %pause;
end

if(tipoALS == 3)
    disp('\n\nStandard options for SEC-SAXS experiment:');
    disp('- Non-negativity (fnnls) for concentrations and spectra');
    disp('- No closure. No Normalization');
    %pause;
    
    moduloSECSAXS;
end

% Equality constrain (optional)

% fijar = input('\nDo you want to fix any species? Y/[N] ','s');
% if isempty(fijar)
    fijar = 'n';
    curvasFijadas = 0;
% end

if (fijar == 'y' || fijar == 'Y')
    j=1;
    for i=1:numeroEspecies
        fprintf('\nInitial estimation #%s - curve %s (%s)\n', num2str(i), num2str(curvasEstInic(i)), OutFiles(curvasEstInic(i)).name);
%         fijarYN = input('Do you want to fix this species? Y/[N] ','s');
        fijarYN = 'N';
        if isempty(fijarYN)
            fijarYN = 'n';
        end
        if (fijarYN == 'y' || fijarYN == 'Y')
            curvasFijadas = 1;
            curvasFijadas2(i) = 1;
            curvasFijadas3(j) = curvasEstInic(i);
            j = j+1;
        else
            curvasFijadas2(i) = 0;
        end
    end
%     curvasFijadas = input('Which curve(s) do you want to fix (e.g. sp 1 and 12, [1 12])?? [0 for none, by default]: ');
    if isempty(curvasFijadas)
        curvasFijadas = 0;
    end
else
 curvasFijadas = 0;
end


% Convergence criterion

disp(' ');
num = input('Convergence criterion (0.1% by default): ');
if isempty(num)
    tolsigma = 0.1;
else
    tolsigma = num;
end

% Number of iterations

disp(' ');
num = input('Maximum number of iterations (50 by default): ');
if isempty(num)
    numIterations = 1000;
else
    numIterations = num;
end

% GRAPHICAL OUTPUT OPTIONAL
% gr=input('\nDo you want graphic output during the ALS optimization [Y]/N)? ','s');
gr = 'y';
if gr=='Y' | gr =='y'
    gr='y';
end
if isempty(gr)
    gr = 'y';
end

close all;

% % CLEAR
% 
% listaVar = who;
% 
% close all;

%% Import Closure Restraints
InFilesList = [];
pattern = '*.txt';
for a = 1:numexpv
    filePattern = fullfile(foldernames{a}, pattern);
    InFiles = dir(filePattern);
    InFilesList = [InFilesList; InFiles];
end
ConcConsAll = [];
for k = 1:length(InFilesList)
    basefilename = InFilesList(k).name;
    %foldername = InFilesList(k).folder;
    foldername = foldernames{k};
    fullFileName = fullfile(foldername,basefilename);
    delimiter = ' ';
    formatSpec = '%f%f%f%[^\n\r]';
    fileID = fopen(fullFileName,'rt');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);
    fclose(fileID);
    ConcCons = [dataArray{1:end-1}]; 
    ConcConsAll = vertcat(ConcConsAll, ConcCons);
    vclos1((ExperimentRanges(k,1)+1):ExperimentRanges(k+1,1),k) = sum(ConcCons,2);
    clearvars filename delimiter formatSpec fileID dataArray ans;
end


%% ALS optimization (workflow)

close all;


combinacionesWF = [1 0 0 0; % A
                    1 1 0 0; % A+R
                    1 0 1 0; % A+K
                    1 0 0 1; % A+P   
                    1 1 1 0; % A+R+K
                    1 1 0 1; % A+R+P
                    1 0 1 1; % A+K+P
                    1 1 1 1]; % A+K+R+P
          
                
                
matricesUsadasWF = combinacionesWF;      

clear -global statisticsWorkflow;
clear -global coptWorkflow;
clear -global speciesWorkflow;
clear -global curvasEliminadasWF;
clear -global estimInicialesWF;
clear -global curvasEstimInicWF;
clear -global reconstWorkflow;
clear -global fitWorkflow;
clear -global coptElimWorkflow;
clear -global diferenciaXiWF;

field = 'Test';

global statisticsWorkflow;
statisticsWorkflow = struct('lackOfFit_PCA', 1, 'lackOfFit_Exp', 1, 'percentR2', 1, 'iteration', 1,  'solution',1);
global coptWorkflow;
coptWorkflow = struct(field, 1);
global speciesWorkflow;
speciesWorkflow = struct(field, 1);
global curvasEliminadasWF;
curvasEliminadasWF = struct(field, 1);
global estimInicialesWF;
estimInicialesWF = struct(field, 1);
global curvasEstimInicWF;
curvasEstimInicWF = struct(field, 1);
global reconstWorkflow;
reconstWorkflow = struct(field, 1);
global fitWorkflow;
fitWorkflow = struct(field, 1);
global coptElimWorkflow;
coptElimWorkflow = struct(field, 1);




test = 1;
testActual = 1;
continuar = 'y';
curvasEliminadas = 0;
innn = 1;



while (continuar=='y');
    
    % ALS optimization for all combinations
    
    for test=testActual:((testActual+length(combinacionesWF(:,1)))-1)
        combinacionActual = matricesUsadasWF(test,:);
        matrizAbsolute = valoresIActual';
        
        if exist('ssel', 'var') clear ssel, end;

        
        %crearMat_avg;
        crearMat_manualpure_rank;
        
        comb = 'Absolute';
        if (matricesUsadasWF(test,2)==1),  comb = [comb, ' + Holtzer']; end; 
        if (matricesUsadasWF(test,3)==1),  comb = [comb, ' + Kratky']; end;
        if (matricesUsadasWF(test,4)==1),  comb = [comb, ' + Porod']; end;
        
        nexp = length(puntosMatrices)*numexp;
        isp = ones (numexp,numeroEspecies);
        csel = NaN (length(Intensities(1,:)),numeroEspecies);
        csel(1,:) = ConcConsAll(1,:)/scalefactor(1,1);
        csel(ExperimentRanges(2)+1,:) = ConcConsAll(ExperimentRanges(2)+1,:)/scalefactor(ExperimentRanges(2)+1,1);
        vclos1 = zeros(ExperimentRanges(3),2);
        for i = 1:ExperimentRanges(2)
            vclos1(i,1)= ConcConsAll(1,1)/scalefactor(i,1);
        end
        for i = ExperimentRanges(2)+1:ExperimentRanges(3)
            vclos1(i,2) = ConcConsAll(ExperimentRanges(2)+1,2)/scalefactor(i,1);
        end
        ssel = NaN(numeroEspecies,size(matrizInput,2));
        
        % Call autoALS.m (text version modified for automatize)
        if tipoALS == 1
            %[coptAUTO,soptAUTO, sdopt, r2opt, itopt, solution]=als_6_auto_csel(matrizInput,estimInicialesInput,nexp,puntosMatrices,numexp,ExperimentRanges,numIterations,tolsigma,gr,isp,csel,ssel);
            [coptAUTO,soptAUTO, sdopt, r2opt, itopt, solution]=als_closure2(matrizInput,estimInicialesInput,nexp,puntosMatrices,numexp,ExperimentRanges, numIterations,tolsigma,gr,isp,csel,ssel,vclos1,valoresS,test);
        end
        if tipoALS == 2
            [coptAUTO,soptAUTO, sdopt, r2opt, itopt, solution]=autoALS2(matrizInput,estimInicialesInput,puntosMatrices, numIterations,tolsigma,comb, gr, ssel);
        end
        if tipoALS == 3
                
            [coptAUTO,soptAUTO, sdopt, r2opt, itopt, solution]=autoALS3(matrizInput,estimInicialesInput,puntosMatrices, numIterations,tolsigma,comb, gr, csel, ssel, expAbs, coeff);
                        
        end
        speciesCol = soptAUTO';
        species = soptAUTO;
        copt = coptAUTO;
        
               
        checkTest_scaled;
        
        
        if  curvasEliminadas ~= 0
           
            for ii = 1:length(chiSquareAll(:,1))
            
                diferenciaXiWF(ii,testActual) = chiSquareAll(ii,testActual) - chiSquareAll(ii,select);
            end
            diferenciaXiWF(curvasEliminadas(1,end),testActual) = 0;
            
        end
        
        
        if curvasEliminadas == 0
            curvasElimLogic(testActual,1) = 0 ;
        else
            curvasElimLogic(testActual,1) = 1 ;
        end
        
        
        
        testActual = testActual+1;
        
    end

    % Obtain smaller chi square from all test we have until now

    for i=(testActual-length(matricesUsadasWF(:,1))):length(matricesUsadasWF(:,1))
        if statisticsWorkflow(i).solution == 1
            temporalChiAverage(i,1) = chiAverageWorkflow(i,1);
        else
            temporalChiAverage(i,1) = 1e20;
        end
    end
    [minChi,combMenorXi] = min(temporalChiAverage);

    % Show solutions: number of test, combination and chi square
    
    disp(' ');
    for il=(testActual-length(matricesUsadasWF(:,1))):length(matricesUsadasWF(:,1))
        comb2 = 'A';
        if (matricesUsadasWF(il,2)==1),  comb2 = [comb2,'+H']; end; 
        if (matricesUsadasWF(il,3)==1),  comb2 = [comb2,'+K']; end;
        if (matricesUsadasWF(il,4)==1),  comb2 = [comb2,'+P']; end;

        fprintf('Test %d: xi2 (%s) = %4.2f', il, comb2 , chiAverageWorkflow(il,1));
        if statisticsWorkflow(il).solution == 1
            fprintf('\n');
        end
        if statisticsWorkflow(il).solution == 2
            fprintf(2,' DIVERGENCE!\n');
        end
        if statisticsWorkflow(il).solution == 3
            fprintf(2,' NUMBER OF ITERATIONS EXCEEDED THE ALLOWED!\n');
        end
    end
    clear il;

    % Plot solutions
    
    close all
     str = 'Do you want to see the plots?? Y/N [Y]: ';
     respuesta = input(str,'s');
     if isempty(respuesta)
         respuesta = 'y';
     end
     if respuesta =='y'
         close all;
        for il=(testActual-length(matricesUsadasWF(:,1))):length(matricesUsadasWF(:,1))
            comb = 'Absolute';
            if (matricesUsadasWF(il,2)==1),  comb = [comb, ' + Holtzer']; end; 
            if (matricesUsadasWF(il,3)==1),  comb = [comb, ' + Kratky']; end;
            if (matricesUsadasWF(il,4)==1),  comb = [comb, ' + Porod']; end;
            figure();
            
            
            
            if curvasElimLogic(il,1) == 0
                subplot(4,1,1:2);
            else
                subplot(5,1,1:2);
            end
            semilogy(speciesWorkflow(il).Test(1:pointsAbs,:));
            texto = ['xi2 = ', num2str(chiAverageWorkflow(il,1))];
            text(0.8, 0.8, texto, 'Units', 'normalized');
            title(comb);
            if statisticsWorkflow(il).solution ~= 1
                text(0.8, 0.7, 'DIVERGENCE!', 'Units', 'normalized','Color','Red');
            end
            
            if curvasElimLogic(il,1) == 0
                subplot(4,1,3:4);
            else
                subplot(5,1,3:4);
            end
            plot(coptElimWorkflow(il).Test);
            title('Concentrations');
            
            
            if curvasElimLogic(il,1) ~= 0
                subplot(5,1,5);
                bar(diferenciaXiWF(:,il));
            end
            
           
            
        end
        fprintf('\n\n');
        disp('Press any key to continue');
        %pause;  
        
        close all;
        clear il; 
%     end

    % Select best combination

    comb = 'Absolute';
    if (matricesUsadasWF(combMenorXi,2)==1),  comb = [comb, ' + Holtzer']; end; 
    if (matricesUsadasWF(combMenorXi,3)==1),  comb = [comb, ' + Kratky']; end;
    if (matricesUsadasWF(combMenorXi,4)==1),  comb = [comb, ' + Porod']; end;
    disp(' ');
    disp(['The best combination is ', comb, ' (xi2 = ', num2str(chiAverageWorkflow(combMenorXi,1)), ')']);
    disp(' ');

%     str = 'Do you want use this combination? Y/N [Y]: ';
%     respuesta = input(str,'s');
    respuesta = 'Y';
    if isempty(respuesta)
        respuesta = 'y';
    end
    if respuesta == 'y'
        select = combMenorXi;
    else
        disp(' ');
        for il=(testActual-length(matricesUsadasWF(:,1))):length(matricesUsadasWF(:,1))
            comb2 = 'A';
            if (matricesUsadasWF(il,2)==1),  comb2 = [comb2,'+H']; end; 
            if (matricesUsadasWF(il,3)==1),  comb2 = [comb2,'+K']; end;
            if (matricesUsadasWF(il,4)==1),  comb2 = [comb2,'+P']; end;

            fprintf('Test %d (%s): xi2  = %4.2f', il, comb2 , chiAverageWorkflow(il,1));
            if statisticsWorkflow(il).solution == 1
                fprintf('\n');
            end
            if statisticsWorkflow(il).solution == 2
                fprintf(2,' DIVERGENCE!\n');
            end
            if statisticsWorkflow(il).solution == 3
                fprintf(2,' NUMBER OF ITERATIONS EXCEEDED THE ALLOWED!\n');
            end
        end
        clear il;

        disp(' ');
%         select = input('Which combination do you want?? ');
        select = 1;
    end
    
    
    
    
    % Clean data set

    disp(' ');
    disp('Now we going to clean the data set...');
    disp(' ');

    % Show the chi of the combination selected     
    
    temp = chiSquareAll(:,select);
    chiTemp = chiAverageWorkflow(select,1);

    fprintf('Average xi2: %4.2f', chiTemp);
    if statisticsWorkflow(select).solution == 1
        fprintf('\n');
    end
    if statisticsWorkflow(select).solution == 2
        fprintf(2,' DIVERGENCE!\n');
    end
    if statisticsWorkflow(select).solution == 3
        fprintf(2,' NUMBER OF ITERATIONS EXCEEDED THE ALLOWED!\n');
    end
    
       
    % Select curves with bigger chi square from the selected combination
    % (first 10)
    
    [mayorXi, curvaMayorXi] = max(chiSquareAll(:,select));
    
    disp(' ');
    disp('The 10 worst xi2 are:');
    worstXi = sort(chiSquareAll(:,select),'descend');
    disp(worstXi(1:10,1));
    disp(['The curve with worst xi2 (', num2str(curvaMayorXi), ') - xi2 = ', num2str(mayorXi)]);
    
%     str = 'Do you want remove this curve? Y/N [Y]: ';
%     eliminar = input(str,'s');
    eliminar = 'n';
    if isempty(eliminar)
        eliminar = 'y';
    end
    if (eliminar=='y')
        coincide = find(curvaMayorXi==curvasEstInic);
        coincide2 = isempty(coincide);
        if (coincide2 ~= 1) 
            fprintf(2,' THIS CURVE IS AN INITIAL ESTIMATION!\n');
            
            str = 'Do you still want eliminate this curve? Y/N [Y]: ';
            eliminarC = input(str,'s');
            if isempty(eliminarC)
                eliminarC = 'y';
            end
            if (eliminarC=='y')
                disp('Curve will be remove from data set');
                curvasEliminadas(1,innn) = curvaMayorXi;
                innn = innn+1;
                valoresIActual(:,curvaMayorXi) = 0;
                disp(' ');
                disp('Now we going to select again the initial estimation without that curve');
                
                % Run pure again with the clean data set
                [pureAbs, curvasEstInic(1,:)] = pure(valoresIActual,numeroEspecies,5); 

                str = ['Do you want use the initial stimations selected by pure?? [', num2str(curvasEstInic), '] Y/N [Y]: '];
                respuesta = input(str,'s');
                if isempty(respuesta)
                    respuesta = 'y';  
                    
                end

                if respuesta=='n'
                    disp('Introduce the curves that you want as initial estimations: ');
                    i=1;
                    while i<=numeroEspecies
                       k = input('Initial estimation: '); 
                       if k<=(curvasTotales-length(curvasEliminadas))
                           curvasEstInic(1,i) = k;
                           i=i+1;
                       else
                           fprintf(2,'This curve does not exist!!! Must be equal or smaller than %s\n', num2str(curvasTotales));
                       end    
                    end
                end 
            else
                % If we don't eliminate anymore curves we don't repeat the
                % combinations
                continuar = 'n';
            end
        else
            curvasEliminadas(1,innn) = curvaMayorXi;
            innn = innn+1;
            valoresIActual(:,curvaMayorXi) = 0;
        end
    
        
    % If we remove some curve, repeat all the combinations without the
    % curve
        close all;
        
    fprintf('\n\n');
    disp('Now we going to repeat the different combinations without this curve');
    disp('Press any key to continue');
    %pause;
    
    test = test +1;
    matricesUsadasWF((test:(test+length(combinacionesWF))-1),:) = combinacionesWF;
    
    
    
    else
        % If we don't eliminate anymore curves we don't repeat the
        % combinations
        continuar = 'n';
        
        
        
    end
    
    
        
end % end Workflow


%% Save INFO FILE


rutaMCRALS = [folderOutput, 'info.txt'];
fileID = fopen(rutaMCRALS,'w');

fprintf(fileID, '\nInput files folder: %s\n\n\n', folderOutput);

fprintf(fileID, 'LIST OF FILES\n');

for i=1:length(OutFiles)
    fprintf(fileID, 'Curve %s - %s\n',num2str(i), OutFiles(i).name);
end


fprintf(fileID,'\n\nNumber of species: %2.0f\n\n', numeroEspecies);

fprintf(fileID,'\n\n');
fprintf(fileID,'CONSTRAINTS\n');


if tipoALS == 1
    fprintf(fileID,'\nType of experiment: Fibrillation or folding\n');
    fprintf(fileID,'\n* Non-negativity: Concentrations and spectra');
    fprintf(fileID,'\n\t- Selected algorithm for concentrations: fnnls (all species)');
    fprintf(fileID,'\n\t- Selected algorithm for spectra: fnnls (all species)');
    fprintf(fileID,'\n\n* Closure: Concentrations');
    fprintf(fileID,'\n\t- Number of closure constants: 1');
    fprintf(fileID,'\n\t- Closure constant: 1.0');    
    fprintf(fileID,'\n\t- Closure condition: Equal condition (all species)');
      
    if curvasEliminadas ~= 0
        fprintf(fileID,'\n* Fixed curves: \n');
        for i=1:length(curvasEliminadas)
            fprintf(fileID,'\t* Curve %d\n', curvasEliminadas(i));
        end 
    end
end

if(tipoALS == 2)
    fprintf(fileID,'\nType of experiment: Titration\n');
    fprintf(fileID,'\n* Non-negativity: Concentrations and spectra');
    fprintf(fileID,'\n\t- Selected algorithm for concentrations: fnnls (all species)');
    fprintf(fileID,'\n\t- Selected algorithm for spectra: fnnls (all species)');
    fprintf(fileID,'\n* Unimodality for concentrations');
    fprintf(fileID,'\n* No closure. Normalization of spectra (same length)');
end

fprintf(fileID,'\n\nConvergence criterion: %.1f\n', tolsigma);
fprintf(fileID,'\nMaximum number of iterations: %d\n\n', numIterations);






%% Save Workspace

% rutaMCRALS = [folderOutput, 'workflow.mat'];
% save(rutaMCRALS);
% 

%% Reconstruction

% str = 'Do you want reconstruct the curves?? [Y]/N: ';
% respuesta = input(str,'s');
respuesta = 'y';
if isempty(respuesta)
    respuesta = 'y';
end

if respuesta =='y'

    disp(' ');
%     reconstTest = input('Which combination do you want the reconstruction?? (i.e. [1 5 7]) - 0 for all: ');
    reconstTest = 0;
    if reconstTest == 0
        for i=1:length(matricesUsadasWF(:,1))
            reconstTest(i) = i;
        end
    end
  
    for indice=1:length(reconstTest)
        testActual = reconstTest(indice);
        comb2 = 'A';
        if (matricesUsadasWF(testActual,2)==1),  comb2 = [comb2,'H']; end;
        if (matricesUsadasWF(testActual,3)==1),  comb2 = [comb2,'K']; end;
        if (matricesUsadasWF(testActual,4)==1),  comb2 = [comb2,'P']; end;

        if testActual <= 9
            folder = ['Test0' , num2str(testActual), '_', comb2,'/'];
        else
            if testActual > 9
            folder = ['Test' , num2str(testActual), '_', comb2, '/'];
            end
        end
        
        rut = [folderOutput,folder, 'Reconstruction/'];
        mkdir(rut);

        for j = 1:curvasTotales  
            if(j<=9)
                rutaMCRALS = [rut, 'curveMCRALS_0' , num2str(j) , '.dat'];
                rutaExp = [rut, 'curveExp_0' , num2str(j) , '.dat'];
%                 rutaResiduals = [rut, 'residuals_0' , num2str(j) , '.txt'];
            else
                rutaMCRALS = [rut, 'curveMCRALS_' , num2str(j) , '.dat'];
                rutaExp = [rut, 'curveExp_' , num2str(j) , '.dat'];
%                 rutaResiduals = [rut, 'residuals_' , num2str(j) , '.txt'];
            end
            curvaMCRALSAbs(:,1) = valoresS(1:pointsAbs,1);
            curvaMCRALSAbs(:,2) = fitWorkflow(testActual).Test(:,j);
            curvaExpAbs(:,1) = valoresS(1:pointsAbs,1);
            curvaExpAbs(:,2) = Intensities(1:pointsAbs,j);
            curvaExpAbs(:,3) = Intensities(1:pointsAbs,j)*0.02;
            save(rutaMCRALS, 'curvaMCRALSAbs', '-ASCII');
            save(rutaExp, 'curvaExpAbs', '-ASCII');
        end 
    end
end
%[sband,cband,normband,tband,fband]=mcrbandsAUTO(coptAUTO,soptAUTO,nexp,csel);



%% Publish


publishReport_two;
% 
% 
%         nameReport = 'Report';
%    
%         for i=1:length(matricesUsadasWF(:,1));
%             publishTest(i) = i;
%         end
% 
%     rutapublish = [folderOutput, nameReport];
% 
%     publish('report.m','format','html', 'outputDir', rutapublish, 'showCode', false, 'createThumbnail', false);
% 
% 
%     close all;

    
    %% Monte Carlo

%montecarlo;
clear dir

