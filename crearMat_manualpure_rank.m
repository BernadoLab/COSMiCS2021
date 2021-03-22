 %% Remove curves

if curvasEliminadas ~= 0
    matrizInputElimin = matrizInput;
    matrizAbsoluteElimin = matrizAbsolute;
    %cont = 1;
    for i=1:length(matrizAbsolute(:,1))
        if (any(curvasEliminadas == i)) == 0
            matrizAbsoluteElimin(i,:) = matrizAbsolute(i,:);
           % cont = cont + 1;
        else
            if (any(curvasEliminadas == i)) == 1
                matrizAbsoluteElimin(i,:) = 0;
            end
        end
    end
    clear matrizAbsolute;
    matrizAbsolute = matrizAbsoluteElimin;
    clear matrizAbsoluteElimin;
end



%% Create individual matrices of each representation

for i=1:length(matrizAbsolute(:,1))
    matrizAbsoluteCortada(i,:) = matrizAbsolute(i,1:puntosAbs);
    SAbsolute = transpose(valoresS(1:puntosAbs,1));
    matrizHoltzerCortada(i,:) = matrizAbsolute(i,1:puntosHoltzer);
    SHoltzer = transpose(valoresS(1:puntosHoltzer,1));
    matrizKratkyCortada(i,:) = matrizAbsolute(i,1:puntosKratky);
    SKratky = transpose(valoresS(1:puntosKratky,1));
    SKratkyCuad = SKratky.*SKratky;
    matrizPorodCortada(i,:) = matrizAbsolute(i,1:puntosPorod);
    SPorod = transpose(valoresS(1:puntosPorod,1));
    SPorodCuat = SPorod.^4;
end

clear matrizAbsolute;
matrizAbsolute = matrizAbsoluteCortada;

for i=1:length(matrizHoltzerCortada(:,1))
    matrizHoltzer(i,:) = matrizHoltzerCortada(i,:) .* SHoltzer(1,:);
end
for i=1:length(matrizKratkyCortada(:,1))
    matrizKratky(i,:) = matrizKratkyCortada(i,:) .* SKratkyCuad(1,:);
end
for i=1:length(matrizPorodCortada(:,1))
    matrizPorod(i,:) = matrizPorodCortada(i,:) .* SPorodCuat(1,:);
end   
clear matrizAbsoluteCortada;
clear matrizHoltzerCortada;    
clear matrizKratkyCortada;
clear matrizPorodCortada;  
    
%% Scale matrices

% We do SVD of each matrix to obtain an adequate scale. We do with the
% clean and cutted data set (matrizKratky, matrizHoltzer y matrizPorod)

% SVD Holtzer
% for i=1:min(size(matrizHoltzer))
    [uHoltzer,sHoltzer,vHoltzer,xHoltzer, SigmaHoltzer] = pcarep(matrizHoltzer',min(size(matrizHoltzer)));
% end 
for i=1:length(sHoltzer)
    eigenvaluesHoltzer(i,1) = sHoltzer(i,i);
end

% SVD Kratky
% for i=1:min(size(matrizKratky))
    [uKrat,sKrat,vKrat,xKrat, SigmaKrat] = pcarep(matrizKratky',min(size(matrizKratky)));
% end
for i=1:length(sKrat)
    eigenvaluesKrat(i,1) = sKrat(i,i);
end

% SVD Porod
% for i=1:min(size(matrizPorod))
    [uPorod,sPorod,vPorod,xPorod, SigmaPorod] = pcarep(matrizPorod',min(size(matrizPorod)));
% end   
for i=1:length(sPorod)
    eigenvaluesPorod(i,1) = sPorod(i,i);
end

escaladoHoltzer = eigenvalues(1,1)/eigenvaluesHoltzer(1,1);
escaladoKratky = eigenvalues(1,1)/eigenvaluesKrat(1,1);
escaladoPorod = eigenvalues(1,1)/eigenvaluesPorod(1,1);


%% Create input matrix

if ((combinacionActual(1,1) == 1))
    matrizInput = matrizAbsolute;
    puntosMatrices = puntosAbs;
end
if ((combinacionActual(1,2) == 1))
    % Scale Holtzer
    matrizHoltzer = matrizHoltzer * escaladoHoltzer; 
    if ((combinacionActual(1,1) == 0))
        matrizInput = matrizHoltzer;
        puntosMatrices = puntosHoltzer;
    else
        matrizInput = [matrizInput matrizHoltzer];
        puntosMatrices = [puntosMatrices puntosHoltzer];
    end
end
if ((combinacionActual(1,3) == 1))
    % Scale Kratky
    matrizKratky = matrizKratky * escaladoKratky;
    matrizInput = [matrizInput matrizKratky];
    puntosMatrices = [puntosMatrices puntosKratky];
end

if ((combinacionActual(1,4) == 1))
    % Scale Porod
    matrizPorod = matrizPorod * escaladoPorod;
    matrizInput = [matrizInput matrizPorod];
    puntosMatrices = [puntosMatrices puntosPorod];
end

%% Create initial estimations dataset 
firstcurve = horzcat(valoresS, valoresI(:,1), valoresE(:,1));
for i = 1:length(valoresI(1,:))
curvetocompare = horzcat(valoresS,valoresI(:,i));
chiforcomp1(i,1) = i;
chiforcomp1(i,2) = compare2curves(firstcurve, curvetocompare);
end
[~,idx] = sort(chiforcomp1(:,2));
sortedchi1values = chiforcomp1(idx,:);
secondcurve = horzcat(valoresS, valoresI(:,(ExperimentRanges(2)+1)), valoresE(:,(ExperimentRanges(2)+1)));
for i = 1:length(valoresI(1,:))
curvetocompare = horzcat(valoresS, valoresI(:,i));
chiforcomp2(i,1) = i;
chiforcomp2(i,2) = compare2curves(secondcurve, curvetocompare);
end
[~,idx] = sort(chiforcomp2(:,2));
sortedchi2values = chiforcomp2(idx,:);

for i = 1:length(sortedchi1values)
curvenumber = sortedchi1values(i,1);
[~,idx]=min(abs(curvenumber-sortedchi2values(:,1)));
chirankavglist(i,1) = curvenumber;
chirankavglist(i,2) = (i+idx)/2;
end
[~,idx] = sort(chirankavglist(:,1));
sortedchirankavg = chirankavglist(idx,:);
n = 1;
chiforcompavg = (chiforcomp1 + chiforcomp2)/2;
% [M,I] = max (chiforcompavg);
% chicompcomb(:,1) = chiforcomp1;
% chicompcomb(:,2) = chiforcomp2;
% chistddev(:,1) = std(chicompcomb,0,2);
[M,I] = max (sortedchirankavg(:,2));
ThirdCurveNumber = I;

PureComponents = horzcat(1,(ExperimentRanges(2)+1),I);

curvasEstInicAbs = PureComponents;
curvasEstInicHoltzer = PureComponents;
curvasEstInicKratky = PureComponents;
curvasEstInicPorod = PureComponents;

for i=1:numeroEspecies
    if ((combinacionActual(1,1) == 1))
        estimInicialesAbs(i,:) = matrizAbsolute((curvasEstInicAbs(1,i)),:);
    end
    if ((combinacionActual(1,2) == 1))
        estimInicialesHoltzer(i,:) = matrizHoltzer((curvasEstInicHoltzer(1,i)),:);
    end 
    if ((combinacionActual(1,3) == 1))
        estimInicialesKratky(i,:) = matrizKratky((curvasEstInicKratky(1,i)),:);
    end   
    if ((combinacionActual(1,4) == 1))
        estimInicialesPorod(i,:) = matrizPorod((curvasEstInicPorod(1,i)),:);
    end   
end
if ((combinacionActual(1,1) == 1))
    estimInicialesInput = estimInicialesAbs;
end
if ((combinacionActual(1,2) == 1))
    estimInicialesInput = [estimInicialesInput estimInicialesHoltzer];
end 
if ((combinacionActual(1,3) == 1))
    estimInicialesInput = [estimInicialesInput estimInicialesKratky];
end
if ((combinacionActual(1,4) == 1))
    estimInicialesInput = [estimInicialesInput estimInicialesPorod];
end


%% Create Equality constrain for spectra (if exist)
clear ssel;
if curvasFijadas == 1
    
    for ll=1:numeroEspecies
        if curvasFijadas2(1,ll) == 0
            ssel(ll,1:length(matrizInput)) = NaN;
        end
        if curvasFijadas2(1,ll) == 1
            ssel(ll,:) = matrizInput(curvasEstInic(1,ll),:);
        end
    end
    

else
    ssel = 0;
end




%% Define number of matrices 

numMatrices = 0;

for i=1:length(combinacionActual(1,:))
    if combinacionActual(1,i) == 1
        numMatrices = numMatrices + 1;
    end
end



