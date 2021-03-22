%% Save solutions


field = 'Test';


speciesWorkflow(test) = struct(field,speciesCol);


coptWorkflow(test) = struct(field,copt);


statisticsWorkflow(test) = struct('lackOfFit_PCA', sdopt(1,1), 'lackOfFit_Exp', sdopt(1,2), 'percentR2', r2opt, 'iteration', itopt,  'solution', solution);
curvasEliminadasWF(test) = struct(field,curvasEliminadas);


estimInicialesWF(test) = struct(field,estimInicialesInput);


curvasEstimInicWF(test) = struct(field,curvasEstInic);


%% Create output folder

comb2 = 'A';
if (matricesUsadasWF(test,2)==1),  comb2 = [comb2,'H']; end; 
if (matricesUsadasWF(test,3)==1),  comb2 = [comb2,'K']; end;
if (matricesUsadasWF(test,4)==1),  comb2 = [comb2,'P']; end;



if test <= 9
    folder = ['Test0' , num2str(test), '_', comb2,'/'];

else
    if test > 9
    folder = ['Test' , num2str(test), '_', comb2, '/'];

    end
end

dir = [folderOutput,folder];
mkdir(dir);

%% Save curves of optimized spectra with error (.dat)

% MCRALS does not give error. The error bar will be 3% of the Intensity


for index=1:numeroEspecies
    texto = ['species', num2str(index), 'Abs(:,1) = transpose(SAbsolute);'];
    eval(texto);
    texto = ['species', num2str(index), 'Abs(:,2) = speciesCol(1:puntosAbs,index);'];
    eval(texto);
    
    if columns >=3
        
    if tipoALS == 1
        DivisonFactor = matrizAbsolute(1,1)/pureAbs(1,1);
        texto = ['species', num2str(index), 'Abs(:,3) = valoresE(1:puntosAbs,curvasEstInicAbs(1,index))/DivisonFactor;'];
        eval(texto);
    end
    
%calculate fake error bars for the solutions


ratioErrores(1:puntosAbs,index) = (valoresE(1:puntosAbs,curvasEstInicAbs(1,index))*100)./(valoresI(1:puntosAbs,curvasEstInicAbs(1,index)));
    plot(ratioErrores(:,index));
    hold on;
    
    if tipoALS == 2
        for nn=1:puntosAbs
            texto = ['species', num2str(index), 'Abs(nn,3) = speciesCol(nn,index).* ratioErrores(nn,index);'];
%             disp('HOLAAAA');
            eval(texto);
        end
    end
        eval(texto);
    rutaMCRALS = [dir,'species', num2str(index), '.dat'];
    texto = ['save(rutaMCRALS, ''species', num2str(index), 'Abs'' , ''-ASCII'');'];
    eval(texto);
    
    end
end
RE = 1e-6;

% If Holtzer combination exists we compare Kratky solution with Absolute scale in order to compare if the solution is the same

if (matricesUsadasWF(1,2) == 1)
    for index=1:numeroEspecies
        texto = ['species', num2str(index), 'Holtzer(:,1) = transpose(SHoltzer);'];
        eval(texto);
        texto = ['species', num2str(index), 'Holtzer(:,2) = speciesCol(puntosAbs+1:(puntosAbs+puntosHoltzer),index);'];
        eval(texto);
        texto = ['species', num2str(index), 'Holtzer(:,2) = species', num2str(index), 'Holtzer(:,2)/escaladoHoltzer;'];
        eval(texto);
        texto = ['species', num2str(index), 'Holtzer(:,2) = species', num2str(index), 'Holtzer(:,2)./transpose(SHoltzer(1,:));']; 
        eval(texto);
        error=0;
        for indice=1:puntosHoltzer
            texto = ['esIgual2(indice,1) = species', num2str(index), 'Abs(indice,2) - species', num2str(index), 'Holtzer(indice,2);'];
            eval(texto);
            if abs(esIgual2(indice,1)) <= RE
                comprob(indice,1) = 1;
            else
                comprob(indice,1) = 0;
            end
        end
        numberOfZeros = nnz(comprob==0);
        if numberOfZeros>=1
            error = 1;
            fprintf(2,'WARNING!!! TEST %s. The optimized species %s in Holtzer matrix is different than solution from Absolute matrix\n',num2str(test), num2str(index));
        end
    end
    if error == 1
        disp(' ');
        str = 'The solutions are not the same, do you want to save the solution for the Holtzer matrix (absolute scale)?? Y/N [Y]: ';
        respuesta = input(str,'s');
        if isempty(respuesta)
            respuesta = 'y';
        end
        for index=1:numeroEspecies
            rutaMCRALS = [dir,'species', num2str(index), 'Holtzer.dat'];
            texto = ['save(rutaMCRALS, ''species', num2str(index), 'Holtzer'' , ''-ASCII'');'];
            eval(texto);
        end
    end
end


% If Kratky combination exists we compare Kratky solution with Absolute scale in order to compare if the solution is the same

if (matricesUsadasWF(1,3) == 1)
    if(matricesUsadasWF(1,2) == 1)
        suma = puntosAbs + puntosHoltzer;
    else
        suma = puntosAbs;
    end
    for index=1:numeroEspecies
        texto = ['species', num2str(index), 'Kratky(:,1) = transpose(SKratky);'];
        eval(texto);
        texto = ['species', num2str(index), 'Kratky(:,2) = speciesCol(suma+1:(suma+puntosKratky),index);'];
        eval(texto);
        texto = ['species', num2str(index), 'Kratky(:,2) = species', num2str(index), 'Kratky(:,2)/escaladoKratky;'];
        eval(texto);
        texto = ['species', num2str(index), 'Kratky(:,2) = species', num2str(index), 'Kratky(:,2)./transpose(SKratkyCuad(1,:));']; 
        eval(texto);
        error=0;
        for indice=1:puntosKratky
            texto = ['esIgual2(indice,1) = species', num2str(index), 'Abs(indice,2) - species', num2str(index), 'Kratky(indice,2);'];
            eval(texto);
            if abs(esIgual2(indice,1)) <= RE
                comprob(indice,1) = 1;
            else
                comprob(indice,1) = 0;
            end
        end
        numberOfZeros = nnz(comprob==0);
        if numberOfZeros>=1
            fprintf(2,'WARNING!!! TEST %s. The optimized species %s in Holtzer matrix is different than solution from Absolute matrix\n',num2str(test), num2str(index));
            error = 1;
        end
    end
    if error == 1
        disp(' ');
        str = 'The solutions are not the same, do you want to save the solution for the Kratky matrix (absolute scale)?? Y/N [Y]: ';
        respuesta = input(str,'s');
        if isempty(respuesta)
            respuesta = 'y';
        end
        for index=1:numeroEspecies
            rutaMCRALS = [dir,'species', num2str(index), 'Kratky.dat'];
            texto = ['save(rutaMCRALS, ''species', num2str(index), 'Kratky'' , ''-ASCII'');'];
            eval(texto);
        end
    end
end

% If Porod combination exists we compare Kratky solution with Absolute scale in order to compare if the solution is the same

if (matricesUsadasWF(1,4) == 1)
    suma = length(species) - puntosPorod;    
    for index=1:numeroEspecies
        texto = ['species', num2str(index), 'Porod(:,1) = transpose(SPorod);'];
        eval(texto);
        texto = ['species', num2str(index), 'Porod(:,2) = speciesCol(suma+1:(suma+puntosPorod),index);'];
        eval(texto);
        texto = ['species', num2str(index), 'Porod(:,2) = species', num2str(index), 'Porod(:,2)/escaladoPorod;'];
        eval(texto);
        texto = ['species', num2str(index), 'Porod(:,2) = species', num2str(index), 'Porod(:,2)./transpose(SPorodCuat(1,:));']; 
        eval(texto);
        error=0;
        for indice=1:puntosPorod
            texto = ['esIgual2(indice,1) = species', num2str(index), 'Abs(indice,2) - species', num2str(index), 'Porod(indice,2);'];
            eval(texto);
            if abs(esIgual2(indice,1)) <= RE
                comprob(indice,1) = 1;
            else
                comprob(indice,1) = 0;
            end
        end
        numberOfZeros = nnz(comprob==0);
        if numberOfZeros>=1
            error = 1;
            fprintf(2,'WARNING!!! TEST %s. The optimized species %s in Porod matrix is different than solution from Absolute matrix\n',num2str(test), num2str(index));
        end
    end
    if error == 1
        disp(' ');
        str = 'The solutions are not the same, do you want to save the solution for the Porod matrix (absolute scale)?? Y/N [Y]: ';
        respuesta = input(str,'s');
        if isempty(respuesta)
            respuesta = 'y';
        end
        for index=1:numeroEspecies
            rutaMCRALS = [dir,'species', num2str(index), 'Porod.dat'];
            texto = ['save(rutaMCRALS, ''species', num2str(index), 'Porod'' , ''-ASCII'');'];
            eval(texto);
        end
    end
end


%% Save concentration profile (.txt)
   
clear coptElimin;
cont = 1;
for i=1:length(copt(:,1))
%     coptElimin(cont,1) =i;
    if (any(curvasEliminadas == i)) == 0
        for j=1:numeroEspecies;
            coptElimin(cont,j) =copt(i,j); 
        end
        cont = cont + 1;
    end
end
clearvars ScaledCopt; 
% The file that we save have a first column with a counter
rutaMCRALS = [dir,'concentration.txt'];
save(rutaMCRALS, 'coptElimin', '-ASCII');
ScaledConc = [dir, 'ScaledConc.txt'];
ScaledCopt(:,1) = copt(:,1).*scalefactor(:,1);
ScaledCopt(:,2) = copt(:,2).*scalefactor(:,1);
ScaledCopt(:,3) = copt(:,3).*scalefactor(:,1);
save(ScaledConc, 'ScaledCopt', '-ASCII');
    
% We create a copt structure where we remove the curves so we don't plot it

% global coptElimWorkflow;
% if ~isempty(coptElimWorkflow) && isstruct(coptElimWorkflow)
%    coptElimWorkflow(end+1) = struct(field, 2);% make sure compatibility
% else % initialize
%    coptElimWorkflow = struct(field, 1);
% end
coptElimWorkflow(test) = struct(field,coptElimin);



%% Calculate pair distribution and save

% dirPR = [dir, 'AutoGnom/'];
% mkdir(dirPR);
% 
% 
% for index=1:numeroEspecies
%     filePR = [dir,'species', num2str(index)];
%     [pR, pRNorm, autoRgTemp, autoI0Temp, gnomDmaxTemp, textoRg, autoGnom, stat] = autoPR(filePR);  
%     texto = [filePR, 'pR = struct(''pR'', pR, ''pRNorm'', pRNorm, ''autoRg'', autoRgTemp, ''autoI0'', autoI0Temp, ''gnomDmax'', gnomDmaxTemp);'];
%     filePR = [dirPR,'species', num2str(index)];
%     savePR = [filePR, '_pR.txt'];
%     save(savePR, 'pR', '-ASCII');
%     savePR = [filePR, '_pRNorm.txt'];
%     save(savePR, 'pRNorm', '-ASCII');
%     savePR = [filePR, '_autoGnom.txt'];
%     fileID = fopen(savePR,'w');
%     fprintf(fileID,'Rg \t I(0) \t Dmax\n\n');
%     fprintf(fileID,'%.2f \t %.0f \t %.2f\n', autoRgTemp, autoI0Temp, gnomDmaxTemp); 
%     
%     
%     
% 
% 
%     
%     texto = ['pR_specie', num2str(index), ' = pRNorm;'];
%     eval(texto);
%     
%     texto = ['textRg', num2str(index), ' = textoRg;'];
%     eval(texto);
%     
%     texto = ['textGnom', num2str(index), ' = autoGnom;'];
%     eval(texto);
%     
%     autoRg(index) = autoRgTemp;
%     autoI0(index) = autoI0Temp;
%     gnomDmax(index) = gnomDmaxTemp(1);
% 
% 
%     filePR = [dir,'species', num2str(index)];
%     
%     filePR = [dir,'species', num2str(index), '.out']
%     filePR2 = [dirPR,'species', num2str(index), '.out'];
% 
%     texto = ['system(''mv ',filePR, ' ' , filePR2 ''');']
%     eval(texto);
%     
%     
% end




%% Reconstruct curves

% Absolute scale

[reconstruccionAbs] = reconstCurvas(copt,species(:,1:puntosAbs));

reconstWorkflow(test) = struct(field,reconstruccionAbs);


%% Calcule chi square

if columns >= 3
for index =1: curvasTotales;  
    reconsTemp(:,1) = valoresS(1:puntosAbs,1);
    reconsTemp(:,2) = reconstruccionAbs(1:puntosAbs,index);
    
    curvaExpTemp(:,1) = valoresS(1:puntosAbs,1);
    curvaExpTemp(:,2) = valoresI(1:puntosAbs,index);
    if columns >=3
    curvaExpTemp(:,3) = valoresE(1:puntosAbs,index);
    end
    
    [chiTemp, fitTemp] = compare2curves(curvaExpTemp, reconsTemp);
    
    chi(index,1) = chiTemp;
    fit(:,index) = fitTemp;  
end

if curvasEliminadas ~= 0
    for index = 1:length(curvasEliminadas)
        chi(curvasEliminadas(1,index)) = 0;
        fit(:,(curvasEliminadas(:,index))) = 0;
    end
    curvasTotalesElim =  length(matrizInput(:,1)) - length(curvasEliminadas);
else
        curvasTotalesElim =  curvasTotales;
end

fitWorkflow(test) = struct(field,fit);

clear chiAverage;
chiAverage = 0;

for index = 1:curvasTotales
    chiAverage = chiAverage + chi(index,1);    
end

chiAverage = chiAverage / curvasTotalesElim;
residuals = transpose(matrizAbsolute) - reconstruccionAbs;
chiAverageWorkflow(test,1) = chiAverage;
chiSquareAll(:,test) = chi(:,1);

for i=1:length(chiSquareAll(:,test))
    chiSquareTemp(i,1) =i;
    chiSquareTemp(i,2) =chiSquareAll(i,test); 
end

rutaMCRALS = [dir,'chiSquareCurves.txt'];
save(rutaMCRALS, 'chiSquareTemp' , '-ASCII');

end

 %% Create a fake curve
% 
% unos(:,1) = valoresS;
% unos(:,2) = 1;
% 
% rutaMCRALS = [dir, 'Test0', num2str(test), '.dat'];
% save(rutaMCRALS, 'unos', '-ASCII');
% 
% clear unos;

%% Statistics file

if solution == 1
    solucion = ['CONVERGENCE'];
else if solution == 2
        solucion = ['ITERATIONS EXCEEDED ALLOWED'];
    else if solution == 3
            solucion = ['DIVERGENCE!!'];
        end
    end
end

rutaMCRALS = [folderOutput,folder, 'infoTest' , int2str(test), '.txt'];
fileID = fopen(rutaMCRALS,'w');


stringMat = ['Absolute'];
if(matricesUsadasWF(test,2) == 1)
    stringMat = [stringMat, ' + Holtzer'];
end
if(matricesUsadasWF(test,3) == 1)
    stringMat = [stringMat, ' + Kratky'];
end
if(matricesUsadasWF(test,4) == 1)
    stringMat = [stringMat, ' + Porod'];
end
    
fprintf(fileID,'Matrices used: %s\n', stringMat);
fprintf(fileID,'\nInitial estimations ');
if pureOK==1
    fprintf(fileID,'(PURE)\n');
else
    fprintf(fileID,'(user)\n');
end
    for i=1:length(PureComponents)
        fprintf(fileID,'\t* Curve %d\n', curvasEstInicAbs(i));
    end 
% fprintf(fileID,'\n\n');
fprintf(fileID,'\nPoints removed at the beginning: ');
fprintf(fileID,'%4.0f\n', elimPoints);
fprintf(fileID,'\nTotal curves: ');
fprintf(fileID,'%4.0f\n', curvasTotalesElim);

fprintf(fileID,'\nCut matrices: \n');
if unidades == 'A' || unidades == 'a'
    fprintf(fileID,'\t- Absolute: %4.2f 1/A\n', corteAbs);
    if(matricesUsadasWF(test,2) == 1)
        fprintf(fileID,'\t- Holtzer: %4.2f 1/A\n', corteHoltzer);
    end
    if(matricesUsadasWF(test,3) == 1)
        fprintf(fileID,'\t- Kratky: %4.2f 1/A\n', corteKratky);
    end
    if(matricesUsadasWF(test,4) == 1)
        fprintf(fileID,'\t- Porod: %4.2f 1/A\n', cortePorod);
    end
end

if unidades == 'N' || unidades == 'n'
    fprintf(fileID,'\t- Absolute: %4.2f 1/nm\n', corteAbs);
    if(matricesUsadasWF(test,2) == 1)
        fprintf(fileID,'\t- Holtzer: %4.2f 1/nm\n', corteHoltzer);
    end
    if(matricesUsadasWF(test,3) == 1)
        fprintf(fileID,'\t- Kratky: %4.2f 1/nm\n', corteKratky);
    end
    if(matricesUsadasWF(test,4) == 1)
        fprintf(fileID,'\t- Porod: %4.2f 1/nm\n', cortePorod);
    end
end
% fprintf(fileID,'\n\n');
% fprintf(fileID,'CONSTRAINTS\n');

if tipoALS == 1
    fprintf(fileID,'\nType of experiment: Fibrillation or folding');
end
if tipoALS == 2
    fprintf(fileID,'\nType of experiment:  Titration');
end
    

if curvasEliminadas ~= 0
    fprintf(fileID,'\nRemoved curves: \n');
    for i=1:length(curvasEliminadas)
        fprintf(fileID,'\t* Curve %d\n', curvasEliminadas(i));
    end 
end

if curvasFijadas ~= 0
    fprintf(fileID,'\nFixed curves: \n');
    for i=1:length(curvasFijadas3)
        fprintf(fileID,'\t* Curve %d\n', curvasFijadas3(i));
    end 
end



fprintf(fileID,'\n\n\n');
fprintf(fileID,'RESULTS OF OPTIMIZATION\n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'%s\n\n', solucion);
fprintf(fileID,'Optimum at iteration %1.0f\n\n', itopt);
fprintf(fileID,'Lack of fit (PCA): %8.2e\n', sdopt(1,1));
fprintf(fileID,'Lack of fit (exp): %8.2e\n', sdopt(1,2));
fprintf(fileID,'Percent of variance explained: %5.2f %%\n', r2opt*100);
fprintf(fileID,'\n');
fprintf(fileID,'Chi square reduced: %5.2f\n', chiAverage(1,1));



fclose(fileID);




%% Save workspace
% 
% dir = [folderOutput,folder];
% rutaMCRALS = [dir, test2];
% save(rutaMCRALS);


%% MCR-BANDS (In progress)

% addpath ../../../../../Software_Scripts/MCR-ALS/Software/mcrbands/
