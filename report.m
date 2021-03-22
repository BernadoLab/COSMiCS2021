

%% 1. OPTIONS
%
fprintf('\nInput files folder: %s\n\n\n', folderInput);

fprintf('LIST OF FILES\n');

for i=1:length(OutFiles)
    fprintf('Curve %s - %s\n',num2str(i), OutFiles(i).name);
end

fprintf('\n\nNumber of species: %2.0f\n\n', numeroEspecies);
fprintf( 'Total curves: ');
fprintf( '%4.0f\n\n', curvasTotales);
fprintf( 'Points removed at the beginning: ');
fprintf( '%4.0f\n\n', elimPoints);
fprintf('\n');
fprintf('Cut matrices: \n');
if unidades == 'A' || unidades == 'a'
    fprintf('\t- Absolute: %4.2f 1/A\n', corteAbs);
    if(matricesUsadasWF(test,2) == 1)
        fprintf('\t- Holtzer: %4.2f 1/A\n', corteHoltzer);
    end
    if(matricesUsadasWF(test,3) == 1)
        fprintf('\t- Kratky: %4.2f 1/A\n', corteKratky);
    end
    if(matricesUsadasWF(test,4) == 1)
        fprintf('\t- Porod: %4.2f 1/A\n', cortePorod);
    end
end

if unidades == 'N' || unidades == 'n'
    fprintf('\t- Absolute: %4.2f 1/nm\n', corteAbs);
    if(matricesUsadasWF(test,2) == 1)
        fprintf('\t- Holtzer: %4.2f 1/nm\n', corteHoltzer);
    end
    if(matricesUsadasWF(test,3) == 1)
        fprintf('\t- Kratky: %4.2f 1/nm\n', corteKratky);
    end
    if(matricesUsadasWF(test,4) == 1)
        fprintf('\t- Porod: %4.2f 1/nm\n', cortePorod);
    end
end
%
%% 2. CONSTRAINTS
%
fprintf('CONSTRAINTS\n');

if tipoALS ==1
    

    fprintf('\n* Non-negativity: Concentrations and spectra');
    fprintf('\n\t- Selected algorithm for concentrations: fnnls (all species)');
    fprintf('\n\t- Selected algorithm for spectra: fnnls (all species)');
    fprintf('\n\n* Closure: Concentrations');
    fprintf('\n\t- Number of closure constants: 1');
    fprintf('\n\t- Closure constant: 1.0');    
    fprintf('\n\t- Closure condition: Equal condition (all species)');
     
end

if(tipoALS == 2)
    fprintf('\n');
    disp('Standard options for titration:');
    disp('- Non-negativity (fnnls) for concentrations and spectra');
    disp('- Unimodality for concentrations');
    disp('- No closure. Normalization of spectra (same length)');
    pause;
end

% const = find(constrains.wcons == 4, 1);
% if ~isempty(const)
%     fprintf('\n\n* Equality in concentration profiles\n');
% end
% const = find(constrains.wcons == 5, 1); 
% if ~isempty(const)
%     fprintf('\n\n* Equality in spectra profiles. Fixed curves: %.0f\n', curvasFijadas);
% end

fprintf('\n\n');

fprintf('\n\n* Convergence criterion: %.1f %%\n', tolsigma);

fprintf('\n\n* Maximum number of iterations: %.0f\n', numIterations);





% colorscheme = othercolor('Blues9',length(curvasFibrilacion));
colorscheme = jet(length(curvasFibrilacion));

figure();
l = getframe(gcf);
for i=1:length(curvasFibrilacion)
    semilogy(valoresS(1:puntosAbs,1),matrizInput(i,1:puntosAbs)','Color',colorscheme(i,:));
    hold on;
end
title('Matriz input (Absolute values) in semilogarithmic scale');

figure();
l = getframe(gcf);
for i=1:length(curvasFibrilacion)
    plot(valoresS(1:puntosHoltzer,1),matrizHoltzer(i,1:puntosHoltzer)','Color',colorscheme(i,:));
    hold on;
end
title('Matriz input (Holtzer values)');

figure();
l = getframe(gcf);
for i=1:length(curvasFibrilacion)
    plot(valoresS(1:puntosKratky,1),matrizKratky(i,1:puntosKratky)','Color',colorscheme(i,:));
    hold on;
end
title('Matriz input (Kratky values)');



figure();
l = getframe(gcf);
for i=1:length(curvasFibrilacion)
    plot(valoresS(1:puntosPorod,1),matrizPorod(i,1:puntosPorod)','Color',colorscheme(i,:));
    hold on;
end
title('Matriz input (Porod values)');



% Plot eigenvalues
    figure();
    l = getframe(gcf);
    bar(eigenvalues(1:10,1));
    title('Eigenvalues');
    text(9,eigenvalues(1,1), num2str(eigenvalues(1:10,:)),'HorizontalAlignment','left','VerticalAlignment','baseline');
    
%     % Plot sigma
%     figure();
%     l = getframe(gcf);
%     bar(sigmaAbs(1:10,1));
%     title('Sigma');
%     text(9,sigmaAbs(1,1), num2str(sigmaAbs(1:10,:)),'HorizontalAlignment','left','VerticalAlignment','baseline');
    
    % Plot eigenvectors
    figure();
    l = getframe(gcf);
    subplot(5,1,1);
    title('Eigenvectors');
    for i= 1:5
        subplot(5,1,i);
        plot(eigenvectors(:,i));
    end
    subplot(5,1,1);
    title('Eigenvectors');

 %<sprintf('\L')>
    
%% 3. RESULTS


fprintf('Test \t Combination \t Rem.curves \t Init.estim. \t l.o.f (PCA) \t l.o.f (exp) \t Variance \t X^2 \t Convergence\n');
fprintf('-------------------------------------------------------------------------------------------------------------------------------------\n');
    
for il=1:length(matricesUsadasWF(:,1))
   
    comb2 = 'A';
    if (matricesUsadasWF(il,2)==1)   comb2 = [comb2,'+H']; end; 
    if (matricesUsadasWF(il,3)==1),  comb2 = [comb2,'+K']; end;
    if (matricesUsadasWF(il,4)==1),  comb2 = [comb2,'+P']; end;

    numMatrices = 0;

    for i=1:length(matricesUsadasWF(il,:))
        if matricesUsadasWF(il,i) == 1
            numMatrices = numMatrices + 1;
        end
    end
    
    if ((numMatrices == 1) || (numMatrices == 2) || (numMatrices == 3))
        comb2 = [comb2, '     '];
    end
%     if comb2 == 'A+K+R+P'; comb2 = 'A+K+R+P';end
    

    
    if (statisticsWorkflow(il).solution ==1)
        solucion = 'Convergence';
    end
    if (statisticsWorkflow(il).solution ==2)
        solucion = 'Divergence';
    end
    if (statisticsWorkflow(il).solution ==3)
        solucion = 'Num. iteration exceeded';
    end
    
    for in = length(curvasEstimInicWF(il).Test)
        textEstim = [num2str(curvasEstimInicWF(il).Test), ' '];
    end
    
    for in = length(curvasEliminadasWF(il).Test)
        textElim = [num2str(curvasEliminadasWF(il).Test), ' '];
    end
    
    fprintf('%d \t %s \t %s \t\t %s \t %f \t %f \t %f \t %.2f \t %s \n', il, comb2, textElim, textEstim, statisticsWorkflow(il).lackOfFit_PCA, statisticsWorkflow(il).lackOfFit_Exp, statisticsWorkflow(il).percentR2, chiAverageWorkflow(il,1), solucion);
    
    if(mod(il, length(combinacionesWF(:,1))) == 0)
        fprintf('\n');
    end
    
    
    clear comb;
end



%% 3. SOLUTIONS
%<p style="page-break-before: always">


for indice=1:length(publishTest)
    
    numeroTest = publishTest(indice);
    
    testActual = publishTest(indice);

    
    
    %%
    
    fprintf('\n\n\n\n\n\n');
    
    fprintf('**********\n');
    fprintf('* Test %d *\n', numeroTest);
    fprintf('**********\n\n');
    
    comb = ['Test ' num2str(testActual), ': Absolute'];
    if (matricesUsadasWF(testActual,2)==1),  comb = [comb, ' + Holtzer']; end; 
    if (matricesUsadasWF(testActual,3)==1),  comb = [comb, ' + Kratky']; end;
    if (matricesUsadasWF(testActual,4)==1),  comb = [comb, ' + Porod']; end;
    figure();
    subplot(2,1,1);
    semilogy(valoresS(1:puntosAbs),speciesWorkflow(testActual).Test(1:puntosAbs,:));
    texto = ['xi2 = ', num2str(chiAverageWorkflow(testActual,1))];
    text(0.8, 0.8, texto, 'Units', 'normalized');
    title(comb);
    if statisticsWorkflow(testActual).solution ~= 1
        text(0.8, 0.7, 'DIVERGENCE!', 'Units', 'normalized','Color','Red');
    end
    subplot(2,1,2);
    plot(coptElimWorkflow(testActual).Test);
    title('Concentrations');

    figure();
    plot(chiSquareAll(:,numeroTest), ':*');
    
    if curvasEliminadasWF(testActual).Test ~= 0
        figure();
        bar(diferenciaXiWF(:,testActual));
    end
%%
end




