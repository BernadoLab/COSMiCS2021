%% Create folder for the report

colorscheme = jet(length(valoresIActual(1,:)));

createFolder = [folderOutput,'/Report/'];

mkdir(createFolder);



%fid = fopen(fileReport, 'w+' );


%disp(' ');
%str =   'Name of the experiment \n ';
%experimentName = input(str,'s');
experimentName = (OutputFolder);
if isempty(experimentName)
    experimentName = 'Experiment 1';
end
%%

set(gcf,'Visible','off');
fig = figure('visible', 'off');
for i=1:length(curvasFibrilacion)
    fig = semilogy(valoresS(1:puntosAbs,1),matrizInput(i,1:puntosAbs)','Color',colorscheme(i,:));
    hold on;
end
title('Matriz input (Absolute values) in semilogarithmic scale');

fileFig = [folderOutput, '/Report/datasetSemilog.png'];
saveas(fig, fileFig);



fig = figure('visible', 'off');
for i=1:length(curvasFibrilacion)
    fig = plot(valoresS(1:puntosHoltzer,1),matrizHoltzer(i,1:puntosHoltzer)','Color',colorscheme(i,:));
    hold on;
end
title('Matriz input (Holtzer values)');

fileFig = [folderOutput, '/Report/datasetHoltzer.png'];
saveas(fig, fileFig);

fig = figure('visible', 'off');;
for i=1:length(curvasFibrilacion)
    fig = plot(valoresS(1:puntosKratky,1),matrizKratky(i,1:puntosKratky)','Color',colorscheme(i,:));
    hold on;
end
title('Matriz input (Kratky values)');

fileFig = [folderOutput, '/Report/datasetKratky.png'];
saveas(fig, fileFig);

fig = figure('visible', 'off');
for i=1:length(curvasFibrilacion)
    fig = plot(valoresS(1:puntosPorod,1),matrizPorod(i,1:puntosPorod)','Color',colorscheme(i,:));
    hold on;
end
title('Matriz input (Porod values)');

fileFig = [folderOutput, '/Report/datasetPorod.png'];
saveas(fig, fileFig);


fig = figure('visible', 'off');
    
    fig = bar(eigenvalues(1:10,1));
    title('Eigenvalues');
    text(9,eigenvalues(1,1), num2str(eigenvalues(1:10,:)),'HorizontalAlignment','left','VerticalAlignment','baseline');
    
    fileFig = [folderOutput, '/Report/eigenvalues.png'];
    saveas(fig, fileFig);
    
    
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

    fileFig = [folderOutput, '/Report/eigenvectors.png'];
    fig = gcf;
    saveas(fig, fileFig);
    

close all;




%%

createFile = [createFolder,'style.css'];

copyfile('/CMPT/people/sagar/softwares/Cosmics-multi/style.css', createFile);

fileReport = [createFolder,'report.html'];


fid = fopen(fileReport, 'w+' );
  
fprintf( fid, '<!DOCTYPE html>');fprintf(fid, '\n');
fprintf( fid, '<html><head>');fprintf(fid, '\n');


fprintf( fid, '<link href=''http://fonts.googleapis.com/css?family=Roboto:400,300,100,500'' rel=''stylesheet'' type=''text/css''>');fprintf(fid, '\n');
fprintf( fid,	'<link href=''http://fonts.googleapis.com/css?family=Roboto+Slab:400,300,100,500'' rel=''stylesheet'' type=''text/css''>');fprintf(fid, '\n');
fprintf( fid,'<link rel="stylesheet" href="style.css">');fprintf(fid, '\n');

fprintf( fid,	'</head>');fprintf(fid, '\n');
fprintf( fid,	'<body>');fprintf(fid, '\n');
		
fprintf( fid,	'<header id="fh5co-header" role="banner">');fprintf(fid, '\n');
fprintf( fid,	'	<nav class="navbar navbar-default" role="navigation">');fprintf(fid, '\n');
fprintf( fid,	'		<div class="container">');fprintf(fid, '\n');
fprintf( fid,	'			<div class="row">');fprintf(fid, '\n');
fprintf( fid,	'				<div class="col-md-10 col-md-offset-1">');fprintf(fid, '\n');
fprintf( fid,	'					<div class="navbar-header"> ');fprintf(fid, '\n');

fprintf( fid,	'					<a href="file:///Users/shelkie/Desktop/index2.html#" class="js-fh5co-nav-toggle fh5co-nav-toggle visible-xs-block" data-toggle="collapse" data-target="#fh5co-navbar" aria-expanded="false" aria-controls="navbar"><i></i></a>');fprintf(fid, '\n');
fprintf( fid,	'					<a class="navbar-brand" href="https://www.ncbi.nlm.nih.gov/pubmed/27889205" target="_blank">COSMiCS</a>');fprintf(fid, '\n');
fprintf( fid,	'					</div>');fprintf(fid, '\n');
fprintf( fid,	'					<div id="fh5co-navbar" class="navbar-collapse collapse">');fprintf(fid, '\n');
fprintf( fid,	'						<ul class="nav navbar-nav navbar-right">');fprintf(fid, '\n');
fprintf( fid,	'							<li class="active"><a href="#dataInfo"><span>Data info <span class="border"></span></span></a></li>');fprintf(fid, '\n');
fprintf( fid,	'							<li><a href="#analysisOptions"><span>Analysis options <span class="border"></span></span></a></li>');fprintf(fid, '\n');
fprintf( fid,	'							<li><a href="#PCA"><span>PCA <span class="border"></span></span></a></li>');fprintf(fid, '\n');
fprintf( fid,	'							<li><a href="#Results"><span>Results <span class="border"></span></span></a></li>');fprintf(fid, '\n');
fprintf( fid,	'							</ul>');fprintf(fid, '\n');
fprintf( fid,	'						</div>');fprintf(fid, '\n');
fprintf( fid,	'					</div>');fprintf(fid, '\n');
fprintf( fid,	'				</div>');fprintf(fid, '\n');
fprintf( fid,	'			</div>');fprintf(fid, '\n');
fprintf( fid,	'		</nav>');fprintf(fid, '\n');
fprintf( fid,	'	</header>');fprintf(fid, '\n');
	
fprintf( fid,	'<div id="fh5co-main">');fprintf(fid, '\n');

fprintf( fid,	'		<div class="fh5co-intro text-center">');fprintf(fid, '\n');
fprintf( fid,	'			<div class="container">');fprintf(fid, '\n');
fprintf( fid,	'				<div class="row">');fprintf(fid, '\n');
fprintf( fid,	'					<div class="col-md-8 col-md-offset-2">');fprintf(fid, '\n');
						
fprintf( fid,	'						<h1 class="intro-lead"> %s </h1>', experimentName );fprintf(fid, '\n');

fprintf( fid,	'					</div>');fprintf(fid, '\n');
fprintf( fid,	'				</div>');fprintf(fid, '\n');
fprintf( fid,	'			</div>');fprintf(fid, '\n');
fprintf( fid,	'		</div>');fprintf(fid, '\n');

% Seccion: Data information
fprintf( fid,	'			<div class="container">');fprintf(fid, '\n');
fprintf( fid,	'				<div class="row">');fprintf(fid, '\n');
fprintf( fid,	'					<div class="col-md-10 col-md-offset-1">');fprintf(fid, '\n');
fprintf( fid,	'						<div class="row">');fprintf(fid, '\n');
fprintf( fid,	'							<div class="col-md-3 animate-box">');fprintf(fid, '\n');
fprintf( fid,	'								<a name="dataInfo"></a>');fprintf(fid, '\n');
fprintf( fid,	'								<h3>Data information</h3>');fprintf(fid, '\n');
%fprintf( fid,	'								<p>Data information about the dataset, folder of the input data and name of the files</p>');fprintf(fid, '\n');

fprintf( fid,	'							</div>');fprintf(fid, '\n');
fprintf( fid,	'							<div class="col-md-9">');fprintf(fid, '\n');
% Aqui va la informacion de la seccion Data information
fprintf( fid,	'								<p class="animate-box"><b>Input folder: </b> %s </p>', OutputFolder);fprintf(fid, '\n');
fprintf( fid,	'								<p class="animate-box"><b>Output folder: </b> %s </p>', folderOutput);fprintf(fid, '\n');
	fprintf( fid,	'<br>');
    fprintf( fid,	'								<p class="animate-box"><b>LIST OF FILES: </b> </p>');
    for i=1:length(OutFiles)
        fprintf(fid, 'Curve %s - %s',num2str(i), OutFiles(i).name);
        fprintf(fid, '<br>');
    end 
 fprintf( fid,	'								<br><p class="animate-box"><b>Total files: </b> %d </p>', curvasTotales);fprintf(fid, '\n');   
 
 
 fprintf( fid,	'<img class="small" src="datasetSemilog.png" alt="Dataset semilogaritmic scale">');fprintf(fid, '\n');
 fprintf( fid,	'<img class="small" src="datasetHoltzer.png" alt="Dataset Holtzer scale">');fprintf(fid, '\n');
 fprintf( fid,	'<img class="small" src="datasetKratky.png" alt="Dataset Kratky scale">');fprintf(fid, '\n');
 fprintf( fid,	'<img class="small" src="datasetPorod.png" alt="Dataset Porod scale">');fprintf(fid, '\n');
 
 
fprintf( fid,	'						</div>');fprintf(fid, '\n');    
fprintf( fid,	'					</div>');fprintf(fid, '\n'); 
fprintf( fid,	'				</div>');fprintf(fid, '\n');    
 

fprintf( fid,	'<br>');fprintf(fid, '\n');    
fprintf( fid,	'<br>');fprintf(fid, '\n');    


colorscheme = jet(length(curvasFibrilacion));

fprintf( fid,	'<br>');fprintf(fid, '\n');    
fprintf( fid,	'<br>');fprintf(fid, '\n');    



 % Seccion: Analysis options
fprintf( fid,	'			<div class="container">');fprintf(fid, '\n');
fprintf( fid,	'				<div class="row">');fprintf(fid, '\n');
fprintf( fid,	'					<div class="col-md-10 col-md-offset-1">');fprintf(fid, '\n');
fprintf( fid,	'						<div class="row">');fprintf(fid, '\n');
fprintf( fid,	'							<div class="col-md-3 animate-box">');fprintf(fid, '\n');
fprintf( fid,	'								<a name="analysisOptions"></a>');fprintf(fid, '\n');
fprintf( fid,	'								<h3>Analysis options</h3>');fprintf(fid, '\n');
%fprintf( fid,	'								<p>Data information about the dataset, folder of the input data and name of the files</p>');fprintf(fid, '\n');

fprintf( fid,	'							</div>');fprintf(fid, '\n');
fprintf( fid,	'							<div class="col-md-9">');fprintf(fid, '\n');
% Aqui va la informacion de la seccion Analysis options
fprintf( fid,	'								<p class="animate-box"><b>Points removed at the beginning: </b> %d </p>', elimPoints);fprintf(fid, '\n');
fprintf( fid,	'								<p class="animate-box"><b>Momentum transfer range: </b></p>');fprintf(fid, '\n');

if unidades == 'A' || unidades == 'a'
    fprintf(fid, '<p class="animate-box">- Absolute: %4.2f &Aring;<sup>-1</sup></p>', corteAbs);
    fprintf(fid, '<p class="animate-box">- Holtzer: %4.2f &Aring;<sup>-1</sup></p>', corteHoltzer);
    fprintf(fid, '<p class="animate-box">- Kratky: %4.2f &Aring;<sup>-1</sup></p>', corteKratky);
    fprintf(fid, '<p class="animate-box">- Porod: %4.2f &Aring;<sup>-1</sup></p>', cortePorod);
end


if unidades == 'N' || unidades == 'n'
    fprintf(fid, '<p class="animate-box">- Absolute: %4.2f 1/nm\n</p>', corteAbs);
    fprintf(fid, '<p class="animate-box">- Holtzer: %4.2f 1/nm\n</p>', corteHoltzer);
    fprintf(fid, '<p class="animate-box">- Kratky: %4.2f 1/nm\n</p>', corteKratky)
    fprintf(fid, '<p class="animate-box">- Porod: %4.2f 1/nm\n</p>', cortePorod);
    
end

fprintf(fid, '<p class="animate-box"><b>Number of species:</b> %d </p>', numeroEspecies);
 
  %fprintf(fid, '<br>');
  
  fprintf( fid,	'			<br><p class="animate-box"><b>CONSTRAINTS </b> %d </p>');fprintf(fid, '\n');   
  
  if tipoALS ==1
    
      fprintf(fid,'<p class="animate-box"><b>Constraints for fibrillation experiment</b></p>');fprintf(fid, '\n'); 
fprintf(fid, '<ul>');fprintf(fid, '\n'); 
    fprintf(fid, '<li>Non-negativity: Concentrations and spectra</li>');fprintf(fid, '\n'); 
    fprintf(fid, '<ul>');fprintf(fid, '\n'); 
    fprintf(fid, '<li>Selected algorithm for concentrations: fnnls (all species)</li>');fprintf(fid, '\n'); 
    fprintf(fid, '<li>Selected algorithm for spectra: fnnls (all species)</li>');fprintf(fid, '\n'); 
    fprintf(fid, '</ul>');fprintf(fid, '\n'); 
    fprintf(fid, '<li>Closure: Concentrations</li>');fprintf(fid, '\n'); 
    fprintf(fid, '<ul>');fprintf(fid, '\n'); 
    fprintf(fid, '<li>Number of closure constants: 1</li>');fprintf(fid, '\n'); 
    fprintf(fid, '<li>Closure constant: 1.0</li>');fprintf(fid, '\n');     
    fprintf(fid, '<li>Closure condition: Equal condition (all species)</li>');fprintf(fid, '\n'); 
   fprintf(fid, '</ul>');fprintf(fid, '\n'); 
    fprintf(fid, '</ul>');fprintf(fid, '\n'); 
end

if(tipoALS == 2)
    
    fprintf(fid,'<p class="animate-box"><b>Constraints for titration experiment</b></p>');fprintf(fid, '\n'); 
fprintf(fid, '<ul>');fprintf(fid, '\n'); 
fprintf(fid, '<li>Non-negativity: Concentrations and spectra</li>');fprintf(fid, '\n'); 
    fprintf(fid, '<ul>');fprintf(fid, '\n'); 
    fprintf(fid, '<li>Selected algorithm for concentrations: fnnls (all species)</li>');fprintf(fid, '\n'); 
    fprintf(fid, '<li>Selected algorithm for spectra: fnnls (all species)</li>');fprintf(fid, '\n'); 
    fprintf(fid, '</ul>');fprintf(fid, '\n'); 
    fprintf(fid, '<li>Unimodality: Concentrations (averaged)</li>');fprintf(fid, '\n'); 
    
    fprintf(fid, '</ul>');fprintf(fid, '\n'); 


end

fprintf( fid,	'<br>');fprintf(fid, '\n');    
fprintf( fid,	'<br>');fprintf(fid, '\n');    


fprintf(fid,'<p class="animate-box"><b>Convergence criterion:</b> %.1f </p>', tolsigma);

fprintf(fid,'<p class="animate-box"><b>Maximum number of iterations:</b> %.0f </p>', numIterations);





fprintf( fid,	'							</div>');fprintf(fid, '\n'); 
fprintf( fid,	'						</div>');fprintf(fid, '\n');    
fprintf( fid,	'					</div>');fprintf(fid, '\n'); 
fprintf( fid,	'				</div>');fprintf(fid, '\n');   



fprintf( fid,	'<br>');fprintf(fid, '\n');    
fprintf( fid,	'<br>');fprintf(fid, '\n');    

% Seccion: PCA
fprintf( fid,	'			<div class="container">');fprintf(fid, '\n');
fprintf( fid,	'				<div class="row">');fprintf(fid, '\n');
fprintf( fid,	'					<div class="col-md-10 col-md-offset-1">');fprintf(fid, '\n');
fprintf( fid,	'						<div class="row">');fprintf(fid, '\n');
fprintf( fid,	'							<div class="col-md-3 animate-box">');fprintf(fid, '\n');
fprintf( fid,	'								<a name="PCA"></a>');fprintf(fid, '\n');
fprintf( fid,	'								<h3>Principal Component Analysis</h3>');fprintf(fid, '\n');
%fprintf( fid,	'								<p>Data information about the dataset, folder of the input data and name of the files</p>');fprintf(fid, '\n');

% fprintf( fid,	'							</div>');fprintf(fid, '\n');
% fprintf( fid,	'							<div class="col-md-9">');fprintf(fid, '\n');

fprintf( fid,	'<div class="fh5co-portfolio-item ">');fprintf(fid, '\n');
fprintf( fid,	'				<div class="fh5co-portfolio-figure animate-box" style="background-image: url(images/work_1.jpg);"></div>');fprintf(fid, '\n');
fprintf( fid,	'				<div class="fh5co-portfolio-description">');fprintf(fid, '\n');



fprintf( fid,	'<img class="small" src="eigenvalues.png" alt="Eigenvalues" >');fprintf(fid, '\n');
 fprintf( fid,	'<img class="small" src="eigenvectors.png" alt="Eigenvectors">');fprintf(fid, '\n');

 
  

fprintf( fid,	'						</div>');fprintf(fid, '\n');    
fprintf( fid,	'					</div>');fprintf(fid, '\n'); 
fprintf( fid,	'				</div>');fprintf(fid, '\n'); 

fprintf( fid,	'						</div>');fprintf(fid, '\n');    
fprintf( fid,	'					</div>');fprintf(fid, '\n'); 
fprintf( fid,	'				</div>');fprintf(fid, '\n'); 



% Seccion: Results
fprintf( fid,	'			<div class="container">');fprintf(fid, '\n');
fprintf( fid,	'				<div class="row">');fprintf(fid, '\n');
fprintf( fid,	'					<div class="col-md-10 col-md-offset-1">');fprintf(fid, '\n');
fprintf( fid,	'						<div class="row">');fprintf(fid, '\n');
fprintf( fid,	'							<div class="col-md-3 animate-box">');fprintf(fid, '\n');




fprintf( fid,	'								<a name="Results"></a>');fprintf(fid, '\n');
fprintf( fid,	'								<h3>Analysis results</h3>');fprintf(fid, '\n');

% fprintf( fid,	'							</div>');fprintf(fid, '\n');
% fprintf( fid,	'							<div class="col-md-9">');fprintf(fid, '\n');

fprintf( fid,	'<div class="fh5co-portfolio-item ">');fprintf(fid, '\n');
fprintf( fid,	'				<div class="fh5co-portfolio-figure animate-box" style="background-image: url(images/work_1.jpg);"></div>');fprintf(fid, '\n');
fprintf( fid,	'				<div class="fh5co-portfolio-description">');fprintf(fid, '\n');


%fprintf(fid, '<pre>');fprintf(fid, '\n');
%fprintf(fid, 'Test \t Combination \t Rem.curves \t Init.estim. \t l.o.f (PCA) \t l.o.f (exp) \t Variance \t X^2 \t Convergence');fprintf(fid, '\n');
%fprintf(fid, '-------------------------------------------------------------------------------------------------------------------------------------');

fprintf(fid, '<table>');fprintf(fid, '\n');
fprintf(fid, '  <tr>');fprintf(fid, '\n');
fprintf(fid, '    <th>Test</th>');fprintf(fid, '\n');
fprintf(fid, '    <th>Combination</th> ');fprintf(fid, '\n');
fprintf(fid, '    <th>Rem.curves</th>');fprintf(fid, '\n');
fprintf(fid, '    <th>Init.estim.</th>');fprintf(fid, '\n');
fprintf(fid, '    <th>l.o.f (PCA)</th>');fprintf(fid, '\n');
fprintf(fid, '    <th>l.o.f (exp)</th>');fprintf(fid, '\n');
fprintf(fid, '    <th>Variance</th>');fprintf(fid, '\n');
fprintf(fid, '    <th>chi^2</th>');fprintf(fid, '\n');
fprintf(fid, '    <th>Convergence</th>');fprintf(fid, '\n');
fprintf(fid, '  </tr>');fprintf(fid, '\n');


    
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


    
    if (statisticsWorkflow(il).solution ==1)
        solucion = 'Convergence';
    end
    if (statisticsWorkflow(il).solution ==2)
        solucion = 'Divergence';
    end
    if (statisticsWorkflow(il).solution ==3)
        solucion = 'Num. iteration exceeded';
    end
    
    textEstim = [];
    
    for in = 1:length(curvasEstimInicWF(il).Test)
        textEstim = [ textEstim, num2str(curvasEstimInicWF(il).Test(in))];
        if (in<length(curvasEstimInicWF(il).Test) && length(curvasEstimInicWF(il).Test) > 1)
            textEstim = [ textEstim, ' , '];  
        end
    end
    
    textElim = [];
    
    for in = 1:length(curvasEliminadasWF(il).Test)
        textElim = [textElim, num2str(curvasEliminadasWF(il).Test(in))];
        if (in<length(curvasEliminadasWF(il).Test) && length(curvasEliminadasWF(il).Test) > 1)
            textElim = [ textElim, ' , '];  
        end
    end
    
 %   fprintf('%d \t %s \t %s \t\t %s \t %f \t %f \t %f \t %.2f \t %s \n', il, comb2, textElim, textEstim, statisticsWorkflow(il).lackOfFit_PCA, statisticsWorkflow(il).lackOfFit_Exp, statisticsWorkflow(il).percentR2, chiAverageWorkflow(il,1), solucion);
    
 


 
fprintf(fid, '  <tr>');fprintf(fid, '\n');
fprintf(fid, '    <td><b>%d</b></td>', il);fprintf(fid, '\n');
fprintf(fid, '    <td>%s</td> ', comb2);fprintf(fid, '\n');
fprintf(fid, '    <td>%s</td>', textElim);fprintf(fid, '\n');
fprintf(fid, '    <td>%s</td>', textEstim);fprintf(fid, '\n');
fprintf(fid, '    <td>%f</td>', statisticsWorkflow(il).lackOfFit_PCA);fprintf(fid, '\n');
fprintf(fid, '    <td>%f</td>', statisticsWorkflow(il).lackOfFit_Exp);fprintf(fid, '\n');
fprintf(fid, '    <td>%f</td>', statisticsWorkflow(il).percentR2);fprintf(fid, '\n');
fprintf(fid, '    <td>%.2f</td>', chiAverageWorkflow(il,1));fprintf(fid, '\n');
fprintf(fid, '    <td>%s</td>',solucion );fprintf(fid, '\n');
fprintf(fid, '  </tr>');fprintf(fid, '\n');


 
 
    if(mod(il, length(combinacionesWF(:,1))) == 0)
        fprintf(fid, '<br>');fprintf(fid, '\n');
    end
    
    
    clear comb;
end
     
fprintf(fid, '</table>');fprintf(fid, '\n');





for indice=1:length(matricesUsadasWF(:,1))
    
    fprintf(fid, '<br>');fprintf(fid, '\n');
    
    fprintf(fid, '<h4>Test %d</h4>', indice);fprintf(fid, '\n');
    
    
    
    comb = ['Test ' num2str(indice), ': Absolute'];
    if (matricesUsadasWF(indice,2)==1),  comb = [comb, ' + Holtzer']; end; 
    if (matricesUsadasWF(indice,3)==1),  comb = [comb, ' + Kratky']; end;
    if (matricesUsadasWF(indice,4)==1),  comb = [comb, ' + Porod']; end;
    
    figure(indice);
    subplot(2,1,1);
    semilogy(valoresS(1:puntosAbs),speciesWorkflow(indice).Test(1:puntosAbs,:));
    texto = ['xi2 = ', num2str(chiAverageWorkflow(indice,1))];
    text(0.8, 0.8, texto, 'Units', 'normalized');
    title(comb);
    if statisticsWorkflow(indice).solution ~= 1
        text(0.8, 0.7, 'DIVERGENCE!', 'Units', 'normalized','Color','Red');
    end
    subplot(2,1,2);
    plot(coptElimWorkflow(indice).Test);
    title('Concentrations');
    
    fileFig = [folderOutput, '/Report/solution', num2str(indice), '.png'];
    
    fig = gcf;
    saveas(fig, fileFig);
    
    close(fig);
    
    fprintf( fid,	'<img class="grande" src="solution%d.png" alt="Solution%d" >', indice, indice);fprintf(fid, '\n');

    figure(indice);
    plot(chiSquareAll(:,indice), ':*');
    


    fileFig = [folderOutput, '/Report/chiTest', num2str(indice), '.png'];
    fig = gcf;
    saveas(fig, fileFig);
    close(fig);
    
    fprintf( fid,	'<img class="grande" src="chiTest%d.png" alt="Eigenvalues" >', indice);fprintf(fid, '\n');
    
    
    if curvasEliminadasWF(indice).Test ~= 0
        figure(indice);
        bar(diferenciaXiWF(:,indice));
        
        fileFig = [folderOutput, '/Report/diffChi', num2str(indice), '.png'];
        fig = gcf;
        saveas(fig, fileFig);
        close(fig);
        fprintf( fid,	'<img class="grande" src="diffChi%d.png" alt="Eigenvalues" >', indice);fprintf(fid, '\n');
        
    end
    
end
  %fprintf(fid, '</pre>');fprintf(fid, '\n');

fprintf( fid,	'						</div>');fprintf(fid, '\n');    
fprintf( fid,	'					</div>');fprintf(fid, '\n'); 
fprintf( fid,	'				</div>');fprintf(fid, '\n'); 

fprintf( fid,	'						</div>');fprintf(fid, '\n');    
fprintf( fid,	'					</div>');fprintf(fid, '\n'); 
fprintf( fid,	'				</div>');fprintf(fid, '\n'); 



  
  
  

   
    
   
    
    

fprintf( fid,	'							<p><a href="report.html" class="btn btn-primary">Back to top</a></p>');fprintf(fid, '\n');
fprintf( fid,	'						</div>');fprintf(fid, '\n');
fprintf( fid,	'					</div>');fprintf(fid, '\n');
fprintf( fid,	'				</div>');fprintf(fid, '\n');
fprintf( fid,	'			</div>');fprintf(fid, '\n');
fprintf( fid,	'		</div>');fprintf(fid, '\n');







%fprintf( fid,	'</div>');

fprintf( fid,	'	<footer id="fh5co-footer">');fprintf(fid, '\n');    
fprintf( fid,	'		<div class="container">');fprintf(fid, '\n');    
fprintf( fid,	'			<div class="row">');fprintf(fid, '\n');    
fprintf( fid,	'				<div class="col-md-10 col-md-offset-1 text-center">');fprintf(fid, '\n');    
fprintf( fid,	'					<p>COSMiCS <br>Created by <a href="mailto:fatimaht@gmail.com" target="_blank">Fatima Herranz-Trillo</a></p>');fprintf(fid, '\n');    
fprintf( fid,	'				</div>');fprintf(fid, '\n');    
fprintf( fid,	'			</div>');fprintf(fid, '\n');    
fprintf( fid,	'		</div>');fprintf(fid, '\n');    
fprintf( fid,	'	</footer>');fprintf(fid, '\n');    

fprintf( fid,	'</body></html>');


fclose(fid);
