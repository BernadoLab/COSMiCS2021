cc = jet(length(curvasFibrilacion));

if seleccion == 1;

    figure();
    subplot(2,2,1);
    for i=1:curvasTotales
        plot(valoresS(1:puntosAbs,1), valoresI(1:puntosAbs,i),'Color',cc(i,:));
        hold on;
    end
    title('Absolute scale');

    subplot(2,2,2);
    for i=1:curvasTotales
        plot(valoresS(1:puntosHoltzer,1), valoresI(1:puntosHoltzer,i).*valoresS(1:puntosHoltzer,1),'Color',cc(i,:));
        hold on;
    end
    title('Holtzer plot');

    subplot(2,2,3);
    valoresSsquare = valoresS.*valoresS;
        for i=1:curvasTotales
        plot(valoresS(1:puntosKratky,1), valoresI(1:puntosKratky,i).*valoresSsquare(1:puntosKratky,1),'Color',cc(i,:));
        hold on;
    end
    title('Kratky plot');

    subplot(2,2,4);
    valoresScuad = valoresS.^4;
    for i=1:curvasTotales
        plot(valoresS(1:puntosPorod,1), valoresI(1:puntosPorod,i).*valoresScuad(1:puntosPorod,1),'Color',cc(i,:));
        hold on;
    end
    title('Porod plot');

end

if seleccion == 2
    
    
    
    figure();
    for i=1:length(curvasFibrilacion)
        semilogy(valoresS(:,1),valoresI(:,i),'Color',cc(i,:));
        hold on;
    end
    title('Data set in semilogarithmic scale');


end


if seleccion == 3
    
    
    
    figure();
    for i=1:length(extraData)
        plot(intensidad(:,i),'Color',cc(i,:));
        hold on;
    end
    title('Additional Data Set');


end



