function [reconstruccion] = reconstCurvas(conc,spectra)

    reconstruccion=zeros(length(spectra),length(conc));   

    for ii = 1:length(conc) % Curvas totales
        for kk = 1:length(spectra) % Puntos totales
            reconstruccion(kk,ii) = 0;
            for jj = 1:length(conc(1,:)) % Numero especies
                reconstruccion(kk,ii) = reconstruccion(kk,ii) + conc(ii,jj).*spectra(jj,kk);
            end
            
        end
    end
    
end