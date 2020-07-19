function [RSUp, RPUp, RSDn, RPDn, RS21, RP21, RS12, RP12]=Cal_RSP_kxy(num_dl, num_kx, num_ky, Eplist, dl, klz, nUp, nDn)

 %The reflection coefficients
    RP21 = zeros(num_kx, num_ky, num_dl);
    RP12 = zeros(num_kx, num_ky, num_dl);
    RS21 = zeros(num_kx, num_ky, num_dl);
    RS12 = zeros(num_kx, num_ky, num_dl);

    for l=1:num_dl
        RP21(:, :, l) = (Eplist(l) * klz(:, :, l + 1) - Eplist(l + 1) * klz(:, :, l))./(Eplist(l + 1) * klz(:, :, l) + Eplist(l) * klz(:, :, l + 1));
        RP12(:, :, l) = -RP21(:, :, l);
        RS21(:, :, l) = (klz(:, :, l + 1) - klz(:, :, l))./ (klz(:, :, l) + klz(:, :, l + 1));
        RS12(:, :, l) = -RS21(:, :, l);
    end
    % For The Upper layer
    RPUp = 0;
    RSUp = 0;

    for m=1:nUp
        l = num_dl - m+1;
        exp11 = exp(-1i * (klz(:, :, l+ 1) - klz(:, :, l)) * dl(l));
        exp12 = exp(1i * (klz(:, :, l+1 ) + klz(:, :, l)) * dl(l));
        exp21 = exp(-1i * (klz(:, :, l+1) + klz(:, :, l)) * dl(l));
        exp22 = exp(1i * (klz(:, :, l+1 ) - klz(:, :, l)) * dl(l));

        RPUp = (RPUp.* exp11 + RP12(:, :, l).* exp12)./(RPUp.* RP12(:, :, l).* exp21 + exp22);
        RSUp = (RSUp.* exp11 + RS12(:, :, l).* exp12)./(RSUp.* RS12(:, :, l).* exp21 + exp22);
    end
    %For the lower space
    RPDn = 0;
    RSDn = 0;
    
    for m=1:nDn
        l = m+1;
        exp11 = exp(1i * (klz(:, :, l - 1) - klz(:, :, l)) * dl(l - 1));
        exp12 = exp(-1i * (klz(:, :, l - 1) + klz(:, :, l)) * dl(l - 1));
        exp21 = exp(+1i * (klz(:, :, l - 1) + klz(:, :, l)) * dl(l - 1));
        exp22 = exp(-1i * (klz(:, :, l - 1) - klz(:, :, l)) * dl(l - 1));
        

        RPDn = (RPDn.* exp11 + RP21(:, :, l - 1).* exp12)./(RPDn.* RP21(:, :, l - 1).* exp21 + exp22);
        RSDn = (RSDn.* exp11 + RS21(:, :, l - 1).* exp12)./(RSDn.* RS21(:, :, l - 1).* exp21 + exp22);
    end
end