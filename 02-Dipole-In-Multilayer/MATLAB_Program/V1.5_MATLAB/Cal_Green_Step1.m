% To calculate the Green tensor's first part
function [PPUpA,PPUpB,SSUpA,SSUpB]=Cal_Green_1(kx,ky,kz)
    SSUpA=zeros(3,3);
    SSUpB=zeros(3,3);

    PPUpA=zeros(3,3);
    PPUpB=zeros(3,3);

    SSUpA(1, 1) = ky^2;
    SSUpA(1, 2) = -kx* ky;
    SSUpA(1, 3) = 0;
    SSUpA(2, 1) = -kx* ky;
    SSUpA(2, 2) = kx^2;
    SSUpA(2, 3) = 0;
    SSUpA(3, 1) = 0;
    SSUpA(3, 2) = 0;
    SSUpA(3, 3) = 0;

    SSUpA = SSUpA./ (kx.^2 + ky.^2)/kz;

    SSUpB(:, :) = SSUpA;


    PPUpA(1, 1) = kx^2 * kz^2;
    PPUpA(1, 2) = kx* ky * kz^2;
    PPUpA(1, 3) = -kx* kz*(kx^2 + ky^2);
    PPUpA(2, 1) = kx*ky*kz^2;
    PPUpA(2, 2) = (ky^2)* kz^2;
    PPUpA(2, 3) = -ky*kz*(kx^2 + ky^2);
    PPUpA(3, 1) = -kx*kz*(kx^2 + ky^2);
    PPUpA(3, 2) = -ky*kz*(kx^2 + ky^2);
    PPUpA(3, 3) = (kx^2 + ky^2)^2;
    PPUpA(:, :) = PPUpA(:, :)./(kx^2 + ky^2 + kz^2) / (kx^2 + ky^2)/kz;

    PPUpB(:,:) = PPUpA;
    PPUpB(1:2, :) = -PPUpB(1:2,:);
end