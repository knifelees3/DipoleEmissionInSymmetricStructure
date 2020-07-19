% To transform the data in xy expression into rho,theta phi expression.
% this program will judge the input's krho kphi coordinate
function [Pattern_theta,Pattern_phi]=Transform_RhoPhi(num_kx,num_l,krho_grid0,kphi_grid,krho0,kphi,PatternXY)

Pattern_theta=zeros(num_l,1);
Pattern_phi=zeros(num_l,1);

for l=1:num_l
    for m=1:num_kx
        for n=1:num_kx
            if l~=num_l
                if (krho0(l+1)>abs(krho_grid0(m,n))&& abs(krho_grid0(m,n))>=krho0(l))
                        Pattern_theta(l+1)=Pattern_theta(l+1)+PatternXY(m,n);
                end
                if (kphi(l+1)>abs(kphi_grid(m,n))&& abs(kphi_grid(m,n))>=kphi(l))
                    Pattern_phi(l)=Pattern_phi(l)+PatternXY(m,n);
                end
            else
               if (krho0(l)==abs(krho_grid0(m,n)))
                     Pattern_theta(l)=Pattern_theta(l)+PatternXY(m,n);
               end
                if (kphi(l)==abs(kphi_grid(m,n)))
                        Pattern_phi(l)=Pattern_phi(l)+PatternXY(m,n);
                end
            end
        end
    end
end

%% Plot part as a test. Should be commented in using
% figure(1)
% plot(krho,Pattern_theta,'s')
% title('Along \rho')
% xlabel('k_{\rho}/k_{0}')
% ylabel('Intensity')
% %%
% figure(2)
% plot(kphi,Pattern_phi,'s')
% title('Along \phi')
% xlabel('\phi')
% ylabel('Intensity')
end