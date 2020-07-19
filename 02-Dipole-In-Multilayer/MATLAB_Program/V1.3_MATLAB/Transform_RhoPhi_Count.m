% Transform the data from Cartesian coordinate into rho phi coordinate
% the data in different range will be count
% Actually this code is useless
function [Pattern_rho,Pattern_phi]=Transform_RhoPhi_Count(num_kx,num_l,krho_grid0,kphi_grid,krho0,kphi,PatternXY)

Pattern_rho=zeros(num_l,1);
Pattern_phi=zeros(num_l,1);

Count_Rho=zeros(num_l,1);
Count_Phi=zeros(num_l,1);

for l=1:num_l
    for m=1:num_kx
        for n=1:num_kx
            if l~=num_l
                if (krho0(l+1)>abs(krho_grid0(m,n))&& abs(krho_grid0(m,n))>=krho0(l))
                        Pattern_rho(l+1)=Pattern_rho(l+1)+PatternXY(m,n);
                        Count_Rho(l)=Count_Rho(l)+1;
                end
                if (kphi(l+1)>abs(kphi_grid(m,n))&& abs(kphi_grid(m,n))>=kphi(l))
                    Pattern_phi(l)=Pattern_phi(l)+PatternXY(m,n);
                    Count_Phi(l)=Count_Phi(l)+1;
                end
            else
               if (krho0(l)==abs(krho_grid0(m,n)))
                     Pattern_rho(l)=Pattern_rho(l)+PatternXY(m,n);                
                     Count_Rho(l)=Count_Rho(l)+1;
               end
                if (kphi(l)==abs(kphi_grid(m,n)))
                      Pattern_phi(l)=Pattern_phi(l)+PatternXY(m,n);
                      Count_Phi(l)=Count_Phi(l)+1;
                end
            end
        end
    end
end

Pattern_rho=Pattern_rho.*Count_Rho;
Pattern_phi=Pattern_phi.*Count_Phi;
%
% figure(1)
% plot(krho,Pattern_rho,'s')
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