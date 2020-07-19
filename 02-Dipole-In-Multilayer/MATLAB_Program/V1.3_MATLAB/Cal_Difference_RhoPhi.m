function [Diff_rho,Diff_phi]=Cal_Difference_RhoPhi(P_rho,P_phi,Exp_rho,Exp_phi)

            % Normalize the value
            n_Exp_rho=Exp_rho/sum(Exp_rho,'all');
            n_Exp_phi=Exp_phi/sum(Exp_phi,'all');
            
            nP_rho=P_rho/sum(P_rho,'all');
            nP_phi=P_phi/sum(P_phi,'all');
            % Calculate the difference
            Diff_rho=sum((nP_rho-n_Exp_rho).^2,'all');
            Diff_phi=sum((nP_phi-n_Exp_phi).^2,'all');
end
