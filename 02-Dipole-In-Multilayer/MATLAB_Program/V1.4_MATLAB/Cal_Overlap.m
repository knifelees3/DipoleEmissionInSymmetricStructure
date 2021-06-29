function Sum_Overlap=Cal_Overlap(Z_Mat_1,Z_Mat_2)
    % Normalize the data
    nZ_Mat_1=Z_Mat_1/sqrt(sum(abs(Z_Mat_1).^2,'all'));
    nZ_Mat_2=Z_Mat_2/sqrt(sum(abs(Z_Mat_2).^2,'all'));
    Overlap = abs(nZ_Mat_1.*nZ_Mat_2);
    Sum_Overlap = sum(Overlap,'all');
end