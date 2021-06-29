function Sum_Diff=Cal_Difference(Z_Mat_1,Z_Mat_2)
    % Normalize the data
    nZ_Mat_1=Z_Mat_1/sqrt(sum(abs(Z_Mat_1).^2,'all'));
    nZ_Mat_2=Z_Mat_2/sqrt(sum(abs(Z_Mat_2).^2,'all'));
    Difference = (nZ_Mat_1 - nZ_Mat_2).^2;
    Sum_Diff = sum(abs(Difference),'all');
end