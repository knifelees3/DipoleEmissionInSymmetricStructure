function Sum_Diff=Cal_Difference(Z_Mat_1,Z_Mat_2)
    % Normalize the data
    nZ_Mat_1=Z_Mat_1/sum(Z_Mat_1,'all');
    nZ_Mat_2=Z_Mat_2/sum(Z_Mat_2,'all');
    Difference = (nZ_Mat_1 - nZ_Mat_2).^2;
    Sum_Diff = sum(sum(abs(Difference)));
end