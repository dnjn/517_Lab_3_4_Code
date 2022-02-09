function V_ext = calcVext(currentTrace, currentXYZ, electrodeXYZ)
    sigma = 0.3333e-6;
    V_ext = zeros(758,1);
    for i = 1:758
        for j = 1:531
            r = sqrt((currentXYZ(j,1)-electrodeXYZ(1))^2+ ...
                     (currentXYZ(j,2)-electrodeXYZ(2))^2+ ...
                     (currentXYZ(j,3)-electrodeXYZ(3))^2);
            V_ext(i) = V_ext(i) + currentTrace(j,i)/(4*pi*sigma*r);
        end
    end
end