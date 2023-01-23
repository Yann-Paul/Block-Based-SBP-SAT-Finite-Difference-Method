function Z = sin_cos(X,Y,a)
%     Z = zeros(size(X));
%     for j=1:100
%         Z = Z + sin(k*2*pi*X-k)/100 + cos(k*2*pi*Y-k)/100;
%     end
    
    n = max(size(a));
    Z = zeros(size(X));
    for j = 1:n
        %Z= Z + a(j)/n*(sin(2*pi*j*X-j.^2)+cos(2*pi*j*Y-j));
        Z= Z + a(j)/n*(sin(2*pi*j*X)+cos(2*pi*j*Y-j.^2));
    end
end