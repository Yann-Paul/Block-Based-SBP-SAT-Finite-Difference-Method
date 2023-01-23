function Z = gauss(X,Y,sigma,mu)
    X_periodic = min(min(abs(X-mu(1)),abs(X+1-mu(1))),abs(X-1-mu(1)));
    Y_periodic = min(min(abs(Y-mu(2)),abs(Y+1-mu(2))),abs(Y-1-mu(2)));
    Z = exp(-((X_periodic).^2/(2*sigma(1)^2)+(Y_periodic).^2/(2*sigma(2))^2));
end