Niter = 100;
Nmax = 20;
cap_av = zeros(Nmax,1);
cap = zeros(Nmax,Niter);
for N = 1:Nmax
    for iter = 1:Niter
        chMax = (random('Rician',0,2,N,1)+1i*random('Rician',0,2,N,1));
%         chMax = (raylrnd(1,[N 1])+1i*raylrnd(1,[N 1]));
        cap(N,iter) = log2(1 + sum(abs(chMax)));
    end
    cap_av(N) = mean(cap(N,:));
end

figure;
plot(1:1:N,cap_av,'lineWidth',2);
grid minor;