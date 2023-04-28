
f(x)=A*exp(-(x-mu)**2/(2*sigma**2))
fit f(x) 'r16_acc.dat' using 1:2 via A,mu,sigma

#plot 'r12_acc_N30.dat' u 1:2, f(x)
