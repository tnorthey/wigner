
% load data
r56 = load('../distances/r56.dat', '-ascii')
chi2 = load('chi2_sample.dat', '-ascii')

% initial guess for centre:
sigma0 = mean(r56)
fprintf( 'initial guess for sigma: %f', sigma0 );
mu0 = std(r56)
fprintf( 'initial guess for mu: %f', mu0 );

[sigma, mu] = gaussfit( r56, 1./chi2, sigma0, mu0 )

fprintf( 'final sigma: %f', sigma );
fprintf( 'final mu: %f', mu );
