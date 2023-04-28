
% load data
r56 = load('../distances/r56.dat', '-ascii')
chi2 = load('chi2_sample.dat', '-ascii')

x = r56;
y = chi2;

[sigmaNew,muNew,Anew]=mygaussfit(x,y);

fprintf( 'final sigma: %f', sigmaNew );
fprintf( 'final mu: %f', muNew );
fprintf( 'final A: %f', Anew );

%y=Anew*exp(-(x-muNew).^2/(2*sigmaNew^2));
%hold on; plot(x,y,'.r');

