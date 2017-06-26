function fco = equalXbmBands(fmin, fmax, N)
% EQUAL_XBM_BANDS - Divide frequency interval into N bands of equal width
% 		    along the human basilar membrane.
% Based on M.C. Liberman's cochlear frequency map for the cat 
% scaled to match human frequency range of hearing.
%	
% Usage: fco = equal_xbm_bands(fmin, fmax, N)
%	fmin	minimum frequency in Hz
%	fmax	maximum frequency in Hz
%	N	number of frequency bands
%	fco	Vector of band cutoff frequencies in Hz (size [1 X (N+1)])
%
%	Copyright Bertrand Delgutte, 1999-2000
%
xmin = Library.cochlearMap(fmin, 20000);
xmax = Library.cochlearMap(fmax, 20000);

dx = (xmax-xmin)/N;
x = xmin:dx:xmax;

fco = Library.invCochlearMap(x, 20000);