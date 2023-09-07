function time_top88(nbglr, nbgtb, nelx, nely, rmin, eta, vmax_c, vmax_t, max_iter, ...
    loadcase, vsm, szmv, Rc, cScale, cover)

% Start timing the entire function
totalTimeStart = cputime;

% Time MATERIAL PROPERTIES section
start = cputime;
top88_TCM(nbglr, nbgtb, nelx, nely, rmin, eta, vmax_c, vmax_t, max_iter, ...
    loadcase, vsm, szmv, Rc, cScale, cover);
elapsedTime_materialProperties = cputime - start;

% ... you'll repeat the above two lines for each of the major sections in your code

% End timing for the entire function
totalTimeElapsed = cputime - totalTimeStart;

% Display the timings
fprintf('\n*** Timing Report ***\n');
fprintf('MATERIAL PROPERTIES: %f seconds\n', elapsedTime_materialProperties);
% ... you'll add similar print statements for each section
fprintf('Total time: %f seconds\n', totalTimeElapsed);

end