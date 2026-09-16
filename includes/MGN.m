function [answerX, answerY] = MGN(spectrogram_mag)
    % MGN  Obtain numerical zeros via the Minimal Grid Neighbors (MGN) method.
    %     A grid point is flagged as a zero when |V f|^2 attains a local minimum
    %     over its 8 neighbours; see Section 5.5.1.
    %
    %   Usage:  [answerX, answerY] = MGN(spectrogram_mag)
    %
    %   Input:
    %   spectrogram_mag :   a matrix with the spectrogram magnitude |V f|^2
    %                       (or, more generally, samples of |GAF|) on the grid.
    %
    %   Output:
    %   answerX         :   an array containing the first coordinate of the numerical zeros.
    %   answerY         :   an array containing the second coordinate of the numerical zeros.
    %
    % ---------------------------------------------------------
    rowsNumber          =   size(spectrogram_mag, 1);
    columnsNumber       =   size(spectrogram_mag, 2);

    % We set the shifts to compare a point 'X' with its 1-delta 'A' neighbors:
    % A A A
    % A X A
    % A A A

    distanceDelta = 1;
    shifts = shiftsGrid(distanceDelta);

    % First comparison using 2D matrix indexing
    [rowsInitialSet, columnsInitialSet]     =   find(abs(spectrogram_mag(2:rowsNumber - 1, 2:columnsNumber - 1)) <= abs(spectrogram_mag(2 + shifts(1, 1):rowsNumber - 1 + shifts(1, 1), 2 + shifts(1, 2):columnsNumber - 1 + shifts(1, 2))));
    matrixResult                            =   sub2ind([rowsNumber, columnsNumber], 1 + rowsInitialSet, 1 + columnsInitialSet);

    % Keep only points that stay minimal against each remaining neighbour (linear indexing).
    for j = 2:size(shifts, 1)
        linear_shift = shifts(j, 1) + (shifts(j, 2) * rowsNumber);

        matrixResult = matrixResult(abs(spectrogram_mag(matrixResult)) <= abs(spectrogram_mag(matrixResult + linear_shift)));
    end

    [answerX, answerY]                      =   ind2sub([rowsNumber, columnsNumber], matrixResult);

end

function shifts = shiftsGrid(distanceDelta)
    range = -distanceDelta:distanceDelta;
    [dr, dc] = ndgrid(range, range);
    shifts = [dr(:), dc(:)];
    center_idx = (size(shifts, 1) + 1) / 2;
    shifts(center_idx, :) = [];
end
