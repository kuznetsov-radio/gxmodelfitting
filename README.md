This code searches for the coronal heating model parameters that provide the best agreement between the model and observed radio maps.

The two central routines are:

a) MultiScanAB.pro: it searches for the heating rate value Q0 that provides the best agreement between the model and observed radio maps, for the specified parameters a and b of the coronal heating model. An example of calling the code is presented in the file examples\TestCHMP1.pro.

b) SearchForLocalMinimumAB.pro: it searches for the parameters of the coronal heating model (a, b, Q0) that provide the best agreement between the model and observed radio maps. The search provides a local minimum of the selected model-to-observations comparison metric. The program also (optionally) determines the region of "good agreement" in the (a, b) space, where the model-to-observations comparison metric is below a certain threshold (relative to the minimum one). An example of calling the code is presented in the file examples\TestCHMP2.pro.

The data (GX Simulator model, etc.) required to run the examples can be found at https://doi.org/10.5281/zenodo.13999701

Here, we repeat the headers of the above-mentioned routines with the descriptions of the parameters:

pro MultiScanAB, RefDir, ModelFileName, EBTELfileName, LibFileName, OutDir, alist, blist, xc, yc, dx, dy, Nx, Ny [, RefFiles=RefFiles, Q0start=Q0start, threshold=threshold, metric=metric, MultiThermal=MultiThermal, ObsDateTime=ObsDateTime, noMultiFreq=noMultiFreq, DEM=DEM, DDM=DDM, Qstep=Qstep, xy_shift=xy_shift, loud=loud, SHtable=SHtable, Nthreads=Nthreads]

Input parameters:<br/>
RefDir - the directory where the observed radio maps/profiles are stored. If the parameter RefFiles is omitted, the program loads all \*.sav files in the RefDir directory. Otherwise, the program loads the file(s) specified by RefDir+RefFiles.<br/> 
For general 2D radio maps, each .sav file should contain a 'ref' map object with three maps:<br/>
I_obs=ref.getmap(0) - the observed radio map (in terms of brightness temperature in K), with the tags I_obs.freq specifying the emission frequency in GHz, and I_obs.id specifying the map title,<br/>
sigma=ref.getmap(1) - the corresponding instrumental noise (with the same dimensions as I_obs),<br/>
beam =ref.getmap(2) - the instrument beam (point-spread function), with the tags beam.a_beam and beam.b_beam specifying the beam half-widths at 1/e level in two ortogonal directions, in arcseconds.<br/>
Other required tags of these maps are standard for the SSW map structure.<br/>
For 1D intensity profiles observed by RATAN-600, each .sav file should contain two fields:<br/>
instrument='RATAN' - the label to identify the data format,<br/>
ref={flux, x, freq, time, rot} - the structure specifying the data, i.e.,<br/> 
flux - the intensity profile, 1D array, in units of sfu/arcsec,<br/> 
x - the corresponding coordinates, 1D array, in arcsec,<br/>
freq - the emission frequency, in GHz,<br/>
time - the observation time, string (the same as in SSW map structures),<br/>
rot - the RATAN positional angle, in degrees.<br/>
Note that 1D profiles and 2D maps cannot be mixed together.<br/>
ModelFileName - name of the .sav file that contains the GX Simulator model.<br/>
EBTELfileName - name of the .sav file that contains the EBTEL table.<br/>
LibFileName - name of the appropriate executable library that computes the radio emission (see https://github.com/kuznetsov-radio/gximagecomputing).<br/>
OutDir - the directory where the results will be stored.<br/>
alist, blist - a and b parameters of the coronal heating models (scalars or 1D arrays).<br/>
xc, yc - x and y coordinates of the model map center, in arcseconds.<br/>
dx, dy - x and y resolutions (pixel sizes) of the model map, in arcseconds.<br/>
Nx, Ny - x and y sizes of the model map, in pixels.<br/>

The optional parameters include:<br/>
RefFiles - if specified, the program loads the observed radio maps/profiles from the files given by RefDir+RefFiles. Default: all \*.sav files in the RefDir directory.<br/>
Q0start - the initial value of Q0. It can be either: a scalar value (applied to all a an b), or a 2D N_a\*N_b array, where N_a and N_b are the sizes of alist and blist arrays. Default: the best-fit Q0 parameters for AR 12924 (extrapolated to the specified a and b) will be used.<br/>
threshold - the threshold value to compute the image mask. Comparison of the model and observed radio maps is performed in the area where (I_obs gt threshold\*max(I_obs)) || (I_mod gt threshold\*max(I_mod)). Default: 0.1.<br/>
metric - the metric to minimize. It can be one of the following three options:<br/>
'rho': rho^2=mean(((I_obs-I_mod)/I_obs)^2),<br/>
'chi': chi^2=mean(((I_obs-I_mod)/sigma)^2),<br/>
'eta': eta^2=mean(((I_obs-I_mod)/mean(I_obs))^2).<br/>
All averagings are performed over the above-mentioned sub-region determined by the threshold. Default: 'eta'.<br/>
MultiThermal - if set, the formulae for the free-free and gyroresonance emissions from multithermal plasma (described by the DEM and DDM, respectively) are used, see Fleishman, Kuznetsov & Landi (2021). If not set, the isothermal approximation with the plasma density and temperature derived from the DDM or DEM (from the DDM, if both are specified) is used. Default: 0 (isothermal approximation is assumed).<br/>
ObsDateTime - an additional string added to the names of the resulting files. Default: ''<br/>
noMultiFreq - if not set (by default), the code optimizes the computations by computing the metrics at all specified frequencies simultaneously. Although the minimization is performed frequency-by-frequency, the Q0 grid and the corresponding metrics computed during the minimization at lower frequencies are then used as pre-computed data at higher frequencies. If set, all frequencies are processed independently; this can be slower, but sometimes more reliable.<br/>
DEM, DDM - these keywords are only applicable if the chosen EBTELfileName .sav file contains both the DEM and DDM tables. In this case, if the /DEM keyword is set, the code loads the DEM table only (the DDM table is ignored). Similarly, if the /DDM keyword is set, the code loads the DDM table only (the DEM table is ignored). If both /DEM and /DDM keywords (or none of them) are set, the code loads both tables. If the chosen file contains only one EBTEL table (either DEM or DDM), the code loads that table; the /DEM and /DDM keywords are ignored.<br/>
Qstep - the initial relative step over Q0 to search for the optimal heating rate value (must be >1). Default: the golden ratio value (1.6180339).<br/>
xy_shift - shift applied to the observed microwave maps/profiles, a 2-element vector (for 2D maps) in the form of xy_shift=[dx, dy], or a scalar value (for 1D profiles), in arcseconds. If this parameter is not specified (by default), the shift is computed automatically each time (i.e., for each frequency and a, b, and Q0 values) to provide the maximum correlation between the observed and model images/profiles.<br/>
loud - if set, the code displays more detailed information when it fails to find a solution (e.g., when the minimization procedure goes beyond the EBTEL table).<br/>
SHtable - a 7*7 table specifying the selective heating coefficients applied to the field lines with different footpoint combinations. Default: no selective heating (all elements of the table equal 1).<br/>
Nthreads - number of processor threads used for computing the model microwave images. Cannot exceed the number of available processors. Default: a system-defined value (typically, the number of available processors).

Results:<br/>
As the result, for each (a, b) combination the program creates in the OutDir directory a .sav file with the name starting with 'fit' and including the used metric, threshold, indicator of the multithermal approach, a and b values, and (optionally) the ObsDateTime string. These .sav files contain the following fields:<br/>
freqList - array of the emission frequencies, in GHz.<br/>
bestQarr - array of the obtained best-fit heating rates Q0 at different frequencies.<br/>
rhoArr, chiArr, etaArr - arrays of the obtained rho^2, chi^2, and eta^2 metrics at different frequencies. Note that only one of those metrics (defined by the 'metric' keyword) is actually minimized; two other metrics correspond to the obtained best-fit Q0 values.<br/>
modImageConvArr - if the input is in the form of 2D maps, this is a (multi-frequency) map object containing the above-mentioned best-fit model radio maps convolved with the instrument beam.<br/>
obsImageArr - if the input is in the form of 2D maps, this is a (multi-frequency) map object containing the observed radio maps rebinned and shifted to match the best-fit model maps at the corresponding frequencies.<br/>
If the input is in the form of 1D profiles, the fields modImageConvArr and obsImageArr are lists of structures (in the same format as described above) specifying the model 1D scans and the observed 1D scans rebinned and shifted to match the models, respectively.<br/>
If the algorithm failed to find the best-fit heating rate at a certain frequency (e.g., the used metric has no minimum within the valid Q0 range, or has more than one local minimum), the corresponding bestQarr, rhoArr, chiArr, and etaArr are set to NaN, and the corresponding image maps contain all zeros. Note: the program does not overwrite the existing fit*.sav files. If the program is interrupted, on the next launch it will compute the results only for those (a, b) values that have not been processed before.<br/>
Also, the program creates in the OutDir directory a .sav file with the name starting with 'Summary' and including the used metric, threshold, indicator of the multithermal approach, and (optionally) the ObsDateTime string. This .sav file contains the following fields:<br/>
alist, blist - the input alist and blist parameters.<br/>
freqList - array of the emission frequencies, in GHz.<br/>
bestQ - 3D array (N_a\*N_b\*N_freq, where N_a, N_b, and N_freq are the sizes of the alist, blist, and freqList arrays, respectively) of the obtained best-fit heating rates Q0 at different values of a, b, and frequency.<br/>
Iobs, Imod - 3D arrays (N_a\*N_b\*N_freq) of the total observed and model radio fluxes at different values of a, b, and frequency. The fluxes correspond to the obtained best-fit Q0 values.<br/>
CC - 3D array (N_a\*N_b\*N_freq) of the correlation coefficients of the observed and model radio maps at different values of a, b, and frequency. The coefficients correspond to the obtained best-fit Q0 values.<br/>
shiftX, shiftY - 3D arrays (N_a\*N_b\*N_freq) of the shifts (in arcseconds) applied to the observed radio maps to obtain the best correlation with the model maps, at different values of a, b, and frequency. The shifts correspond to the obtained best-fit Q0 values. If the input is in the form of 1D profiles, shiftY is always zero.<br/>
rho, chi, eta - 3D arrays (N_a\*N_b\*N_freq) of the obtained rho^2, chi^2, and eta^2 metrics at different values of a, b, and frequency. Note that only one of those metrics (defined by the 'metric' keyword) is actually minimized; two other metrics correspond to the obtained best-fit Q0 values.<br/>
If the Summary\*.sav file exists, it will be overwritten.

pro SearchForLocalMinimumAB, RefFileName, ModelFileName, EBTELfileName, LibFileName, OutDir, xc, yc, dx, dy, Nx, Ny, a_start, b_start, da, db [, Q0start=Q0start, metric=metric, threshold_img=threshold_img, threshold_metric=threshold_metric, MultiThermal=MultiThermal, ObsDateTime=ObsDateTime, ObsFreq=ObsFreq, DEM=DEM, DDM=DDM, a_range=a_range, b_range=b_range, noArea=noArea, Qstep=Qstep, xy_shift=xy_shift, loud=loud, SHtable=SHtable, Nthreads=Nthreads]

Input parameters:<br/>
RefFileName - name of the .sav file that contains the observed radio map/profile (at a single frequency).<br/> 
The file should contain a 'ref' map object with three maps:<br/>
I_obs=ref.getmap(0) - the observed radio map (in terms of brightness temperature in K), with the tags I_obs.freq specifying the emission frequency in GHz, and I_obs.id specifying the map title,<br/>
sigma=ref.getmap(1) - the corresponding instrumental noise (with the same dimensions as I_obs),<br/>
beam =ref.getmap(2) - the instrument beam (point-spread function), with the tags beam.a_beam and beam.b_beam specifying the beam half-widths at 1/e level in two ortogonal directions, in arcseconds.<br/>
Other required tags of these maps are standard for the SSW map structure.<br/>
For 1D intensity profiles observed by RATAN-600, the .sav file should contain two fields:<br/>
instrument='RATAN' - the label to identify the data format,<br/>
ref={flux, x, freq, time, rot} - the structure specifying the data, i.e.,<br/>
flux - the intensity profile, 1D array, in units of sfu/arcsec,<br/>
x - the corresponding coordinates, 1D array, in arcsec,<br/>
freq - the emission frequency, in GHz,<br/>
time - the observation time, string (the same as in SSW map structures),<br/>
rot - the RATAN positional angle, in degrees.<br/>
ModelFileName - name of the .sav file that contains the GX Simulator model.<br/>
EBTELfileName - name of the .sav file that contains the EBTEL table.<br/>
LibFileName - name of the appropriate executable library that computes the radio emission (see https://github.com/kuznetsov-radio/gximagecomputing).<br/>
OutDir - the directory where the results will be stored.<br/>
xc, yc - x and y coordinates of the model map center, in arcseconds.<br/>
dx, dy - x and y resolutions (pixel sizes) of the model map, in arcseconds.<br/>
Nx, Ny - x and y sizes of the model map, in pixels.<br/>
a_start, b_start - the starting point of the metric minimization algorithm in the (a, b) space.<br/>
da, db - the grid sizes of the metric minimization algorithm in the a and b directions, respectively.<br/>

The optional parameters include:<br/>
Q0start - the starting point of the metric minimization algorithm in the Q0 space. This parameter is applied at the point (a_start, b_start) only; at other (a, b) points, the best-fit value of Q0 at the nearest adjacent point is used. Default: the best-fit Q0 parameters for AR 12924 (extrapolated to the specified a_start and b_start) will be used.<br/>
metric - the metric to minimize. It can be one of the following three options:<br/>
'rho': rho^2=mean(((I_obs-I_mod)/I_obs)^2),<br/>
'chi': chi^2=mean(((I_obs-I_mod)/sigma)^2),<br/>
'eta': eta^2=mean(((I_obs-I_mod)/mean(I_obs))^2).<br/>
All averagings are performed over the sub-region determined by the threshold_img threshold. Default: 'eta'.<br/>
threshold_img - the threshold value to compute the image mask. Comparison of the model and observed radio maps is performed in the area where (I_obs gt threshold_img\*max(I_obs)) || (I_mod gt threshold_img\*max(I_mod)). Default: 0.1.<br/>
threshold_metric - the threshold value that determines the area of "good agreement" between the model and observations in the (a, b) space. The algorithm will scan the (a, b) space with the da and db resolutions, until it finds all points where the condition metric(a, b) < threshold_metric\*min(metric) is satisfied. Default: 2.<br/>
MultiThermal - if set, the formulae for the free-free and gyroresonance emissions from multithermal plasma (described by the DEM and DDM, respectively) are used, see Fleishman, Kuznetsov & Landi (2021). If not set, the isothermal approximation with the plasma density and temperature derived from the DDM or DEM (from the DDM, if both are specified) is used. Default: 0 (isothermal approximation is assumed).<br/>
ObsDateTime - an additional string added to the names of the resulting files. Default: ''.<br/>
ObsFreq - an additional string added to the names of the resulting files. Default: ''.<br/>
DEM, DDM - these keywords are only applicable if the chosen EBTELfileName .sav file contains both the DEM and DDM tables. In this case, if the /DEM keyword is set, the code loads the DEM table only (the DDM table is ignored). Similarly, if the /DDM keyword is set, the code loads the DDM table only (the DEM table is ignored). If both /DEM and /DDM keywords (or none of them) are set, the code loads both tables. If the chosen file contains only one EBTEL table (either DEM or DDM), the code loads that table; the /DEM and /DDM keywords are ignored.<br/>
a_range, b_range - the allowed ranges of the a and b indices, 2-element arrays. The search for the best fit is performed in the range of a_range[0] < a < a_range[1], b_range[0] < b < b_range[1].  Default: a_range=[-10, 10], b_range=[-10, 10]. Note that you cannot extend the a and b ranges beyond the default values, due to the limitations related to the file-naming convention.<br/>
noArea - if set, the area of good agreement (within the threshold_metric threshold) is not computed, i.e., the code stops after finding the local minimum. Default: 0 (the area of good agreement is computed).<br/>
Qstep - the initial relative step over Q0 to search for the optimal heating rate value (must be >1). Default: the golden ratio value (1.6180339).<br/>
xy_shift - shift applied to the observed microwave map/profile, a 2-element vector in the form of xy_shift=[dx, dy] (for 2D maps), or a scalar value (for 1D profiles), in arcseconds.  If this parameter is not specified (by default), the shift is computed automatically each time (i.e., for each frequency and a, b, and Q0 values) to provide the maximum correlation between the observed and model images/profiles. <br/>
loud - if set, the code displays more detailed information when it fails to find a solution (e.g., when the minimization procedure goes beyond the EBTEL table).<br/>
SHtable - a 7*7 table specifying the selective heating coefficients applied to the field lines with different footpoint combinations. Default: no selective heating (all elements of the table equal 1).<br/>
Nthreads - number of processor threads used for computing the model microwave images. Cannot exceed the number of available processors. Default: a system-defined value (typically, the number of available processors).

Results:<br/>
The output of the program is similar to that of the MultiScanAB.pro, with the difference that only one frequency is considered. The program creates in the OutDir directory a .sav file with the name starting with 'Summary' and including the used metric, threshold, indicator of the multithermal approach, and (optionally) the ObsDateTime and ObsFreq strings. This .sav file contains the following fields:<br/>
alist, blist - arrays of a and b parameters covering the range searched by the program during the metric minimization process, with da and db steps.<br/>
freqList - 1-element array containing the emission frequency, in GHz.<br/>
bestQ - 3D array (N_a\*N_b\*1, where N_a and N_b are the sizes of the alist and blist arrays, respectively) of the obtained best-fit heating rates Q0 at different values of a and b.<br/>
ItotalObs, ItotalMod - 3D arrays (N_a\*N_b\*1) of the total observed and model radio fluxes at different values of a and b. The fluxes correspond to the obtained best-fit Q0 values.<br/>
CC - 3D array (N_a\*N_b\*1) of the correlation coefficients of the observed and model radio maps at different values of a and b. The coefficients correspond to the obtained best-fit Q0 values.<br/>
shiftX, shiftY - 3D arrays (N_a\*N_b\*1) of the shifts (in arcseconds) applied to the observed radio maps to obtain the best correlation with the model maps, at different values of a, b, and frequency. The shifts correspond to the obtained best-fit Q0 values. If the input is in the form of 1D profiles, shiftY is always zero.<br/>
rho, chi, eta - 3D arrays (N_a\*N_b\*1) of the obtained rho^2, chi^2, and eta^2 metrics at different values of a and b. Note that only one of those metrics (defined by the 'metric' keyword) is actually minimized; two other metrics correspond to the obtained best-fit Q0 values.<br/>
Note, that if the algorithm has failed to find the best-fit heating rate Q0 at a certain combination of a and b (e.g., the used metric has no minimum within the valid Q0 range, or has more than one local minimum), the corresponding bestQ, ItotalMod, CC, rho, chi, and eta are set to NaN.<br/>
If the data for a certain combination of a and b are missing (because the search for the best-fit Q0 is performed only within a subset of the rectangular area determined by the regular grids alist and blist), the corresponding bestQ, ItotalMod, CC, rho, chi, and eta are set to -1. Therefore, when  analyzing the results, only the points where bestQ is finite and bestQ>0 should be considered. If the Summary\*.sav file exists, it will be overwritten.<br/>
The program saves the temporary progress: for each (a, b) combination it creates in the OutDir directory a .sav file with the name starting with 'fit' and including the used metric, threshold, indicator of the multithermal approach, a and b values, and (optionally) the ObsDateTime and ObsFreq strings. These .sav files contain the following fields:<br/>
freqList - 1-element array containing the emission frequency, in GHz.<br/>
bestQarr - 1-element array containing the obtained best-fit heating rate Q0.<br/>
rhoArr, chiArr, etaArr -1-element arrays containing the obtained rho^2, chi^2, and eta^2 metrics. Note that only one of those metrics (defined by the 'metric' keyword) is actually minimized; two other metrics correspond to the obtained best-fit Q0 values.<br/>
modImageArr - map object containing the best-fit model radio map (corresponding to the best-fit heating rate Q0). The map is not convolved with the instrument beam.<br/>
modImageConvArr - if the input is in form of 2D map, this is a map object containing the above-mentioned best-fit model radio map convolved with the instrument beam.<br/>
obsImageArr - if the input is in form of 2D map, this is a map object containing the observed radio map rebinned and shifted to match the best-fit model map.<br/>
If the input is in the form of 1D profile, the fields modImageConvArr and obsImageArr are (1-element) lists of structures (in the same format as described above) specifying the model 1D scans and the observed 1D scans rebinned and shifted to match the models, respectively.<br/>
If the algorithm failed to find the best-fit heating rate (e.g., the used metric has no minimum within the valid Q0 range, or has more than one local minimum), the corresponding bestQarr, rhoArr, chiArr, and etaArr are set to NaN, and the corresponding image maps contain all zeros.<br/>
Note: the program does not overwrite the existing fit\*.sav files. If the program is interrupted, on the next launch it will compute the results only for those (a, b) values that have not been processed before.
