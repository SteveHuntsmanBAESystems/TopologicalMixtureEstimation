**README**

**OVERVIEW**
This is research code implementing topological mixture estimation (https://arxiv.org/abs/1712.04487) in MATLAB. The key file is tme1d.m, but this calls other files, including for the simpler technique of topological density estimation.

**CONTENTS**

    ColorIndices.txt    (sample data for an interesting example)
    jsd.m               (Jensen-Shannon divergence [JSD])
    README-TME.txt      (this file)
    tde1d.m             (topological density estimation [TDE])
    tight_subplot.m     (third party file, used for nice-looking plots)
    tme1d.m             (topological mixture estimation [TME])
    unidec.m            (unimodal decomposition)
    unimixmaxjsd.m      (unimodal mixture maximization of the JSD)

**EXAMPLE** 
To run TME on an interesting example, save ColorIndices.txt to a directory whose string representation is D, put the .m files in your MATLAB path, and run the commands
    
        X = getfield(importdata([D,filesep,'ColorIndices.txt'],'',10),'data');
        T = tme1d(X,1); % second argument is a plot flag

This example will take a few seconds to run before displaying figures. The argument '10' in importdata refers to the header lines of the color index file that describe how to get this data from its original source.

**TIPS**
For sample data with n >> 1000, it is generally a good idea to sort and downsample to quantiles, i.e., if X is the sample data and 
    
    n = numel(X), 
    
then take, e.g., 

    X2 = sort(X(:)); X2 = X2(1:(n/1000):end) 

and work with X2. If this gives unsatisfactory results, constant-bandwidth kernel methods are probably not a great idea in the first place, but you might persist and decompose the sample data to allow the bandwidth to vary in a piecewise constant manner. Indeed, applying TME recursively to mixture components is a good idea. A fast Gauss transform would probably yield significant improvements in runtime. This is left as an exercise.

**CITE**
If you use TME in your work, please cite the ICML 2018 proceedings: besides PMLR, the paper is at https://arxiv.org/abs/1712.04487. I would also personally like to hear about your application, evaluation, etc.    

**CONTACT**

    steve.huntsman@baesystems.com

**ACKNOWLEDGEMENTS**
TME code is based upon work supported by the Defense Advanced Research Projects Agency (DARPA) and the Air Force Research Laboratory (AFRL). Any opinions, findings and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of DARPA or AFRL.

**LICENSE INFORMATION** 
TME files (including this one) are distributed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA. TDE files are distributed under a BSD license. The third-party file tight_subplot.m is distributed under its own license. All licenses are included in the files themselves, and this paragraph is for informational purposes only.

Copyright (c) 2018, BAE Systems. All rights reserved.
