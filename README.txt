To convert from .ms to .mzML use proteowizard
    there is also mzXML, but this is a dated file format and only to be used if you are using
    legacy software
    proteowizard was downloaded from github and is written in C++ so you have to use VScode
    to write it to a library so you can use it in the project
    I do not have the C++ knowledge to do this, so instead we will download the end user app
    version of proteowizard and call it as a subprocess of our application
        this will slightly reduce efficiency but due to small data load this will likely not
        be significant

pyOpenMS lecture: https://www.youtube.com/watch?v=XWDFZXYp1qM
    pre-processing (feature detection, smoothing etc...)
    start by using proteowizard to centroid data then use these 3 methods to detect features:
        1. MassTraceDetection
        2. ElutionPeakDetection
        3. FeatureFinderMetabo
            these three methods will result in a feature file that is cleaner than the origonal raw data
        once you have a featurefile you can use
            MapAligner
        to account for run to run retention time variance
        chooses a feature file (with max features, make it a pooled QC sample) as a reference that the rest will be
        fitted to
        FeatureLinking: clusters feature info into a
            ConsensusFeature
        that will have features from different features with similar m/z and RT clustered
        convert feature map to a dataframe that can be easily displayed as a table
        IDMapper: annotates data in files, allows for filtering of non-fragmented features?
        Additional pre-processing suggestions:
            SiriusAdapter: wrapper for SIRIUS algorithm and CSI:FingerID that can generate formula and structure
            predictions
            MetaboliteAdductDecharger in pyopenms: assigns adducts and reconstructs netural masses?
            GNPSExport: generates all nescessary files for Feature-Based Molecular Networking (FBMN) and Ion Identity
            Molecular Networking (IIMN)
                this is explained in another lecture, there are many MS seminars associated with this youtube channel
                it is a very good resource
