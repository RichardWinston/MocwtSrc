In the traditional method of characteristics for groundwater solute-transport models, advective transport is represented by moving particles that track concentration. This approach can lead to global mass-balance problems because in models of aquifers having complex boundary conditions and heterogeneous properties, particles can originate in cells having different pore volumes and (or) be introduced (or removed) at cells representing fluid sources (or sinks) of varying strengths. Use of volume-weighted particles means that each particle tracks solute mass. In source/sink cells, the changes in particle weights will match the volume of water added or removed through external fluxes. This enables the new method to conserve mass in source/sink cells as well as globally. This approach also leads to potential efficiencies by allowing the number of particles per cell to vary spatially—using more particles where concentration gradients are high and fewer where gradients are low. The approach also eliminates the need for the model user to have to distinguish between “weak” and “strong” fluid source (or sink) cells. The new model determines whether solute mass added by fluid sources in a cell should be represented by (1) new particles having weights representing appropriate fractions of the volume of water added by the source, or (2) distributing the solute mass added over all particles already in the source cell. The first option is more appropriate for the condition of a strong source. The latter option is more appropriate for a weak source. At sinks, decisions whether or not to remove a particle are replaced by a reduction in particle weight in proportion to the volume of water removed. A number of test cases demonstrate that the new method works well and conserves mass. The method is incorporated into a new version of the U.S. Geological Survey's MODFLOW-GWT solute-transport model. 

All the models in this archive are hypothetical models.

The models in this archive are used to help document the new version of MODFLOW-GWT using volume-weighted particles.

The digital object identifier for this archive is doi:10.5066/F78050RV

References:
Hsieh, P.A., and Winston, R.B., 2002, User’s Guide To Model Viewer, A Program For Three-Dimensional Visualization of Ground-water Model Results: U.S. Geological Survey Open-File Report 02-106, 18 p.

Konikow, L.F., 1977, Modeling chloride movement in the alluvial aquifer at the Rocky Mountain Arsenal, Colorado: U.S. Geological Survey Water-Supply Paper 2044, 43 p.

Konikow, L.F., Goode, D.J., and Hornberger, G.Z., 1996, A Three- Dimensional Method-of-Characteristics Solute-Transport Model (MOC3D): U.S. Geological Survey Water-Resources Investigations Report 96-4267, 87 p.

Wexler, E.J., 1992, Analytical solutions for one-, two-, and three-dimensional solute transport in ground-water systems with uniform flow: U.S. Geological Survey Techniques of Water-Resources Investigations, book 3, chap. B7, 190 p.

Winston, R.B., Konikow, L.F., and Hornberger, G.Z. 2017. Volume-Weighted Particle-Tracking Method for Solute-Transport Modeling: Implementation in MODFLOW-GWT. U.S. Geological Survey, Techniques and Methods 6-A##.
    Files:
    -----

        model.Finite.1D
            -----------
            Finite.1D is a one-dimensional MODFLOW-GWT model. Flow in the model
            is steady-state. Transport is transient. It has three observation
            locations. The purpose of the model is to test the MOCWT algorithm
            by comparing the simulated concentrations at the observation
            locations with an analytical solution. The analytical solution is
            included in this archive as Finite1. Another similar MODFLOW-GWT
            model in this archive is FiniteAlpha1 which has a higher
            longitudinal dispersivity.
            This folder contains input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

        model.rmamocel
            -----------
            Rmamocel is a MODFLOW-GWT model using the ELLAM method. It is one of
            several models used to compare results between MODFLOW-GWT and
            MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

        model.rmamocim
            -----------
            Rmamocim is a MODFLOW-GWT model using the MOCIMP method. It is one
            of several models used to compare results between MODFLOW-GWT and
            MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

        model.rmamocwt_dist
            -----------
            Rmamocwt_dist is a MODFLOW-GWT model using the MOCWT method. It is
            one of several models used to compare results between MODFLOW-GWT
            and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

        model.rma_moc
            -----------
            Rma_moc is a MODFLOW-GWT model using the MOC method. It is one of
            several models used to compare results between MODFLOW-GWT and
            MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

        model.fd
            -----------
            Fd is a MODFLOW-2000 model. It is used in conjunction with fd_MT3D.
            It is one of several models used to compare results between
            MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mf2k.exe from the bin directory into the this
            directory and double-click on it.

        model.fd_MT3D
            -----------
            Fd_MT3D is an MT3DMS model using the Finite-Difference algorithm. It
            is used in conjunction with fd. It is one of several models used to
            compare results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mt3dms5b.exe from the bin directory into the
            this directory and the flow-transport link file (*.ftl) from the
            corresponding MODFLOW-2000 output directory and double-click on
            mt3dms5b.exe. At the prompt, enter the name of the MT3DMS name file
            (*.mnm).

        model.hmoc
            -----------
            Hmoc is a MODFLOW-2000 model. It is used in conjunction with
            hmoc_MT3D. It is one of several models used to compare results
            between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mf2k.exe from the bin directory into the this
            directory and double-click on it.

        model.hmoc_MT3D
            -----------
            Hmoc_MT3D is an MT3DMS model using the Hybrid-MOC algorithm. It is
            used in conjunction with hmoc. It is one of several models used to
            compare results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mt3dms5b.exe from the bin directory into the
            this directory and the flow-transport link file (*.ftl) from the
            corresponding MODFLOW-2000 output directory and double-click on
            mt3dms5b.exe. At the prompt, enter the name of the MT3DMS name file
            (*.mnm).

        model.mmoc
            -----------
            Mmoc is a MODFLOW-2000 model. It is used in conjunction with
            mmoc_MT3D. It is one of several models used to compare results
            between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mf2k.exe from the bin directory into the this
            directory and double-click on it.

        model.mmoc_MT3D
            -----------
            Mmoc_MT3D is an MT3DMS model using the Modified MOC algorithm. It is
            used in conjunction with mmoc. It is one of several models used to
            compare results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mt3dms5b.exe from the bin directory into the
            this directory and the flow-transport link file (*.ftl) from the
            corresponding MODFLOW-2000 output directory and double-click on
            mt3dms5b.exe. At the prompt, enter the name of the MT3DMS name file
            (*.mnm).

        model.moc
            -----------
            Moc is a MODFLOW-2000 model. It is used in conjunction with
            moc_MT3D. It is one of several models used to compare results
            between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mf2k.exe from the bin directory into the this
            directory and double-click on it.

        model.moc_MT3D
            -----------
            Moc_MT3D is an MT3DMS model using the MOC algorithm. It is used in
            conjunction with moc. It is one of several models used to compare
            results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mt3dms5b.exe from the bin directory into the
            this directory and the flow-transport link file (*.ftl) from the
            corresponding MODFLOW-2000 output directory and double-click on
            mt3dms5b.exe. At the prompt, enter the name of the MT3DMS name file
            (*.mnm).

        model.tvd
            -----------
            Tvd is a MODFLOW-2000 model. It is used in conjunction with
            tvd_MT3D. It is one of several models used to compare results
            between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mf2k.exe from the bin directory into the this
            directory and double-click on it.

        model.tvd_MT3D
            -----------
            Tvd_MT3D is an MT3DMS model using the total-variation-diminishing
            algorithm. It is used in conjunction with tvd. It is one of several
            models used to compare results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mt3dms5b.exe from the bin directory into the
            this directory and the flow-transport link file (*.ftl) from the
            corresponding MODFLOW-2000 output directory and double-click on
            mt3dms5b.exe. At the prompt, enter the name of the MT3DMS name file
            (*.mnm).

        model.rmamocwti
            -----------
            Rmamocwti is a MODFLOW-GWT model using the MOCWTI method. It is one
            of several models used to compare results between MODFLOW-GWT and
            MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

        model.DiracUnidirectional
            -----------
            DiracUnidirectional is a MODFLOW-GWT model using the MOCWTI method.
            Its purpose is to compare the numerical solution with an analytical
            solution to the same problem.
            This folder contains input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

        model.DiracAngle
            -----------
            DiracAngle is a MODFLOW-GWT model using the MOCWTI method. Its
            purpose is to compare the numerical solution with an analytical
            solution to the same problem. It differs from DiracUnidirectional in
            that flow is at a 45 degree angle to the grid instead of parallel to
            the rows.
            This folder contains input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

        model.Point3
            -----------
            Point3 is a MODFLOW-GWT model of a continuous source in a uniform
            three-dimensional flow field. It is used for testing MODFLOW-GWT by
            comparing its results with the Point3Analytical model in this
            archive. This folder contains the input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

        model.Finite1
            -----------
            This model is an analytical model documented in Wexler (1992) and
            available from http://water.usgs.gov/software/ANALGWST/. The program
            models a continuous solute source in a one-dimensional uniform flow
            field. The particular program in Wexler (1992) used for this model
            is named Finite. This folder contains the input files for the model.
            The output of this file is used for comparison with the Finite.1D
            MODFLOW-GWT model in this archive. Another similar analytical Finite
            model in this archive is Finite2 which has a higher longidudinal
            dispersivity.
            To run the model, copy Finite.exe from the bin directory into this
            directory and double-click on it. Then specify the Finite1.dat as
            the input file.

        model.Finite2
            -----------
            This model is an analytical model documented in Wexler (1992) and
            available from http://water.usgs.gov/software/ANALGWST/. The program
            models a continuous solute source in a one-dimensional uniform flow
            field. The particular program in Wexler (1992) used for this model
            is named Finite. This folder contains the input files for the model.
            The output of this file is used for comparison with the Finite.1D
            MODFLOW-GWT model in this archive. Another similar analytical Finite
            model in this archive is Finite1 which has a lower longidudinal
            dispersivity.
            To run the model, copy Finite.exe from the bin directory into this
            directory and double-click on it. Then specify the Finite2.dat as
            the input file.

        model.Point3Analytical
            -----------
            This model is an analytical model documented in Wexler (1992) and
            available from http://water.usgs.gov/software/ANALGWST/. The program
            models a continuous solute source in a three-dimensional uniform
            flow field. The particular program in Wexler (1992) used for this
            model is named Point3. This folder contains the input files for the
            model. The output of this file is used for comparison with the
            Point3 MODFLOW-GWT model in this archive.
            To run the model, copy Point3.exe from the bin directory into this
            directory and double-click on it. Then specify the moc9.z025.dat as
            the input file.

        model.Dirac
            -----------
            Dirac contains the input for an analytical solution to the problem
            of an instantaneous solute source in a uniform 3D flow field. The
            code was originally described in Konikow and others (1996). The
            analytical solution is from Wexler (1992). The purpose of the model
            is for comparison with the DiracUnidirectional and DiracAngle
            MODFLOW-GWT models.
            To run the model, copy Dirac.exe into this folder and double-click
            it. select uni.3d-120D.a-L1z2.5.dat as the input file.

        model.FiniteAlpha1
            -----------
            FiniteAlpha1 is a one-dimensional MODFLOW-GWT model. Flow in the
            model is steady-state. Transport is transient. It has three
            observation locations. The purpose of the model is to test the MOCWT
            algorithm by comparing the simulated concentrations at the
            observation locations with an analytical solution. The analytical
            solution is included in this archive as Finite2. Another similar
            MODFLOW-GWT model in this archive is Finite.1D which has a lower
            longitudinal dispersivity.
            This folder contains input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

        output.Finite.1D
            -----------
            Finite.1D is a one-dimensional MODFLOW-GWT model. Flow in the model
            is steady-state. Transport is transient. It has three observation
            locations. The purpose of the model is to test the MOCWT algorithm
            by comparing the simulated concentrations at the observation
            locations with an analytical solution.
            This folder contains output files for the model.

        output.rmamocel
            -----------
            Rmamocel is a MODFLOW-GWT model using the ELLAM method. It is one of
            several models used to compare results between MODFLOW-GWT and
            MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

        output.rmamocim
            -----------
            Rmamocim is a MODFLOW-GWT model using the MOCIMP method. It is one
            of several models used to compare results between MODFLOW-GWT and
            MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

        output.rmamocwti
            -----------
            Rmamocwti is a MODFLOW-GWT model using the MOCWTI method. It is one
            of several models used to compare results between MODFLOW-GWT and
            MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

        output.rmamocwt_dist
            -----------
            Rmamocwt_dist is a MODFLOW-GWT model using the MOCWT method. It is
            one of several models used to compare results between MODFLOW-GWT
            and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

        output.rma_moc
            -----------
            Rma_moc is a MODFLOW-GWT model using the MOC method. It is one of
            several models used to compare results between MODFLOW-GWT and
            MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

        output.fd
            -----------
            Fd is a MODFLOW-2000 model. It is used in conjunction with fd_MT3D.
            It is one of several models used to compare results between
            MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

        output.fd_MT3D
            -----------
            Fd_MT3D is an MT3DMS model using the Finite-Difference algorithm. It
            is used in conjunction with fd. It is one of several models used to
            compare results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

        output.hmoc
            -----------
            Hmoc is a MODFLOW-2000 model. It is used in conjunction with
            hmoc_MT3D. It is one of several models used to compare results
            between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

        output.hmoc_MT3D
            -----------
            Hmoc_MT3D is an MT3DMS model using the Hybrid-MOC algorithm. It is
            used in conjunction with hmoc. It is one of several models used to
            compare results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

        output.mmoc
            -----------
            Mmoc is a MODFLOW-2000 model. It is used in conjunction with
            mmoc_MT3D. It is one of several models used to compare results
            between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.

        output.mmoc_MT3D
            -----------
            Mmoc_MT3D is an MT3DMS model using the Modified MOC algorithm. It is
            used in conjunction with mmoc. It is one of several models used to
            compare results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

        output.moc
            -----------
            Moc is a MODFLOW-2000 model. It is used in conjunction with
            moc_MT3D. It is one of several models used to compare results
            between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

        output.moc_MT3D
            -----------
            Moc_MT3D is an MT3DMS model using the MOC algorithm. It is used in
            conjunction with moc. It is one of several models used to compare
            results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

        output.tvd
            -----------
            Tvd is a MODFLOW-2000 model. It is used in conjunction with
            tvd_MT3D. It is one of several models used to compare results
            between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

        output.tvd_MT3D
            -----------
            Tvd_MT3D is an MT3DMS model using the total-variation-diminishing
            algorithm. It is used in conjunction with tvd. It is one of several
            models used to compare results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

        output.DiracUnidirectional
            -----------
            DiracUnidirectional is a MODFLOW-GWT model using the MOCWTI method.
            Its purpose is to compare the numerical solution with an analytical
            solution to the same problem.
            This folder contains output files for the model.

        output.DiracAngle
            -----------
            DiracAngle is a MODFLOW-GWT model using the MOCWTI method. Its
            purpose is to compare the numerical solution with an analytical
            solution to the same problem. It differs from DiracUnidirectional in
            that flow is at a 45 degree angle to the grid instead of parallel to
            the rows.
            This folder contains output files for the model.

        output.Point3
            -----------
            Point3 is a MODFLOW-GWT model of a continuous source in a uniform
            three-dimensional flow field. It is used for testing MODFLOW-GWT by
            comparing its results with the Point3Analytical model in this
            archive. This folder contains the out
            put files for the model.

        output.Finite1
            -----------
            This model is an analytical model documented in Wexler (1992) and
            available from http://water.usgs.gov/software/ANALGWST/. The program
            models a continuous solute source in a one-dimensional uniform flow
            field. The particular program in Wexler (1992) used for this model
            is named Finite. This folder contains the output file for the model.
            The output of this file is used for comparison with the Finite.1D
            MODFLOW-GWT model in this archive. Another similar analytical Finite
            model in this archive is Finite2 which has a higher longidudinal
            dispersivity.

        output.Finite2
            -----------
            This model is an analytical model documented in Wexler (1992) and
            available from http://water.usgs.gov/software/ANALGWST/. The program
            models a continuous solute source in a one-dimensional uniform flow
            field. The particular program in Wexler (1992) used for this model
            is named Finite. This folder contains the output file for the model.
            The output of this file is used for comparison with the Finite.1D
            MODFLOW-GWT model in this archive. Another similar analytical Finite
            model in this archive is Finite1 which has a lower longidudinal
            dispersivity.

        output.Point3Analytical
            -----------
            This model is an analytical model documented in Wexler (1992) and
            available from http://water.usgs.gov/software/ANALGWST/.  The
            particular program in Wexler (1992) used for this model is named
            Point3. This folder contains the output files for the model. The
            output of this file is used for comparison with the Point3
            MODFLOW-GWT model in this archive.

        output.Dirac
            -----------
            Dirac contains the output for an analytical solution to the problem
            of an instantaneous solute source in a uniform 3D flow field. The
            code was originally described in Konikow and others (1996). The
            analytical solution is from Wexler (1992). The purpose of the model
            is for comparison with the DiracUnidirectional and DiracAngle
            MODFLOW-GWT models.

        output.FiniteAlpha1
            -----------
            FiniteAlpha1 is a one-dimensional MODFLOW-GWT model. Flow in the
            model is steady-state. Transport is transient. It has three
            observation locations. The purpose of the model is to test the MOCWT
            algorithm by comparing the simulated concentrations at the
            observation locations with an analytical solution. The analytical
            solution is included in this archive as Finite2. Another similar
            MODFLOW-GWT model in this archive is Finite.1D which has a lower
            longitudinal dispersivity.
            This folder contains output files for the model.

        Point3
            -----------
            Point3 contains files related to the Point3 MODFLOW-GWT model.
    Files:
    -----
        readme.txt:
            Readme.txt is a human readable file documenting the archive.

        modelgeoref.txt:
            modelgeoref.txt is a file that documents the location of the model
            in digital degrees. Because all the models in this archive are
            hypothetical models, they do not have a true location. For this
            archive, modelgeoref.txt gives the coordinates of the water science
            center.


    \bin\
        Description:
        -----------
        This directory contains the executables used to do the analysis

        Files:
        -----
        bin\MOCWT.exe:
            MODFLOW-GWT with the addition of the MOCWT and MOCWTI solute
            transport algorithms.

        bin\Dirac.exe:
            Dirac is a compiled version of a code for an analytical solution to
            the problem of an instantaneous solute source in a uniform 3D flow
            field. The code was originally described in Konikow and others
            (1996). The analytical solution is from Wexler (1992). Dirac.exe is
            compiled for the Windows operating system.

        bin\Finite.exe:
            Finite.exe is a compiled version of the Finite program (Wexler,
            1992) for the Windows operating system.

        bin\Point3.exe:
            Point3.exe is a compiled version of the Point3 program (Wexler,
            1992) for the Windows operating system.  The program gives an
            analytical solution to a continuous solute source in a
            three-dimensional uniform flow field.

        bin\mt3dms5b.exe:
            MT3DMS version 5.3

        bin\mf2k.exe:
            MODFLOW-2000 version 1.19


    \georef\
        Description:
        -----------
        This directory contains a shape file defining the active and inactive
        areas of the model used for benchmarking MODFLOW-GWT for comparison with
        MT3DMS.
        The test problem represents an analog based on, and greatly simplified
        from, the groundwater contamination problem at the Rocky Mountain
        Arsenal, Colorado (see Konikow, 1977).

        Files:
        -----
        georef\MOCWT.dbf:
            Shapefile attribute database

        georef\MOCWT.prj:
            Projection file

        georef\MOCWT.shp:
            Shapefile shapes file

        georef\MOCWT.shx:
            Shapefile shape index file


        model\model.Finite.1D\
            Description:
            -----------
            Finite.1D is a one-dimensional MODFLOW-GWT model. Flow in the model
            is steady-state. Transport is transient. It has three observation
            locations. The purpose of the model is to test the MOCWT algorithm
            by comparing the simulated concentrations at the observation
            locations with an analytical solution. The analytical solution is
            included in this archive as Finite1. Another similar MODFLOW-GWT
            model in this archive is FiniteAlpha1 which has a higher
            longitudinal dispersivity.
            This folder contains input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

            Files:
            -----
            model\model.Finite.1D\finite.bas:
                MODFLOW Basic Package input file

            model\model.Finite.1D\finite.bcf:
                MODFLOW Block-Centered Flow Package input file

            model\model.Finite.1D\finite.dis:
                MODFLOW Discretization file

            model\model.Finite.1D\finite.gwt:
                MODFLOW-GWT Transport name file

            model\model.Finite.1D\finite.mob:
                MODFLOW-GWT Obsesrvation input file

            model\model.Finite.1D\finite.moc:
                MODFLOW-GWT Transport input file

            model\model.Finite.1D\finite.nam:
                MODFLOW Name file

            model\model.Finite.1D\finite.oc:
                MODFLOW Output Control input file

            model\model.Finite.1D\finite.sip:
                MODFLOW Strongly Implicit Procedure Package input file

            model\model.Finite.1D\MF2K_GWT.BAT:
                Batch file used to run a model

            model\model.Finite.1D\modflow.bf:
                MODFLOW-2000 batch input file


        model\model.rmamocel\
            Description:
            -----------
            Rmamocel is a MODFLOW-GWT model using the ELLAM method. It is one of
            several models used to compare results between MODFLOW-GWT and
            MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

            Files:
            -----
            model\model.rmamocel\MF2K_GWT.BAT:
                Batch file used to run a model

            model\model.rmamocel\modflow.bf:
                MODFLOW-2000 batch input file

            model\model.rmamocel\rmamocel.bas:
                MODFLOW Basic Package input file

            model\model.rmamocel\rmamocel.dis:
                MODFLOW Discretization file

            model\model.rmamocel\rmamocel.gwt:
                MODFLOW-GWT Transport name file

            model\model.rmamocel\rmamocel.lpf:
                MODFLOW Layer Property Flow Package input file

            model\model.rmamocel\rmamocel.mob:
                MODFLOW-GWT Obsesrvation input file

            model\model.rmamocel\rmamocel.moc:
                MODFLOW-GWT Transport input file

            model\model.rmamocel\rmamocel.nam:
                MODFLOW Name file

            model\model.rmamocel\rmamocel.oc:
                MODFLOW Output Control input file

            model\model.rmamocel\rmamocel.pcg:
                MODFLOW Preconditioned Conjugate Gradient Package input file

            model\model.rmamocel\rmamocel.wel:
                MODFLOW Well Package input file


        model\model.rmamocim\
            Description:
            -----------
            Rmamocim is a MODFLOW-GWT model using the MOCIMP method. It is one
            of several models used to compare results between MODFLOW-GWT and
            MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

            Files:
            -----
            model\model.rmamocim\MF2K_GWT.BAT:
                Batch file used to run a model

            model\model.rmamocim\modflow.bf:
                MODFLOW-2000 batch input file

            model\model.rmamocim\rmamocim.bas:
                MODFLOW Basic Package input file

            model\model.rmamocim\rmamocim.dis:
                MODFLOW Discretization file

            model\model.rmamocim\rmamocim.gwt:
                MODFLOW-GWT Transport name file

            model\model.rmamocim\rmamocim.lpf:
                MODFLOW Layer Property Flow Package input file

            model\model.rmamocim\rmamocim.mob:
                MODFLOW-GWT Obsesrvation input file

            model\model.rmamocim\rmamocim.moc:
                MODFLOW-GWT Transport input file

            model\model.rmamocim\rmamocim.nam:
                MODFLOW Name file

            model\model.rmamocim\rmamocim.oc:
                MODFLOW Output Control input file

            model\model.rmamocim\rmamocim.pcg:
                MODFLOW Preconditioned Conjugate Gradient Package input file

            model\model.rmamocim\rmamocim.wel:
                MODFLOW Well Package input file


        model\model.rmamocwt_dist\
            Description:
            -----------
            Rmamocwt_dist is a MODFLOW-GWT model using the MOCWT method. It is
            one of several models used to compare results between MODFLOW-GWT
            and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

            Files:
            -----
            model\model.rmamocwt_dist\MF2K_GWT.BAT:
                Batch file used to run a model

            model\model.rmamocwt_dist\modflow.bf:
                MODFLOW-2000 batch input file

            model\model.rmamocwt_dist\rmamocwt.bas:
                MODFLOW Basic Package input file

            model\model.rmamocwt_dist\rmamocwt.dis:
                MODFLOW Discretization file

            model\model.rmamocwt_dist\rmamocwt.gwt:
                MODFLOW-GWT Transport name file

            model\model.rmamocwt_dist\rmamocwt.ipda:
                MODFLOW-GWT Initial Particle Density File-Array-Based input file

            model\model.rmamocwt_dist\rmamocwt.lpf:
                MODFLOW Layer Property Flow Package input file

            model\model.rmamocwt_dist\rmamocwt.mob:
                MODFLOW-GWT Obsesrvation input file

            model\model.rmamocwt_dist\rmamocwt.moc:
                MODFLOW-GWT Transport input file

            model\model.rmamocwt_dist\rmamocwt.nam:
                MODFLOW Name file

            model\model.rmamocwt_dist\rmamocwt.oc:
                MODFLOW Output Control input file

            model\model.rmamocwt_dist\rmamocwt.pcg:
                MODFLOW Preconditioned Conjugate Gradient Package input file

            model\model.rmamocwt_dist\rmamocwt.prtp:
                MODFLOW-GWT PTRP input file

            model\model.rmamocwt_dist\rmamocwt.wel:
                MODFLOW Well Package input file

            model\model.rmamocwt_dist\rmamocwt_dist.bas:
                MODFLOW Basic Package input file

            model\model.rmamocwt_dist\rmamocwt_dist.dis:
                MODFLOW Discretization file

            model\model.rmamocwt_dist\rmamocwt_dist.gwt:
                MODFLOW-GWT Transport name file

            model\model.rmamocwt_dist\rmamocwt_dist.ipda:
                MODFLOW-GWT Initial Particle Density File-Array-Based input file

            model\model.rmamocwt_dist\rmamocwt_dist.lpf:
                MODFLOW Layer Property Flow Package input file

            model\model.rmamocwt_dist\rmamocwt_dist.mob:
                MODFLOW-GWT Obsesrvation input file

            model\model.rmamocwt_dist\rmamocwt_dist.moc:
                MODFLOW-GWT Transport input file

            model\model.rmamocwt_dist\rmamocwt_dist.nam:
                MODFLOW Name file

            model\model.rmamocwt_dist\rmamocwt_dist.oc:
                MODFLOW Output Control input file

            model\model.rmamocwt_dist\rmamocwt_dist.pcg:
                MODFLOW Preconditioned Conjugate Gradient Package input file

            model\model.rmamocwt_dist\rmamocwt_dist.prtp:
                MODFLOW-GWT PTRP input file

            model\model.rmamocwt_dist\rmamocwt_dist.wel:
                MODFLOW Well Package input file


        model\model.rma_moc\
            Description:
            -----------
            Rma_moc is a MODFLOW-GWT model using the MOC method. It is one of
            several models used to compare results between MODFLOW-GWT and
            MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

            Files:
            -----
            model\model.rma_moc\MF2K_GWT.BAT:
                Batch file used to run a model

            model\model.rma_moc\modflow.bf:
                MODFLOW-2000 batch input file

            model\model.rma_moc\rmamocex.bas:
                MODFLOW Basic Package input file

            model\model.rma_moc\rmamocex.dis:
                MODFLOW Discretization file

            model\model.rma_moc\rmamocex.gwt:
                MODFLOW-GWT Transport name file

            model\model.rma_moc\rmamocex.lpf:
                MODFLOW Layer Property Flow Package input file

            model\model.rma_moc\rmamocex.mob:
                MODFLOW-GWT Obsesrvation input file

            model\model.rma_moc\rmamocex.moc:
                MODFLOW-GWT Transport input file

            model\model.rma_moc\rmamocex.nam:
                MODFLOW Name file

            model\model.rma_moc\rmamocex.oc:
                MODFLOW Output Control input file

            model\model.rma_moc\rmamocex.pcg:
                MODFLOW Preconditioned Conjugate Gradient Package input file

            model\model.rma_moc\rmamocex.wel:
                MODFLOW Well Package input file

            model\model.rma_moc\rmamocwt.ipda:
                MODFLOW-GWT Initial Particle Density File-Array-Based input file


        model\model.fd\
            Description:
            -----------
            Fd is a MODFLOW-2000 model. It is used in conjunction with fd_MT3D.
            It is one of several models used to compare results between
            MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mf2k.exe from the bin directory into the this
            directory and double-click on it.

            Files:
            -----
            model\model.fd\MODFLOW.BAT:
                Batch file used to run a model

            model\model.fd\modflow.bf:
                MODFLOW-2000 batch input file

            model\model.fd\MT3DMS.BAT:
                Batch file used to run a model

            model\model.fd\rmamtfd.bas:
                MODFLOW Basic Package input file

            model\model.fd\rmamtfd.dis:
                MODFLOW Discretization file

            model\model.fd\rmamtfd.lmt:
                MODFLOW Link to MT3DMS or MT3D-USGS input file

            model\model.fd\rmamtfd.lpf:
                MODFLOW Layer Property Flow Package input file

            model\model.fd\rmamtfd.nam:
                MODFLOW Name file

            model\model.fd\rmamtfd.oc:
                MODFLOW Output Control input file

            model\model.fd\rmamtfd.pcg:
                MODFLOW Preconditioned Conjugate Gradient Package input file

            model\model.fd\rmamtfd.wel:
                MODFLOW Well Package input file


        model\model.fd_MT3D\
            Description:
            -----------
            Fd_MT3D is an MT3DMS model using the Finite-Difference algorithm. It
            is used in conjunction with fd. It is one of several models used to
            compare results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mt3dms5b.exe from the bin directory into the
            this directory and the flow-transport link file (*.ftl) from the
            corresponding MODFLOW-2000 output directory and double-click on
            mt3dms5b.exe. At the prompt, enter the name of the MT3DMS name file
            (*.mnm).

            Files:
            -----
            model\model.fd_MT3D\rmamtfd.adv:
                MT3DMS or MT3D-USGS Advection Package input file

            model\model.fd_MT3D\rmamtfd.btn:
                MT3DMS or MT3D-USGS Basic Transport Package input file

            model\model.fd_MT3D\rmamtfd.dsp:
                MT3DMS or MT3D-USGS Dispersion Package input file

            model\model.fd_MT3D\rmamtfd.gcg:
                MT3DMS or MT3D-USGS Generalized Conjugate Gradient Solver
                Package input file

            model\model.fd_MT3D\rmamtfd.mnm:
                MT3DMS or MT3D-USGS Name file

            model\model.fd_MT3D\rmamtfd.ssm:
                MT3DMS or MT3D-USGS Sink and Source Mixing Package input file


        model\model.hmoc\
            Description:
            -----------
            Hmoc is a MODFLOW-2000 model. It is used in conjunction with
            hmoc_MT3D. It is one of several models used to compare results
            between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mf2k.exe from the bin directory into the this
            directory and double-click on it.

            Files:
            -----
            model\model.hmoc\MODFLOW.BAT:
                Batch file used to run a model

            model\model.hmoc\modflow.bf:
                MODFLOW-2000 batch input file

            model\model.hmoc\MT3DMS.BAT:
                Batch file used to run a model

            model\model.hmoc\rmamthm.bas:
                MODFLOW Basic Package input file

            model\model.hmoc\rmamthm.dis:
                MODFLOW Discretization file

            model\model.hmoc\rmamthm.lmt:
                MODFLOW Link to MT3DMS or MT3D-USGS input file

            model\model.hmoc\rmamthm.lpf:
                MODFLOW Layer Property Flow Package input file

            model\model.hmoc\rmamthm.nam:
                MODFLOW Name file

            model\model.hmoc\rmamthm.oc:
                MODFLOW Output Control input file

            model\model.hmoc\rmamthm.pcg:
                MODFLOW Preconditioned Conjugate Gradient Package input file

            model\model.hmoc\rmamthm.wel:
                MODFLOW Well Package input file

            model\model.hmoc\rmamthmoc.bas:
                MODFLOW Basic Package input file

            model\model.hmoc\rmamthmoc.dis:
                MODFLOW Discretization file

            model\model.hmoc\rmamthmoc.lmt:
                MODFLOW Link to MT3DMS or MT3D-USGS input file

            model\model.hmoc\rmamthmoc.lpf:
                MODFLOW Layer Property Flow Package input file

            model\model.hmoc\rmamthmoc.nam:
                MODFLOW Name file

            model\model.hmoc\rmamthmoc.oc:
                MODFLOW Output Control input file

            model\model.hmoc\rmamthmoc.pcg:
                MODFLOW Preconditioned Conjugate Gradient Package input file

            model\model.hmoc\rmamthmoc.wel:
                MODFLOW Well Package input file


        model\model.hmoc_MT3D\
            Description:
            -----------
            Hmoc_MT3D is an MT3DMS model using the Hybrid-MOC algorithm. It is
            used in conjunction with hmoc. It is one of several models used to
            compare results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mt3dms5b.exe from the bin directory into the
            this directory and the flow-transport link file (*.ftl) from the
            corresponding MODFLOW-2000 output directory and double-click on
            mt3dms5b.exe. At the prompt, enter the name of the MT3DMS name file
            (*.mnm).

            Files:
            -----
            model\model.hmoc_MT3D\rmamthm.adv:
                MT3DMS or MT3D-USGS Advection Package input file

            model\model.hmoc_MT3D\rmamthm.btn:
                MT3DMS or MT3D-USGS Basic Transport Package input file

            model\model.hmoc_MT3D\rmamthm.dsp:
                MT3DMS or MT3D-USGS Dispersion Package input file

            model\model.hmoc_MT3D\rmamthm.gcg:
                MT3DMS or MT3D-USGS Generalized Conjugate Gradient Solver
                Package input file

            model\model.hmoc_MT3D\rmamthm.mnm:
                MT3DMS or MT3D-USGS Name file

            model\model.hmoc_MT3D\rmamthm.ssm:
                MT3DMS or MT3D-USGS Sink and Source Mixing Package input file

            model\model.hmoc_MT3D\rmamthmoc.adv:
                MT3DMS or MT3D-USGS Advection Package input file

            model\model.hmoc_MT3D\rmamthmoc.btn:
                MT3DMS or MT3D-USGS Basic Transport Package input file

            model\model.hmoc_MT3D\rmamthmoc.dsp:
                MT3DMS or MT3D-USGS Dispersion Package input file

            model\model.hmoc_MT3D\rmamthmoc.gcg:
                MT3DMS or MT3D-USGS Generalized Conjugate Gradient Solver
                Package input file

            model\model.hmoc_MT3D\rmamthmoc.mnm:
                MT3DMS or MT3D-USGS Name file

            model\model.hmoc_MT3D\rmamthmoc.ssm:
                MT3DMS or MT3D-USGS Sink and Source Mixing Package input file


        model\model.mmoc\
            Description:
            -----------
            Mmoc is a MODFLOW-2000 model. It is used in conjunction with
            mmoc_MT3D. It is one of several models used to compare results
            between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mf2k.exe from the bin directory into the this
            directory and double-click on it.

            Files:
            -----
            model\model.mmoc\MODFLOW.BAT:
                Batch file used to run a model

            model\model.mmoc\modflow.bf:
                MODFLOW-2000 batch input file

            model\model.mmoc\MT3DMS.BAT:
                Batch file used to run a model

            model\model.mmoc\rmamtmmoc.bas:
                MODFLOW Basic Package input file

            model\model.mmoc\rmamtmmoc.dis:
                MODFLOW Discretization file

            model\model.mmoc\rmamtmmoc.lmt:
                MODFLOW Link to MT3DMS or MT3D-USGS input file

            model\model.mmoc\rmamtmmoc.lpf:
                MODFLOW Layer Property Flow Package input file

            model\model.mmoc\rmamtmmoc.nam:
                MODFLOW Name file

            model\model.mmoc\rmamtmmoc.oc:
                MODFLOW Output Control input file

            model\model.mmoc\rmamtmmoc.pcg:
                MODFLOW Preconditioned Conjugate Gradient Package input file

            model\model.mmoc\rmamtmmoc.wel:
                MODFLOW Well Package input file


        model\model.mmoc_MT3D\
            Description:
            -----------
            Mmoc_MT3D is an MT3DMS model using the Modified MOC algorithm. It is
            used in conjunction with mmoc. It is one of several models used to
            compare results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mt3dms5b.exe from the bin directory into the
            this directory and the flow-transport link file (*.ftl) from the
            corresponding MODFLOW-2000 output directory and double-click on
            mt3dms5b.exe. At the prompt, enter the name of the MT3DMS name file
            (*.mnm).

            Files:
            -----
            model\model.mmoc_MT3D\rmamtmmoc.adv:
                MT3DMS or MT3D-USGS Advection Package input file

            model\model.mmoc_MT3D\rmamtmmoc.btn:
                MT3DMS or MT3D-USGS Basic Transport Package input file

            model\model.mmoc_MT3D\rmamtmmoc.dsp:
                MT3DMS or MT3D-USGS Dispersion Package input file

            model\model.mmoc_MT3D\rmamtmmoc.gcg:
                MT3DMS or MT3D-USGS Generalized Conjugate Gradient Solver
                Package input file

            model\model.mmoc_MT3D\rmamtmmoc.mnm:
                MT3DMS or MT3D-USGS Name file

            model\model.mmoc_MT3D\rmamtmmoc.ssm:
                MT3DMS or MT3D-USGS Sink and Source Mixing Package input file


        model\model.moc\
            Description:
            -----------
            Moc is a MODFLOW-2000 model. It is used in conjunction with
            moc_MT3D. It is one of several models used to compare results
            between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mf2k.exe from the bin directory into the this
            directory and double-click on it.

            Files:
            -----
            model\model.moc\MODFLOW.BAT:
                Batch file used to run a model

            model\model.moc\modflow.bf:
                MODFLOW-2000 batch input file

            model\model.moc\MT3DMS.BAT:
                Batch file used to run a model

            model\model.moc\rmamtmoc.bas:
                MODFLOW Basic Package input file

            model\model.moc\rmamtmoc.dis:
                MODFLOW Discretization file

            model\model.moc\rmamtmoc.lmt:
                MODFLOW Link to MT3DMS or MT3D-USGS input file

            model\model.moc\rmamtmoc.lpf:
                MODFLOW Layer Property Flow Package input file

            model\model.moc\rmamtmoc.nam:
                MODFLOW Name file

            model\model.moc\rmamtmoc.oc:
                MODFLOW Output Control input file

            model\model.moc\rmamtmoc.pcg:
                MODFLOW Preconditioned Conjugate Gradient Package input file

            model\model.moc\rmamtmoc.wel:
                MODFLOW Well Package input file


        model\model.moc_MT3D\
            Description:
            -----------
            Moc_MT3D is an MT3DMS model using the MOC algorithm. It is used in
            conjunction with moc. It is one of several models used to compare
            results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mt3dms5b.exe from the bin directory into the
            this directory and the flow-transport link file (*.ftl) from the
            corresponding MODFLOW-2000 output directory and double-click on
            mt3dms5b.exe. At the prompt, enter the name of the MT3DMS name file
            (*.mnm).

            Files:
            -----
            model\model.moc_MT3D\rmamtmoc.adv:
                MT3DMS or MT3D-USGS Advection Package input file

            model\model.moc_MT3D\rmamtmoc.btn:
                MT3DMS or MT3D-USGS Basic Transport Package input file

            model\model.moc_MT3D\rmamtmoc.dsp:
                MT3DMS or MT3D-USGS Dispersion Package input file

            model\model.moc_MT3D\rmamtmoc.gcg:
                MT3DMS or MT3D-USGS Generalized Conjugate Gradient Solver
                Package input file

            model\model.moc_MT3D\rmamtmoc.mnm:
                MT3DMS or MT3D-USGS Name file

            model\model.moc_MT3D\rmamtmoc.ssm:
                MT3DMS or MT3D-USGS Sink and Source Mixing Package input file


        model\model.tvd\
            Description:
            -----------
            Tvd is a MODFLOW-2000 model. It is used in conjunction with
            tvd_MT3D. It is one of several models used to compare results
            between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mf2k.exe from the bin directory into the this
            directory and double-click on it.

            Files:
            -----
            model\model.tvd\MODFLOW.BAT:
                Batch file used to run a model

            model\model.tvd\modflow.bf:
                MODFLOW-2000 batch input file

            model\model.tvd\MT3DMS.BAT:
                Batch file used to run a model

            model\model.tvd\rmamttvd.bas:
                MODFLOW Basic Package input file

            model\model.tvd\rmamttvd.dis:
                MODFLOW Discretization file

            model\model.tvd\rmamttvd.lmt:
                MODFLOW Link to MT3DMS or MT3D-USGS input file

            model\model.tvd\rmamttvd.lpf:
                MODFLOW Layer Property Flow Package input file

            model\model.tvd\rmamttvd.nam:
                MODFLOW Name file

            model\model.tvd\rmamttvd.oc:
                MODFLOW Output Control input file

            model\model.tvd\rmamttvd.pcg:
                MODFLOW Preconditioned Conjugate Gradient Package input file

            model\model.tvd\rmamttvd.wel:
                MODFLOW Well Package input file


        model\model.tvd_MT3D\
            Description:
            -----------
            Tvd_MT3D is an MT3DMS model using the total-variation-diminishing
            algorithm. It is used in conjunction with tvd. It is one of several
            models used to compare results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy mt3dms5b.exe from the bin directory into the
            this directory and the flow-transport link file (*.ftl) from the
            corresponding MODFLOW-2000 output directory and double-click on
            mt3dms5b.exe. At the prompt, enter the name of the MT3DMS name file
            (*.mnm).

            Files:
            -----
            model\model.tvd_MT3D\rmamttvd.adv:
                MT3DMS or MT3D-USGS Advection Package input file

            model\model.tvd_MT3D\rmamttvd.btn:
                MT3DMS or MT3D-USGS Basic Transport Package input file

            model\model.tvd_MT3D\rmamttvd.dsp:
                MT3DMS or MT3D-USGS Dispersion Package input file

            model\model.tvd_MT3D\rmamttvd.gcg:
                MT3DMS or MT3D-USGS Generalized Conjugate Gradient Solver
                Package input file

            model\model.tvd_MT3D\rmamttvd.mnm:
                MT3DMS or MT3D-USGS Name file

            model\model.tvd_MT3D\rmamttvd.ssm:
                MT3DMS or MT3D-USGS Sink and Source Mixing Package input file


        model\model.rmamocwti\
            Description:
            -----------
            Rmamocwti is a MODFLOW-GWT model using the MOCWTI method. It is one
            of several models used to compare results between MODFLOW-GWT and
            MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

            Files:
            -----
            model\model.rmamocwti\MF2K_GWT.BAT:
                Batch file used to run a model

            model\model.rmamocwti\modflow.bf:
                MODFLOW-2000 batch input file

            model\model.rmamocwti\rmamocwti.bas:
                MODFLOW Basic Package input file

            model\model.rmamocwti\rmamocwti.dis:
                MODFLOW Discretization file

            model\model.rmamocwti\rmamocwti.gwt:
                MODFLOW-GWT Transport name file

            model\model.rmamocwti\rmamocwti.lpf:
                MODFLOW Layer Property Flow Package input file

            model\model.rmamocwti\rmamocwti.mob:
                MODFLOW-GWT Obsesrvation input file

            model\model.rmamocwti\rmamocwti.moc:
                MODFLOW-GWT Transport input file

            model\model.rmamocwti\rmamocwti.nam:
                MODFLOW Name file

            model\model.rmamocwti\rmamocwti.oc:
                MODFLOW Output Control input file

            model\model.rmamocwti\rmamocwti.pcg:
                MODFLOW Preconditioned Conjugate Gradient Package input file

            model\model.rmamocwti\rmamocwti.prtp:
                MODFLOW-GWT PTRP input file

            model\model.rmamocwti\rmamocwti.wel:
                MODFLOW Well Package input file


        model\model.DiracUnidirectional\
            Description:
            -----------
            DiracUnidirectional is a MODFLOW-GWT model using the MOCWTI method.
            Its purpose is to compare the numerical solution with an analytical
            solution to the same problem.
            This folder contains input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

            Files:
            -----
            model\model.DiracUnidirectional\D-Uni.bas:
                MODFLOW Basic Package input file

            model\model.DiracUnidirectional\D-Uni.dis:
                MODFLOW Discretization file

            model\model.DiracUnidirectional\D-Uni.gwt:
                MODFLOW-GWT Transport name file

            model\model.DiracUnidirectional\D-Uni.ipdl:
                MODFLOW-GWT Initial Particle Density File-List-Based input file

            model\model.DiracUnidirectional\D-Uni.lpf:
                MODFLOW Layer Property Flow Package input file

            model\model.DiracUnidirectional\D-Uni.mob:
                MODFLOW-GWT Obsesrvation input file

            model\model.DiracUnidirectional\D-Uni.moc:
                MODFLOW-GWT Transport input file

            model\model.DiracUnidirectional\D-Uni.nam:
                MODFLOW Name file

            model\model.DiracUnidirectional\D-Uni.oc:
                MODFLOW Output Control input file

            model\model.DiracUnidirectional\D-Uni.pcg:
                MODFLOW Preconditioned Conjugate Gradient Package input file

            model\model.DiracUnidirectional\MF2K_GWT.BAT:
                Batch file used to run a model

            model\model.DiracUnidirectional\modflow.bf:
                MODFLOW-2000 batch input file


        model\model.DiracAngle\
            Description:
            -----------
            DiracAngle is a MODFLOW-GWT model using the MOCWTI method. Its
            purpose is to compare the numerical solution with an analytical
            solution to the same problem. It differs from DiracUnidirectional in
            that flow is at a 45 degree angle to the grid instead of parallel to
            the rows.
            This folder contains input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

            Files:
            -----
            model\model.DiracAngle\D-Ang.bas:
                MODFLOW Basic Package input file

            model\model.DiracAngle\D-Ang.dis:
                MODFLOW Discretization file

            model\model.DiracAngle\D-Ang.gwt:
                MODFLOW-GWT Transport name file

            model\model.DiracAngle\D-Ang.ipdl:
                MODFLOW-GWT Initial Particle Density File-List-Based input file

            model\model.DiracAngle\D-Ang.lpf:
                MODFLOW Layer Property Flow Package input file

            model\model.DiracAngle\D-Ang.mob:
                MODFLOW-GWT Obsesrvation input file

            model\model.DiracAngle\D-Ang.moc:
                MODFLOW-GWT Transport input file

            model\model.DiracAngle\D-Ang.nam:
                MODFLOW Name file

            model\model.DiracAngle\D-Ang.oc:
                MODFLOW Output Control input file

            model\model.DiracAngle\D-Ang.pcg:
                MODFLOW Preconditioned Conjugate Gradient Package input file

            model\model.DiracAngle\MF2K_GWT.BAT:
                Batch file used to run a model

            model\model.DiracAngle\modflow.bf:
                MODFLOW-2000 batch input file


        model\model.Point3\
            Description:
            -----------
            Point3 is a MODFLOW-GWT model of a continuous source in a uniform
            three-dimensional flow field. It is used for testing MODFLOW-GWT by
            comparing its results with the Point3Analytical model in this
            archive. This folder contains the input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

            Files:
            -----
            model\model.Point3\MF2K_GWT.BAT:
                Batch file used to run a model

            model\model.Point3\modflow.bf:
                MODFLOW-2000 batch input file

            model\model.Point3\point3.bas:
                MODFLOW Basic Package input file

            model\model.Point3\point3.bcf:
                MODFLOW Block-Centered Flow Package input file

            model\model.Point3\point3.dis:
                MODFLOW Discretization file

            model\model.Point3\point3.gwt:
                MODFLOW-GWT Transport name file

            model\model.Point3\point3.ipdl:
                MODFLOW-GWT Initial Particle Density File-List-Based input file

            model\model.Point3\point3.mob:
                MODFLOW-GWT Obsesrvation input file

            model\model.Point3\point3.moc:
                MODFLOW-GWT Transport input file

            model\model.Point3\point3.nam:
                MODFLOW Name file

            model\model.Point3\point3.oc:
                MODFLOW Output Control input file

            model\model.Point3\point3.prtp:
                MODFLOW-GWT PTRP input file

            model\model.Point3\point3.sip:
                MODFLOW Strongly Implicit Procedure Package input file

            model\model.Point3\point3.wel:
                MODFLOW Well Package input file


        model\model.Finite1\
            Description:
            -----------
            This model is an analytical model documented in Wexler (1992) and
            available from http://water.usgs.gov/software/ANALGWST/. The program
            models a continuous solute source in a one-dimensional uniform flow
            field. The particular program in Wexler (1992) used for this model
            is named Finite. This folder contains the input files for the model.
            The output of this file is used for comparison with the Finite.1D
            MODFLOW-GWT model in this archive. Another similar analytical Finite
            model in this archive is Finite2 which has a higher longidudinal
            dispersivity.
            To run the model, copy Finite.exe from the bin directory into this
            directory and double-click on it. Then specify the Finite1.dat as
            the input file.

            Files:
            -----
            model\model.Finite1\Finite1.dat:
                Input file for the Finite program.


        model\model.Finite2\
            Description:
            -----------
            This model is an analytical model documented in Wexler (1992) and
            available from http://water.usgs.gov/software/ANALGWST/. The program
            models a continuous solute source in a one-dimensional uniform flow
            field. The particular program in Wexler (1992) used for this model
            is named Finite. This folder contains the input files for the model.
            The output of this file is used for comparison with the Finite.1D
            MODFLOW-GWT model in this archive. Another similar analytical Finite
            model in this archive is Finite1 which has a lower longidudinal
            dispersivity.
            To run the model, copy Finite.exe from the bin directory into this
            directory and double-click on it. Then specify the Finite2.dat as
            the input file.

            Files:
            -----
            model\model.Finite2\Finite2.dat:
                Input file for the Finite program.


        model\model.Point3Analytical\
            Description:
            -----------
            This model is an analytical model documented in Wexler (1992) and
            available from http://water.usgs.gov/software/ANALGWST/. The program
            models a continuous solute source in a three-dimensional uniform
            flow field. The particular program in Wexler (1992) used for this
            model is named Point3. This folder contains the input files for the
            model. The output of this file is used for comparison with the
            Point3 MODFLOW-GWT model in this archive.
            To run the model, copy Point3.exe from the bin directory into this
            directory and double-click on it. Then specify the moc9.z025.dat as
            the input file.

            Files:
            -----
            model\model.Point3Analytical\InputFiles.txt:
                InputFiles.txt is used in RunPoint3.bat to specify the responses
                for file name requests.

            model\model.Point3Analytical\moc9.z025.dat:
                Moc9.z025.dat is the input file for the Point3 analytical model.

            model\model.Point3Analytical\RunPoint3.bat:
                Batch file used to run the Point3 model.


        model\model.Dirac\
            Description:
            -----------
            Dirac contains the input for an analytical solution to the problem
            of an instantaneous solute source in a uniform 3D flow field. The
            code was originally described in Konikow and others (1996). The
            analytical solution is from Wexler (1992). The purpose of the model
            is for comparison with the DiracUnidirectional and DiracAngle
            MODFLOW-GWT models.
            To run the model, copy Dirac.exe into this folder and double-click
            it. select uni.3d-120D.a-L1z2.5.dat as the input file.

            Files:
            -----
            model\model.Dirac\uni.3d.120D.a-L1z2.5.dat:
                uni.3d-120D.a-L1z2.5.dat is the input file for the program Dirac
                which simulates and instantaneous solute point source.


        model\model.FiniteAlpha1\
            Description:
            -----------
            FiniteAlpha1 is a one-dimensional MODFLOW-GWT model. Flow in the
            model is steady-state. Transport is transient. It has three
            observation locations. The purpose of the model is to test the MOCWT
            algorithm by comparing the simulated concentrations at the
            observation locations with an analytical solution. The analytical
            solution is included in this archive as Finite2. Another similar
            MODFLOW-GWT model in this archive is Finite.1D which has a lower
            longitudinal dispersivity.
            This folder contains input files for the model.
            To run the model, copy MOCWT.exe from the bin directory into the
            this directory and double-click on it.

            Files:
            -----
            model\model.FiniteAlpha1\finiteAlpha1.bas:
                MODFLOW Basic Package input file

            model\model.FiniteAlpha1\finiteAlpha1.bcf:
                MODFLOW Block-Centered Flow Package input file

            model\model.FiniteAlpha1\finiteAlpha1.dis:
                MODFLOW Discretization file

            model\model.FiniteAlpha1\finiteAlpha1.gwt:
                MODFLOW-GWT Transport name file

            model\model.FiniteAlpha1\finiteAlpha1.mob:
                MODFLOW-GWT Obsesrvation input file

            model\model.FiniteAlpha1\finiteAlpha1.moc:
                MODFLOW-GWT Transport input file

            model\model.FiniteAlpha1\finiteAlpha1.nam:
                MODFLOW Name file

            model\model.FiniteAlpha1\finiteAlpha1.oc:
                MODFLOW Output Control input file

            model\model.FiniteAlpha1\finiteAlpha1.sip:
                MODFLOW Strongly Implicit Procedure Package input file

            model\model.FiniteAlpha1\MF2K_GWT.BAT:
                Batch file used to run a model

            model\model.FiniteAlpha1\modflow.bf:
                MODFLOW-2000 batch input file


        output\output.Finite.1D\
            Description:
            -----------
            Finite.1D is a one-dimensional MODFLOW-GWT model. Flow in the model
            is steady-state. Transport is transient. It has three observation
            locations. The purpose of the model is to test the MOCWT algorithm
            by comparing the simulated concentrations at the observation
            locations with an analytical solution.
            This folder contains output files for the model.

            Files:
            -----
            output\output.Finite.1D\finite.bud:
                MODFLOW Cell-By-Cell flow file

            output\output.Finite.1D\finite.cna:
                MODFLOW-GWT ASCII concentration output file

            output\output.Finite.1D\finite.fdn:
                MODFLOW Formatted Drawdown file

            output\output.Finite.1D\finite.fhd:
                MODFLOW Formatted Head file

            output\output.Finite.1D\finite.lst:
                MODFLOW Listing file

            output\output.Finite.1D\finite.oba:
                MODFLOW-GWT Obsesrvation output file

            output\output.Finite.1D\finite.out:
                SUTRA main output file

            output\output.Finite.1D\finite.pta:
                MODFLOW-GWT ASCII particle output file

            output\output.Finite.1D\finite.vla:
                MODFLOW-GWT ASCII velocity output file

            output\output.Finite.1D\mf2kerr.p00:
                MODFLOW-2000 error output file

            output\output.Finite.1D\modbatch.rpt:
                MODFLOW-2000 batch report file


        output\output.rmamocel\
            Description:
            -----------
            Rmamocel is a MODFLOW-GWT model using the ELLAM method. It is one of
            several models used to compare results between MODFLOW-GWT and
            MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

            Files:
            -----
            output\output.rmamocel\mf2kerr.p00:
                MODFLOW-2000 error output file

            output\output.rmamocel\modbatch.rpt:
                MODFLOW-2000 batch report file

            output\output.rmamocel\rmamocel.bud:
                MODFLOW Cell-By-Cell flow file

            output\output.rmamocel\rmamocel.cna:
                MODFLOW-GWT ASCII concentration output file

            output\output.rmamocel\rmamocel.fdn:
                MODFLOW Formatted Drawdown file

            output\output.rmamocel\rmamocel.fhd:
                MODFLOW Formatted Head file

            output\output.rmamocel\rmamocel.lst:
                MODFLOW Listing file

            output\output.rmamocel\rmamocel.oba:
                MODFLOW-GWT Obsesrvation output file

            output\output.rmamocel\rmamocel.out:
                SUTRA main output file

            output\output.rmamocel\rmamocel.vla:
                MODFLOW-GWT ASCII velocity output file


        output\output.rmamocim\
            Description:
            -----------
            Rmamocim is a MODFLOW-GWT model using the MOCIMP method. It is one
            of several models used to compare results between MODFLOW-GWT and
            MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

            Files:
            -----
            output\output.rmamocim\mf2kerr.p00:
                MODFLOW-2000 error output file

            output\output.rmamocim\modbatch.rpt:
                MODFLOW-2000 batch report file

            output\output.rmamocim\rmamocim.bud:
                MODFLOW Cell-By-Cell flow file

            output\output.rmamocim\rmamocim.cna:
                MODFLOW-GWT ASCII concentration output file

            output\output.rmamocim\rmamocim.fdn:
                MODFLOW Formatted Drawdown file

            output\output.rmamocim\rmamocim.fhd:
                MODFLOW Formatted Head file

            output\output.rmamocim\rmamocim.lst:
                MODFLOW Listing file

            output\output.rmamocim\rmamocim.oba:
                MODFLOW-GWT Obsesrvation output file

            output\output.rmamocim\rmamocim.out:
                SUTRA main output file

            output\output.rmamocim\rmamocim.vla:
                MODFLOW-GWT ASCII velocity output file


        output\output.rmamocwti\
            Description:
            -----------
            Rmamocwti is a MODFLOW-GWT model using the MOCWTI method. It is one
            of several models used to compare results between MODFLOW-GWT and
            MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

            Files:
            -----
            output\output.rmamocwti\mf2kerr.p00:
                MODFLOW-2000 error output file

            output\output.rmamocwti\modbatch.rpt:
                MODFLOW-2000 batch report file

            output\output.rmamocwti\rmamocwti.bud:
                MODFLOW Cell-By-Cell flow file

            output\output.rmamocwti\rmamocwti.cna:
                MODFLOW-GWT ASCII concentration output file

            output\output.rmamocwti\rmamocwti.fdn:
                MODFLOW Formatted Drawdown file

            output\output.rmamocwti\rmamocwti.fhd:
                MODFLOW Formatted Head file

            output\output.rmamocwti\rmamocwti.lst:
                MODFLOW Listing file

            output\output.rmamocwti\rmamocwti.oba:
                MODFLOW-GWT Obsesrvation output file

            output\output.rmamocwti\rmamocwti.out:
                SUTRA main output file

            output\output.rmamocwti\rmamocwti.pta:
                MODFLOW-GWT ASCII particle output file

            output\output.rmamocwti\rmamocwti.vla:
                MODFLOW-GWT ASCII velocity output file


        output\output.rmamocwt_dist\
            Description:
            -----------
            Rmamocwt_dist is a MODFLOW-GWT model using the MOCWT method. It is
            one of several models used to compare results between MODFLOW-GWT
            and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

            Files:
            -----
            output\output.rmamocwt_dist\mf2kerr.p00:
                MODFLOW-2000 error output file

            output\output.rmamocwt_dist\modbatch.rpt:
                MODFLOW-2000 batch report file

            output\output.rmamocwt_dist\rmamocwt.bud:
                MODFLOW Cell-By-Cell flow file

            output\output.rmamocwt_dist\rmamocwt.cna:
                MODFLOW-GWT ASCII concentration output file

            output\output.rmamocwt_dist\rmamocwt.fdn:
                MODFLOW Formatted Drawdown file

            output\output.rmamocwt_dist\rmamocwt.fhd:
                MODFLOW Formatted Head file

            output\output.rmamocwt_dist\rmamocwt.lst:
                MODFLOW Listing file

            output\output.rmamocwt_dist\rmamocwt.oba:
                MODFLOW-GWT Obsesrvation output file

            output\output.rmamocwt_dist\rmamocwt.out:
                SUTRA main output file

            output\output.rmamocwt_dist\rmamocwt.pta:
                MODFLOW-GWT ASCII particle output file

            output\output.rmamocwt_dist\rmamocwt.vla:
                MODFLOW-GWT ASCII velocity output file

            output\output.rmamocwt_dist\rmamocwt_dist.bud:
                MODFLOW Cell-By-Cell flow file

            output\output.rmamocwt_dist\rmamocwt_dist.cna:
                MODFLOW-GWT ASCII concentration output file

            output\output.rmamocwt_dist\rmamocwt_dist.fdn:
                MODFLOW Formatted Drawdown file

            output\output.rmamocwt_dist\rmamocwt_dist.fhd:
                MODFLOW Formatted Head file

            output\output.rmamocwt_dist\rmamocwt_dist.lst:
                MODFLOW Listing file

            output\output.rmamocwt_dist\rmamocwt_dist.oba:
                MODFLOW-GWT Obsesrvation output file

            output\output.rmamocwt_dist\rmamocwt_dist.out:
                SUTRA main output file

            output\output.rmamocwt_dist\rmamocwt_dist.pta:
                MODFLOW-GWT ASCII particle output file

            output\output.rmamocwt_dist\rmamocwt_dist.vla:
                MODFLOW-GWT ASCII velocity output file


        output\output.rma_moc\
            Description:
            -----------
            Rma_moc is a MODFLOW-GWT model using the MOC method. It is one of
            several models used to compare results between MODFLOW-GWT and
            MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

            Files:
            -----
            output\output.rma_moc\mf2kerr.p00:
                MODFLOW-2000 error output file

            output\output.rma_moc\modbatch.rpt:
                MODFLOW-2000 batch report file

            output\output.rma_moc\rmamocex.bud:
                MODFLOW Cell-By-Cell flow file

            output\output.rma_moc\rmamocex.cna:
                MODFLOW-GWT ASCII concentration output file

            output\output.rma_moc\rmamocex.fdn:
                MODFLOW Formatted Drawdown file

            output\output.rma_moc\rmamocex.fhd:
                MODFLOW Formatted Head file

            output\output.rma_moc\rmamocex.lst:
                MODFLOW Listing file

            output\output.rma_moc\rmamocex.oba:
                MODFLOW-GWT Obsesrvation output file

            output\output.rma_moc\rmamocex.out:
                SUTRA main output file

            output\output.rma_moc\rmamocex.vla:
                MODFLOW-GWT ASCII velocity output file

            output\output.rma_moc\rmamocwt.pta:
                MODFLOW-GWT ASCII particle output file


        output\output.fd\
            Description:
            -----------
            Fd is a MODFLOW-2000 model. It is used in conjunction with fd_MT3D.
            It is one of several models used to compare results between
            MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

            Files:
            -----
            output\output.fd\mf2kerr.p00:
                MODFLOW-2000 error output file

            output\output.fd\modbatch.rpt:
                MODFLOW-2000 batch report file

            output\output.fd\rmamtfd.bud:
                MODFLOW Cell-By-Cell flow file

            output\output.fd\rmamtfd.fhd:
                MODFLOW Formatted Head file

            output\output.fd\rmamtfd.ftl:
                MODFLOW MT3DMS or MT3D-USGS Flow Transport Link output file

            output\output.fd\rmamtfd.lst:
                MODFLOW Listing file


        output\output.fd_MT3D\
            Description:
            -----------
            Fd_MT3D is an MT3DMS model using the Finite-Difference algorithm. It
            is used in conjunction with fd. It is one of several models used to
            compare results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

            Files:
            -----
            output\output.fd_MT3D\rmamtfd.cnf:
                MT3DMS or MT3D-USGS Grid Configuration output file

            output\output.fd_MT3D\rmamtfd.mls:
                MT3DMS or MT3D-USGS listing file

            output\output.fd_MT3D\rmamtfdA.mass:
                MT3DMS or MT3D-USGS mass balance summary file

            output\output.fd_MT3D\rmamtfdA.mto:
                MT3DMS or MT3D-USGS Observation Output file

            output\output.fd_MT3D\rmamtfdA.ucn:
                MT3DMS or MT3D-USGS Concentration output file


        output\output.hmoc\
            Description:
            -----------
            Hmoc is a MODFLOW-2000 model. It is used in conjunction with
            hmoc_MT3D. It is one of several models used to compare results
            between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

            Files:
            -----
            output\output.hmoc\mf2kerr.p00:
                MODFLOW-2000 error output file

            output\output.hmoc\modbatch.rpt:
                MODFLOW-2000 batch report file

            output\output.hmoc\rmamthm.bud:
                MODFLOW Cell-By-Cell flow file

            output\output.hmoc\rmamthm.fhd:
                MODFLOW Formatted Head file

            output\output.hmoc\rmamthm.ftl:
                MODFLOW MT3DMS or MT3D-USGS Flow Transport Link output file

            output\output.hmoc\rmamthm.lst:
                MODFLOW Listing file

            output\output.hmoc\rmamthmoc.bud:
                MODFLOW Cell-By-Cell flow file

            output\output.hmoc\rmamthmoc.fhd:
                MODFLOW Formatted Head file

            output\output.hmoc\rmamthmoc.ftl:
                MODFLOW MT3DMS or MT3D-USGS Flow Transport Link output file

            output\output.hmoc\rmamthmoc.lst:
                MODFLOW Listing file


        output\output.hmoc_MT3D\
            Description:
            -----------
            Hmoc_MT3D is an MT3DMS model using the Hybrid-MOC algorithm. It is
            used in conjunction with hmoc. It is one of several models used to
            compare results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

            Files:
            -----
            output\output.hmoc_MT3D\rmamthm.cnf:
                MT3DMS or MT3D-USGS Grid Configuration output file

            output\output.hmoc_MT3D\rmamthm.mls:
                MT3DMS or MT3D-USGS listing file

            output\output.hmoc_MT3D\rmamthmA.mass:
                MT3DMS or MT3D-USGS mass balance summary file

            output\output.hmoc_MT3D\rmamthmA.mto:
                MT3DMS or MT3D-USGS Observation Output file

            output\output.hmoc_MT3D\rmamthmA.ucn:
                MT3DMS or MT3D-USGS Concentration output file

            output\output.hmoc_MT3D\rmamthmoc.cnf:
                MT3DMS or MT3D-USGS Grid Configuration output file

            output\output.hmoc_MT3D\rmamthmoc.mls:
                MT3DMS or MT3D-USGS listing file


        output\output.mmoc\
            Description:
            -----------
            Mmoc is a MODFLOW-2000 model. It is used in conjunction with
            mmoc_MT3D. It is one of several models used to compare results
            between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains input files for the model.

            Files:
            -----
            output\output.mmoc\mf2kerr.p00:
                MODFLOW-2000 error output file

            output\output.mmoc\modbatch.rpt:
                MODFLOW-2000 batch report file

            output\output.mmoc\rmamtmmoc.bud:
                MODFLOW Cell-By-Cell flow file

            output\output.mmoc\rmamtmmoc.fhd:
                MODFLOW Formatted Head file

            output\output.mmoc\rmamtmmoc.ftl:
                MODFLOW MT3DMS or MT3D-USGS Flow Transport Link output file

            output\output.mmoc\rmamtmmoc.lst:
                MODFLOW Listing file


        output\output.mmoc_MT3D\
            Description:
            -----------
            Mmoc_MT3D is an MT3DMS model using the Modified MOC algorithm. It is
            used in conjunction with mmoc. It is one of several models used to
            compare results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

            Files:
            -----
            output\output.mmoc_MT3D\rmamtmmA.mass:
                MT3DMS or MT3D-USGS mass balance summary file

            output\output.mmoc_MT3D\rmamtmmA.mto:
                MT3DMS or MT3D-USGS Observation Output file

            output\output.mmoc_MT3D\rmamtmmA.ucn:
                MT3DMS or MT3D-USGS Concentration output file

            output\output.mmoc_MT3D\rmamtmmoc.cnf:
                MT3DMS or MT3D-USGS Grid Configuration output file

            output\output.mmoc_MT3D\rmamtmmoc.mls:
                MT3DMS or MT3D-USGS listing file


        output\output.moc\
            Description:
            -----------
            Moc is a MODFLOW-2000 model. It is used in conjunction with
            moc_MT3D. It is one of several models used to compare results
            between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

            Files:
            -----
            output\output.moc\mf2kerr.p00:
                MODFLOW-2000 error output file

            output\output.moc\modbatch.rpt:
                MODFLOW-2000 batch report file

            output\output.moc\rmamtmoc.bud:
                MODFLOW Cell-By-Cell flow file

            output\output.moc\rmamtmoc.fhd:
                MODFLOW Formatted Head file

            output\output.moc\rmamtmoc.ftl:
                MODFLOW MT3DMS or MT3D-USGS Flow Transport Link output file

            output\output.moc\rmamtmoc.lst:
                MODFLOW Listing file


        output\output.moc_MT3D\
            Description:
            -----------
            Moc_MT3D is an MT3DMS model using the MOC algorithm. It is used in
            conjunction with moc. It is one of several models used to compare
            results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

            Files:
            -----
            output\output.moc_MT3D\rmamtmoA.mass:
                MT3DMS or MT3D-USGS mass balance summary file

            output\output.moc_MT3D\rmamtmoA.mto:
                MT3DMS or MT3D-USGS Observation Output file

            output\output.moc_MT3D\rmamtmoA.ucn:
                MT3DMS or MT3D-USGS Concentration output file

            output\output.moc_MT3D\rmamtmoc.cnf:
                MT3DMS or MT3D-USGS Grid Configuration output file

            output\output.moc_MT3D\rmamtmoc.mls:
                MT3DMS or MT3D-USGS listing file


        output\output.tvd\
            Description:
            -----------
            Tvd is a MODFLOW-2000 model. It is used in conjunction with
            tvd_MT3D. It is one of several models used to compare results
            between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

            Files:
            -----
            output\output.tvd\mf2kerr.p00:
                MODFLOW-2000 error output file

            output\output.tvd\modbatch.rpt:
                MODFLOW-2000 batch report file

            output\output.tvd\rmamttvd.bud:
                MODFLOW Cell-By-Cell flow file

            output\output.tvd\rmamttvd.fhd:
                MODFLOW Formatted Head file

            output\output.tvd\rmamttvd.ftl:
                MODFLOW MT3DMS or MT3D-USGS Flow Transport Link output file

            output\output.tvd\rmamttvd.lst:
                MODFLOW Listing file


        output\output.tvd_MT3D\
            Description:
            -----------
            Tvd_MT3D is an MT3DMS model using the total-variation-diminishing
            algorithm. It is used in conjunction with tvd. It is one of several
            models used to compare results between MODFLOW-GWT and MT3DMS. 
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977). 
            This folder contains output files for the model.

            Files:
            -----
            output\output.tvd_MT3D\rmamttvA.mass:
                MT3DMS or MT3D-USGS mass balance summary file

            output\output.tvd_MT3D\rmamttvA.mto:
                MT3DMS or MT3D-USGS Observation Output file

            output\output.tvd_MT3D\rmamttvA.ucn:
                MT3DMS or MT3D-USGS Concentration output file

            output\output.tvd_MT3D\rmamttvd.cnf:
                MT3DMS or MT3D-USGS Grid Configuration output file

            output\output.tvd_MT3D\rmamttvd.mls:
                MT3DMS or MT3D-USGS listing file


        output\output.DiracUnidirectional\
            Description:
            -----------
            DiracUnidirectional is a MODFLOW-GWT model using the MOCWTI method.
            Its purpose is to compare the numerical solution with an analytical
            solution to the same problem.
            This folder contains output files for the model.

            Files:
            -----
            output\output.DiracUnidirectional\D-Uni.bud:
                MODFLOW Cell-By-Cell flow file

            output\output.DiracUnidirectional\D-Uni.cna:
                MODFLOW-GWT ASCII concentration output file

            output\output.DiracUnidirectional\D-Uni.fhd:
                MODFLOW Formatted Head file

            output\output.DiracUnidirectional\D-Uni.lst:
                MODFLOW Listing file

            output\output.DiracUnidirectional\D-Uni.out:
                SUTRA main output file

            output\output.DiracUnidirectional\D-Uni.vla:
                MODFLOW-GWT ASCII velocity output file

            output\output.DiracUnidirectional\mf2kerr.p00:
                MODFLOW-2000 error output file

            output\output.DiracUnidirectional\modbatch.rpt:
                MODFLOW-2000 batch report file


        output\output.DiracAngle\
            Description:
            -----------
            DiracAngle is a MODFLOW-GWT model using the MOCWTI method. Its
            purpose is to compare the numerical solution with an analytical
            solution to the same problem. It differs from DiracUnidirectional in
            that flow is at a 45 degree angle to the grid instead of parallel to
            the rows.
            This folder contains output files for the model.

            Files:
            -----
            output\output.DiracAngle\D-Ang.bud:
                MODFLOW Cell-By-Cell flow file

            output\output.DiracAngle\D-Ang.cna:
                MODFLOW-GWT ASCII concentration output file

            output\output.DiracAngle\D-Ang.fhd:
                MODFLOW Formatted Head file

            output\output.DiracAngle\D-Ang.lst:
                MODFLOW Listing file

            output\output.DiracAngle\D-Ang.out:
                SUTRA main output file

            output\output.DiracAngle\D-Ang.vla:
                MODFLOW-GWT ASCII velocity output file

            output\output.DiracAngle\mf2kerr.p00:
                MODFLOW-2000 error output file

            output\output.DiracAngle\modbatch.rpt:
                MODFLOW-2000 batch report file


        output\output.Point3\
            Description:
            -----------
            Point3 is a MODFLOW-GWT model of a continuous source in a uniform
            three-dimensional flow field. It is used for testing MODFLOW-GWT by
            comparing its results with the Point3Analytical model in this
            archive. This folder contains the out
            put files for the model.

            Files:
            -----
            output\output.Point3\mf2kerr.p00:
                MODFLOW-2000 error output file

            output\output.Point3\modbatch.rpt:
                MODFLOW-2000 batch report file

            output\output.Point3\point3.bud:
                MODFLOW Cell-By-Cell flow file

            output\output.Point3\point3.cna:
                MODFLOW-GWT ASCII concentration output file

            output\output.Point3\point3.fdn:
                MODFLOW Formatted Drawdown file

            output\output.Point3\point3.fhd:
                MODFLOW Formatted Head file

            output\output.Point3\point3.lst:
                MODFLOW Listing file

            output\output.Point3\point3.mbrp:
                MODFLOW-GWT Generalized Mass Balance output file

            output\output.Point3\point3.oba:
                MODFLOW-GWT Obsesrvation output file

            output\output.Point3\point3.out:
                SUTRA main output file

            output\output.Point3\point3.ptb:
                MODFLOW-GWT Binary particle output file

            output\output.Point3\point3.vla:
                MODFLOW-GWT ASCII velocity output file


        output\output.Finite1\
            Description:
            -----------
            This model is an analytical model documented in Wexler (1992) and
            available from http://water.usgs.gov/software/ANALGWST/. The program
            models a continuous solute source in a one-dimensional uniform flow
            field. The particular program in Wexler (1992) used for this model
            is named Finite. This folder contains the output file for the model.
            The output of this file is used for comparison with the Finite.1D
            MODFLOW-GWT model in this archive. Another similar analytical Finite
            model in this archive is Finite2 which has a higher longidudinal
            dispersivity.

            Files:
            -----
            output\output.Finite1\Finite1.prt:
                Finite program output file used for comparison with the
                Finite.1D MODFLOW-GWT model.


        output\output.Finite2\
            Description:
            -----------
            This model is an analytical model documented in Wexler (1992) and
            available from http://water.usgs.gov/software/ANALGWST/. The program
            models a continuous solute source in a one-dimensional uniform flow
            field. The particular program in Wexler (1992) used for this model
            is named Finite. This folder contains the output file for the model.
            The output of this file is used for comparison with the Finite.1D
            MODFLOW-GWT model in this archive. Another similar analytical Finite
            model in this archive is Finite1 which has a lower longidudinal
            dispersivity.

            Files:
            -----
            output\output.Finite2\Finite2.prt:
                Finite program output file used for comparison with the
                Finite.1D MODFLOW-GWT model.


        output\output.Point3Analytical\
            Description:
            -----------
            This model is an analytical model documented in Wexler (1992) and
            available from http://water.usgs.gov/software/ANALGWST/.  The
            particular program in Wexler (1992) used for this model is named
            Point3. This folder contains the output files for the model. The
            output of this file is used for comparison with the Point3
            MODFLOW-GWT model in this archive.

            Files:
            -----
            output\output.Point3Analytical\moc9.z025.out:
                This is the output of the Point3Analytical model.


        output\output.Dirac\
            Description:
            -----------
            Dirac contains the output for an analytical solution to the problem
            of an instantaneous solute source in a uniform 3D flow field. The
            code was originally described in Konikow and others (1996). The
            analytical solution is from Wexler (1992). The purpose of the model
            is for comparison with the DiracUnidirectional and DiracAngle
            MODFLOW-GWT models.

            Files:
            -----
            output\output.Dirac\CDUMP.OUT:
                Binary output file from Dirac4 program.

            output\output.Dirac\clog.out:
                Log output file from Dirac4 program.

            output\output.Dirac\uni.3d.120D.a-L1z2.5.out:
                Main output file from the Dirac program.


        output\output.FiniteAlpha1\
            Description:
            -----------
            FiniteAlpha1 is a one-dimensional MODFLOW-GWT model. Flow in the
            model is steady-state. Transport is transient. It has three
            observation locations. The purpose of the model is to test the MOCWT
            algorithm by comparing the simulated concentrations at the
            observation locations with an analytical solution. The analytical
            solution is included in this archive as Finite2. Another similar
            MODFLOW-GWT model in this archive is Finite.1D which has a lower
            longitudinal dispersivity.
            This folder contains output files for the model.

            Files:
            -----
            output\output.FiniteAlpha1\finiteAlpha1.bud:
                MODFLOW Cell-By-Cell flow file

            output\output.FiniteAlpha1\finiteAlpha1.cna:
                MODFLOW-GWT ASCII concentration output file

            output\output.FiniteAlpha1\finiteAlpha1.fdn:
                MODFLOW Formatted Drawdown file

            output\output.FiniteAlpha1\finiteAlpha1.fhd:
                MODFLOW Formatted Head file

            output\output.FiniteAlpha1\finiteAlpha1.lst:
                MODFLOW Listing file

            output\output.FiniteAlpha1\finiteAlpha1.oba:
                MODFLOW-GWT Obsesrvation output file

            output\output.FiniteAlpha1\finiteAlpha1.out:
                SUTRA main output file

            output\output.FiniteAlpha1\finiteAlpha1.pta:
                MODFLOW-GWT ASCII particle output file

            output\output.FiniteAlpha1\finiteAlpha1.vla:
                MODFLOW-GWT ASCII velocity output file

            output\output.FiniteAlpha1\mf2kerr.p00:
                MODFLOW-2000 error output file

            output\output.FiniteAlpha1\modbatch.rpt:
                MODFLOW-2000 batch report file


    \source\
        Description:
        -----------
        Source code for the models used to run the model simulations.
            source\Dirac\dirac4.f:
                Dirac4.f is Fortran source code for a program giving an
                analytical solution to the problem of an instantaneous solute
                source in a uniform 3D flow field. The code was originally
                described in Konikow and others (1996). The analytical solution
                is from Wexler (1992).

            source\Finite\dimens.inc:
                Include file for Finite as documented in Wexler (1992).

            source\Finite\exerfc.f:
                Fortran Source code file for Finite as documented in Wexler
                (1992).

            source\Finite\finite.f:
                Fortran Source code file for Finite as documented in Wexler
                (1992).

            source\Finite\lenchr.f:
                Fortran Source code file for Finite as documented in Wexler
                (1992).

            source\Finite\locsubs.intel.f:
                Fortran Source code file for Finite as documented in Wexler
                (1992). This file has been modified to include compiler-specific
                functions for date and time for the Intel fortran compiler.

            source\Finite\ofile.f:
                Fortran Source code file for Finite as documented in Wexler
                (1992).

            source\Finite\pltdat.inc:
                Include file for Finite as documented in Wexler (1992).

            source\Finite\subs1.f:
                Fortran Source code file for Finite as documented in Wexler
                (1992).

            source\Finite\title.f:
                Fortran Source code file for Finite as documented in Wexler
                (1992).

            source\Finite\units.inc:
                Include file for Finite as documented in Wexler (1992).

            source\MODFLOW-2000\ccfd.c:
                C source code file for MODFLOW-2000.

            source\MODFLOW-2000\ccfd.h:
                C header file for MODFLOW-2000.

            source\MODFLOW-2000\ctime.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\daf1.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\de45.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\glo1bas6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gmg1.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\ground.com:
                Include file for MODFLOW-2000.

            source\MODFLOW-2000\gutsdaf.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1bas6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1bcf6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1chd6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1drn6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1drt1.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1ets1.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1evt6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1fhb1.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1gag5.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1ghb6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1hfb6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1huf2.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1ibs6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1lak3.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1lpf1.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1mnw1.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1mnw2.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1mnwi.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1rch6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1res1.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1riv6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1sfr2.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1str6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1sub1.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1swt1.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\gwf1wel6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\hufutl2.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\hydmod.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\hydmod.inc:
                Include file for MODFLOW-2000.

            source\MODFLOW-2000\lmg1.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\lmt6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\lmt6.inc:
                Include file for MODFLOW-2000.

            source\MODFLOW-2000\memchk.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\mf2k.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\mf2kgmg.c:
                C source code file for MODFLOW-2000.

            source\MODFLOW-2000\mf2kgmg.h:
                C header file for MODFLOW-2000.

            source\MODFLOW-2000\mhc1.f90:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\obs1adv2.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\obs1bas6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\obs1drn6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\obs1drt1.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\obs1ghb6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\obs1riv6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\obs1str6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\openspec.inc:
                Include file for MODFLOW-2000.

            source\MODFLOW-2000\parallel.inc:
                Include file for MODFLOW-2000.

            source\MODFLOW-2000\param.inc:
                Include file for MODFLOW-2000.

            source\MODFLOW-2000\params.inc:
                Include file for MODFLOW-2000.

            source\MODFLOW-2000\parutl1.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\pcg2.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\pes1bas6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\pes1gau1.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\r_vector.c:
                C source code file for MODFLOW-2000.

            source\MODFLOW-2000\r_vector.h:
                C header file for MODFLOW-2000.

            source\MODFLOW-2000\rtedaf.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\sen1bas6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\sen1chd6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\sen1drn6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\sen1drt1.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\sen1ets1.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\sen1evt6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\sen1ghb6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\sen1hfb6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\sen1huf2.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\sen1lpf1.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\sen1rch6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\sen1riv6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\sen1str6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\sen1wel6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\sip5.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\solvers.c:
                C source code file for MODFLOW-2000.

            source\MODFLOW-2000\solvers.h:
                C header file for MODFLOW-2000.

            source\MODFLOW-2000\sor5.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-2000\startdaf.com:
                Include file for MODFLOW-2000.

            source\MODFLOW-2000\utl6.f:
                Fortran source code file for MODFLOW-2000.

            source\MODFLOW-GWT\CCFD.c:
                C source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\CCFD.h:
                C header file for MODFLOW-GWT.

            source\MODFLOW-GWT\MF2KGMG.c:
                C source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\MF2KGMG.h:
                C header file for MODFLOW-GWT.

            source\MODFLOW-GWT\r_vector.c:
                C source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\r_vector.h:
                C header file for MODFLOW-GWT.

            source\MODFLOW-GWT\SOLVERS.c:
                C source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\SOLVERS.h:
                C header file for MODFLOW-GWT.

            source\MODFLOW-GWT\amg1r6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\ctime.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\daf1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\de45.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\glo1bas6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gmg1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\ground.com:
                Include file used in source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gutsdaf.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1bas6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1bcf6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1chd6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1drn6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1drt1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1ets1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1evt6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1fhb1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1gag5.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1ghb6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1hfb6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1huf2.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1ibs6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1lak3.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1lpf1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1mnw1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1mnw2.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1mnwi.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1rch6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1res1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1riv6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1sfr2.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1str6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1sub1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1swt1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwf1wel6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1age6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1bas6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1bdy6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1bflx6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1ccbd1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1cg6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1chd6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1dk6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1dp6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1drn6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1drt1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1dsp6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1ell6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1evt6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1fhb6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1ghb6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1hfb6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1int6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1lak3.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1main6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1mb6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1mbrp1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1mnw1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1mnw2.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1mnwo1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1mov6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1obs6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1pct1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1ptob.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1rch6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1riv6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1sfr2.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1sgm.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1src6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1sstr6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1utl6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1vbal.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1vel6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\gwt1wel6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\hufutl2.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\hydmod.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\hydmod.inc:
                Include file used in source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\lmg1.f:
                Fortran source code file for MODFLOW-GWT. Either this file or
                lmg1OFF.f must be used when compiling MODFLOW-GWT. This file is
                used when the gmg solver is included in MODFLOW-GWT.

            source\MODFLOW-GWT\lmg1OFF.f:
                Fortran source code file for MODFLOW-GWT. Either this file or
                lmg1.f must be used when compiling MODFLOW-GWT. This file is
                used when the gmg solver is NOT included in MODFLOW-GWT.

            source\MODFLOW-GWT\lmt6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\lmt6.inc:
                Include file used in source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\memchk.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\mf2k.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\mhc1.f90:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\move_weight.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\obs1adv2.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\obs1bas6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\obs1drn6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\obs1drt1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\obs1ghb6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\obs1riv6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\obs1str6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\openspec.inc:
                Include file used in source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\parallel.inc:
                Include file used in source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\param.inc:
                Include file used in source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\params.inc:
                Include file used in source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\para-non.f:
                Fortran source code file for MODFLOW-GWT. This file is included
                when compiling a serial version of MODLOW-GWT.

            source\MODFLOW-GWT\parutl1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\pcg2.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\pes1bas6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\pes1gau1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\ppc_v3.0.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\ptwt.inc:
                Include file used in source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\ptwt1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\rtedaf.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\sen1bas6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\sen1chd6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\sen1drn6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\sen1drt1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\sen1ets1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\sen1evt6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\sen1ghb6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\sen1hfb6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\sen1huf2.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\sen1lpf1.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\sen1rch6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\sen1riv6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\sen1str6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\sen1wel6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\sip5.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\sor5.f:
                Fortran source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\startdaf.com:
                Include file used in source code file for MODFLOW-GWT.

            source\MODFLOW-GWT\utl6.f:
                Fortran source code file for MODFLOW-GWT.

            source\MT3DMS\automake.fig:
                Configuration file for compiling MT3DMS

            source\MT3DMS\filespec.inc:
                Include file for MT3DMS.

            source\MT3DMS\mt_adv5.for:
                Fortran source code file for MT3DMS.

            source\MT3DMS\mt_btn5.for:
                Fortran source code file for MT3DMS.

            source\MT3DMS\mt_dsp5.for:
                Fortran source code file for MT3DMS.

            source\MT3DMS\mt_fmi5.for:
                Fortran source code file for MT3DMS.

            source\MT3DMS\mt_gcg5.for:
                Fortran source code file for MT3DMS.

            source\MT3DMS\mt_hss5.for:
                Fortran source code file for MT3DMS.

            source\MT3DMS\mt_rct5.for:
                Fortran source code file for MT3DMS.

            source\MT3DMS\mt_ssm5.for:
                Fortran source code file for MT3DMS.

            source\MT3DMS\mt_tob5.for:
                Fortran source code file for MT3DMS.

            source\MT3DMS\mt_utl5.for:
                Fortran source code file for MT3DMS.

            source\MT3DMS\mt3dms5.for:
                Fortran source code file for MT3DMS.

            source\Point3\dgdiss.f:
                Fortran Source code file for Point3 as documented in Wexler
                (1992).

            source\Point3\dimens.inc:
                Include file for Point3 as documented in Wexler (1992).

            source\Point3\exerfc.f:
                Fortran Source code file for Point3.

            source\Point3\glqpts.f:
                Fortran Source code file for Point3 as documented in Wexler
                (1992).

            source\Point3\lenchr.f:
                Fortran Source code file for Point3 as documented in Wexler
                (1992).

            source\Point3\locsubs.intel.f:
                Fortran Source code file for Point3 as documented in Wexler
                (1992). This file has been modified to include compiler-specific
                functions for date and time for the Intel fortran compiler.

            source\Point3\ofile.f:
                Fortran Source code file for Point3 as documented in Wexler
                (1992).

            source\Point3\plot2d.f:
                Fortran Source code file for Point3 as documented in Wexler
                (1992).

            source\Point3\pltdat.inc:
                Include file for Point3 as documented in Wexler (1992).

            source\Point3\point3.f:
                Fortran Source code file for Point3 as documented in Wexler
                (1992).

            source\Point3\rdplot.f:
                Fortran Source code file for Point3 as documented in Wexler
                (1992).

            source\Point3\subs3.f:
                Fortran Source code file for Point3 as documented in Wexler
                (1992).

            source\Point3\title.f:
                Fortran Source code file for Point3 as documented in Wexler
                (1992).

            source\Point3\units.inc:
                Include file for Point3 as documented in Wexler (1992).

        webrelease\BrowseGraphic.png:
            Image of the model used for benchmarking MODFLOW-GWT for comparison
            with MT3DMS.
            The test problem represents an analog based on, and greatly
            simplified from, the groundwater contamination problem at the Rocky
            Mountain Arsenal, Colorado (see Konikow, 1977).

        webrelease\Metadata.xml:
            FGDC metadata for this archive

        ancillary\rmamocex.mmb:
            Argus ONE project file for rma_moc model.

        ancillary\rmamocel.mmb:
            Argus ONE project file for rmamocel model.

        ancillary\rmamocim.mmb:
            Argus ONE project file for rmamocim model.

        ancillary\rmamocwt_dist.mmb:
            Argus ONE project file for rmamocwt_dist model.

        ancillary\rmamocwti.mmb:
            Argus ONE project file for rmamocwti model.

        ancillary\rmamtfd.mmb:
            Argus ONE project file for fd and fd_MT3D models.

        ancillary\rmamthmoc.mmb:
            Argus ONE project file for hmoc and hmoc_MT3D models.

        ancillary\rmamtmmoc.mmb:
            Argus ONE project file for mmoc and mmoc_MT3D models.

        ancillary\rmamtmoc.mmb:
            Argus ONE project file for moc and moc_MT3D models.

        ancillary\rmamttvd.mmb:
            Argus ONE project file for tvd and tvd_MT3D models.

        ancillary\finite.mmb:
            Argus ONE project file for Finite.1D model.

        ancillary\finiteAlpha1.mmb:
            Argus ONE project file

        ancillary\DiracUniDirectional.mmb:
            Argus ONE project file for the DiracUnidirectional model

        ancillary\DiracAngle.mmb:
            Argus ONE project file for the DiracAngle model


        ancillary\Point3\
            Description:
            -----------
            Point3 contains files related to the Point3 MODFLOW-GWT model.

            Files:
            -----
            ancillary\Point3\Point3.mmb:
                Argus ONE project file for the Point3 model.

            ancillary\Point3\Point3.mv:
                Model Viewer file for the Point3 model. To use, Point3.mv must
                be in the same directory as the input and output files for the
                model. ModelViewer (Hsieh and Winston, 2002) can be used to open
                Point3.mv.

