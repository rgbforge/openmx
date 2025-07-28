cmake_minimum_required(VERSION 3.12)
project(OpenMX C Fortran)
set(CMAKE_C_COMPILER "mpiicx")
set(CMAKE_Fortran_COMPILER "mpiifx")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -xHOST -fiopenmp -fcommon -Wno-error=implicit-function-declaration -g")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -xHOST -fiopenmp -g")

set(SOURCE_SUBDIR "${CMAKE_SOURCE_DIR}/source")
set(LIBELPADIR "${SOURCE_SUBDIR}/elpa-2018.05.001")
set(FORT_MOD_DIR "${CMAKE_BINARY_DIR}/fortran_modules")
set(INSTALL_DIR "${CMAKE_SOURCE_DIR}/work")

link_directories("$ENV{MKLROOT}/lib/intel64")

set(MKL_LIBS
    mkl_scalapack_lp64
    mkl_intel_lp64
    mkl_intel_thread
    mkl_core
    mkl_blacs_intelmpi_lp64
)
set(SYSTEM_LIBS
    iomp5
    pthread
    m
    dl
    ifcore
)

include_directories("$ENV{MKLROOT}/include")
include_directories("$ENV{MKLROOT}/include/fftw")
include_directories(${LIBELPADIR})

set(OPENMX_C_SOURCES
    source/openmx.c source/openmx_common.c source/Input_std.c source/Inputtools.c source/init.c source/LU_inverse.c source/ReLU_inverse.c
    source/truncation.c source/readfile.c source/FT_PAO.c source/FT_NLP.c source/FT_ProExpn_VNA.c source/FT_VNA.c source/FT_ProductPAO.c
    source/Hamiltonian_Cluster.c source/Hamiltonian_Cluster_Hs.c source/Hamiltonian_Cluster_NC_Hs2.c source/Hamiltonian_Band_NC_Hs2.c
    source/Overlap_Cluster_NC_Ss2.c source/Overlap_Band_NC_Ss2.c source/Overlap_Cluster.c source/Overlap_Cluster_Ss.c
    source/Set_ContMat_Cluster_LNO.c source/Hamiltonian_Band.c source/Matrix_Band_LNO.c source/Overlap_Band.c source/Hamiltonian_Cluster_NC.c
    source/Hamiltonian_Band_NC.c source/Hamiltonian_Cluster_SO.c source/Get_OneD_HS_Col.c source/SetPara_DFT.c source/XC_Ceperly_Alder.c source/XC_CA_LSDA.c
    source/XC_PW92C.c source/XC_PBE.c source/XC_EX.c source/DFT.c source/Mixing_DM.c source/Mixing_H.c source/Mixing_V.c source/Force.c source/Stress.c source/Poisson.c source/Poisson_ESM.c
    source/Cluster_DFT_Col.c source/Cluster_DFT_NonCol.c source/Cluster_DFT_Dosout.c source/Cluster_DFT_ON2.c source/Cluster_DFT_LNO.c
    source/Band_DFT_Col.c source/Band_DFT_NonCol.c source/Band_DFT_NonCol_GB.c source/Band_DFT_kpath.c source/Band_DFT_kpath_LNO.c source/Band_DFT_MO.c
    source/Unfolding_Bands.c source/Band_DFT_Dosout.c source/Set_Density_Grid.c source/Set_Orbitals_Grid.c source/Set_Aden_Grid.c
    source/Gauss_Legendre.c source/zero_cfrac.c source/xyz2spherical.c source/AngularF.c source/RadialF.c source/Dr_RadialF.c source/PhiF.c
    source/VNAF.c source/Dr_VNAF.c source/VH_AtomF.c source/Dr_VH_AtomF.c source/RF_BesselF.c source/QuickSort.c source/Nonlocal_RadialF.c
    source/KumoF.c source/Dr_KumoF.c source/Mulliken_Charge.c source/Occupation_Number_LDA_U.c source/Eff_Hub_Pot.c source/Coulomb_Interaction.c
    source/EulerAngle_Spin.c source/Smoothing_Func.c source/Orbital_Moment.c source/Pot_NeutralAtom.c source/Simple_Mixing_DM.c
    source/DIIS_Mixing_DM.c source/ADIIS_Mixing_DM.c source/GR_Pulay_DM.c source/Kerker_Mixing_Rhok.c source/DIIS_Mixing_Rhok.c
    source/Total_Energy.c source/Contract_Hamiltonian.c source/Contract_iHNL.c source/Cont_Matrix0.c source/Cont_Matrix1.c source/Cont_Matrix2.c
    source/Cont_Matrix3.c source/Cont_Matrix4.c source/Opt_Contraction.c source/Initial_CntCoes.c source/Initial_CntCoes2.c
    source/Set_XC_Grid.c source/Set_XC_NL1_Grid.c source/Get_Orbitals.c source/Get_dOrbitals.c source/Get_Cnt_Orbitals.c
    source/Get_Cnt_dOrbitals.c source/Gaunt.c source/Find_CGrids.c source/MD_pac.c source/RestartFileDFT.c source/Output_CompTime.c source/Merge_LogFile.c source/Make_FracCoord.c
    source/Make_InputFile_with_FinalCoord.c source/Output_Energy_Decomposition.c source/Divide_Conquer.c source/Divide_Conquer_LNO.c source/Krylov.c
    source/Divide_Conquer_Dosout.c source/EGAC_DFT.c source/LNO.c source/Eigen_lapack.c source/Eigen_lapack2.c source/Eigen_lapack3.c source/EigenBand_lapack.c
    source/Eigen_PReHH.c source/BroadCast_ReMatrix.c source/Eigen_PHH.c source/BroadCast_ComplexMatrix.c source/lapack_dstedc1.c source/lapack_dstedc2.c source/lapack_dstedc3.c
    source/lapack_dstegr1.c source/lapack_dstegr2.c source/lapack_dstegr3.c source/lapack_dstevx1.c source/lapack_dstevx2.c source/lapack_dstevx3.c
    source/lapack_dstevx4.c source/lapack_dstevx5.c source/lapack_dsteqr1.c source/Nonlocal_Basis.c source/Set_OLP_Kin.c source/Set_Nonlocal.c source/Set_ProExpn_VNA.c
    source/Set_CoreHoleMatrix.c source/Set_OLP_p.c source/Set_Hamiltonian.c source/Set_Vpot.c source/Voronoi_Charge.c source/Voronoi_Orbital_Moment.c source/Fuzzy_Weight.c
    source/dampingF.c source/deri_dampingF.c source/Spherical_Bessel.c source/iterout.c source/iterout_md.c source/Allocate_Arrays.c source/Free_Arrays.c
    source/Init_List_YOUSO.c source/outputfile1.c source/malloc_multidimarray.c source/PrintMemory.c source/PrintMemory_Fix.c
    source/dtime.c source/OutData.c source/OutData_Binary.c source/init_alloc_first.c source/File_CntCoes.c source/SCF2File.c source/mimic_sse.c source/Make_Comm_Worlds.c
    source/Set_Allocate_Atom2CPU.c source/Cutoff.c source/Generating_MP_Special_Kpt.c source/Maketest.c source/Runtest.c source/Memory_Leak_test.c
    source/Force_test.c source/Stress_test.c source/Show_DFT_DATA.c source/Generate_Wannier.c source/TRAN_Allocate.c source/TRAN_DFT.c source/TRAN_DFT_Dosout.c source/TRAN_Apply_Bias2e.c
    source/TRAN_Deallocate_Electrode_Grid.c source/TRAN_Deallocate_RestartFile.c source/TRAN_RestartFile.c source/TRAN_Calc_CentGreen.c source/TRAN_Input_std.c
    source/TRAN_Set_CentOverlap.c source/TRAN_Calc_CentGreenLesser.c source/TRAN_Input_std_Atoms.c source/TRAN_Set_Electrode_Grid.c
    source/TRAN_Calc_GridBound.c source/TRAN_Set_IntegPath.c source/TRAN_Output_HKS.c source/TRAN_Set_MP.c source/TRAN_Calc_SelfEnergy.c source/TRAN_Output_Trans_HS.c
    source/TRAN_Calc_Hopping_G.c source/TRAN_Calc_SurfGreen.c source/TRAN_Set_SurfOverlap.c source/TRAN_Add_Density_Lead.c source/TRAN_Add_ADensity_Lead.c source/TRAN_Set_Value.c
    source/TRAN_Poisson.c source/TRAN_adjust_Ngrid.c source/TRAN_Print.c source/TRAN_Print_Grid.c source/Lapack_LU_inverse.c source/TRAN_Distribute_Node.c source/TRAN_Output_HKS_Write_Grid.c
    source/TRAN_Credit.c source/TRAN_Check_Region_Lead.c source/TRAN_Check_Region.c source/TRAN_Check_Input.c source/DFTDvdW_init.c source/DFTD3vdW_init.c source/neb.c source/neb_run.c source/neb_check.c
    source/TRAN_Allocate_NC.c source/TRAN_DFT_NC.c source/TRAN_Set_CentOverlap_NC.c source/TRAN_Set_SurfOverlap_NC.c
    source/TRAN_Calc_OneTransmission.c source/TRAN_Main_Analysis.c source/TRAN_Main_Analysis_NC.c source/MTRAN_EigenChannel.c source/TRAN_Channel_Functions.c source/TRAN_Channel_Output.c
    source/TRAN_Calc_CurrentDensity.c source/TRAN_CDen_Main.c source/NBO_Cluster.c source/NBO_Krylov.c source/Population_Analysis_Wannier.c source/Population_Analysis_Wannier2.c
    source/NabraMatrixElements.c source/Set_dOrbitals_Grid.c source/Calc_optical.c source/Band_DFT_NonCol_Optical.c source/Cluster_DFT_Optical.c
    source/Band_DFT_Col_Optical_ScaLAPACK.c source/Cluster_DFT_Optical_ScaLAPACK.c
)
set(ELPA_FORTRAN_SOURCES
    ${LIBELPADIR}/mod_precision.F90 ${LIBELPADIR}/elpa_utilities.F90 ${LIBELPADIR}/elpa1_compute_real.F90
    ${LIBELPADIR}/elpa1_compute_complex.F90 ${LIBELPADIR}/aligned_mem.F90 ${LIBELPADIR}/elpa2_determine_workload.F90
    ${LIBELPADIR}/mod_redist_band_real.F90 ${LIBELPADIR}/mod_redist_band_complex.F90 ${LIBELPADIR}/mod_pack_unpack_cpu_real.F90
    ${LIBELPADIR}/mod_pack_unpack_cpu_complex.F90 ${LIBELPADIR}/real.F90 ${LIBELPADIR}/complex.F90
    ${LIBELPADIR}/mod_single_hh_trafo_real.F90 ${LIBELPADIR}/mod_compute_hh_trafo_real.F90 ${LIBELPADIR}/mod_compute_hh_trafo_complex.F90
    ${LIBELPADIR}/elpa2_compute_real.F90 ${LIBELPADIR}/elpa2_compute_complex.F90
    ${LIBELPADIR}/elpa_solve_evp_real_2stage_double_impl.F90 ${LIBELPADIR}/elpa_solve_evp_complex_2stage_double_impl.F90
)
set(LOCAL_FORTRAN_SOURCES source/elpa1.f90 source/solve_evp_real.f90 source/solve_evp_complex.f90)

add_library(elpa_lib OBJECT ${ELPA_FORTRAN_SOURCES})
set_target_properties(elpa_lib PROPERTIES Fortran_MODULE_DIRECTORY ${FORT_MOD_DIR})

add_library(local_fortran_lib OBJECT ${LOCAL_FORTRAN_SOURCES})
target_include_directories(local_fortran_lib PRIVATE ${FORT_MOD_DIR})
add_dependencies(local_fortran_lib elpa_lib)

add_executable(openmx
    ${OPENMX_C_SOURCES}
    $<TARGET_OBJECTS:elpa_lib>
    $<TARGET_OBJECTS:local_fortran_lib>
)
set_target_properties(openmx PROPERTIES LINKER_LANGUAGE C)

target_link_libraries(openmx PRIVATE ${MKL_LIBS} ${SYSTEM_LIBS})
install(TARGETS openmx RUNTIME DESTINATION ${INSTALL_DIR})

set(GCC_UTILS bandgnu13 bin2txt cube2xsf intensity_map md2axsf gcube2oned)
foreach(util ${GCC_UTILS})
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${util}
        COMMAND gcc ${SOURCE_SUBDIR}/${util}.c -lm -o ${CMAKE_CURRENT_BINARY_DIR}/${util}
        DEPENDS ${SOURCE_SUBDIR}/${util}.c
    )
    add_custom_target(build_${util} ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${util})
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${util} DESTINATION ${INSTALL_DIR} PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)
endforeach()

function(add_mpi_util name)
    add_executable(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE ${MKL_LIBS} ${SYSTEM_LIBS})
    install(TARGETS ${name} RUNTIME DESTINATION ${INSTALL_DIR})
endfunction()

add_mpi_util(DosMain source/DosMain.c source/Inputtools.c source/malloc_multidimarray.c source/Tetrahedron_Blochl.c)
add_mpi_util(jx source/jx.c source/jx_quicksort.c source/jx_LNO.c source/jx_config.c source/jx_tools.c source/jx_cluster.c source/jx_band_psum.c source/jx_band_indiv.c source/Inputtools.c source/read_scfout.c)
add_mpi_util(analysis_example source/analysis_example.c source/read_scfout.c)
add_mpi_util(esp source/esp.c source/Inputtools.c)
add_mpi_util(polB source/polB.c source/read_scfout.c)
add_mpi_util(calB source/calB.c source/read_scfout.c)
add_mpi_util(Z2FH source/Z2FH.c source/read_scfout.c)
add_mpi_util(tp source/tp.c source/read_scfout.c)
add_mpi_util(example_mpi_spawn source/example_mpi_spawn.c)
