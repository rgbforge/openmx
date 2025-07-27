
cmake_minimum_required(VERSION 3.12)
project(OpenMX C Fortran)



set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} $ENV{C_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} $ENV{Fortran_FLAGS}")


set(LIBELPADIR "${CMAKE_CURRENT_SOURCE_DIR}/elpa-2018.05.001")
set(FORT_MOD_DIR "${CMAKE_CURRENT_BINARY_DIR}/fortran_modules")


include_directories("$ENV{MKLROOT}/include")
include_directories("$ENV{MKLROOT}/include/fftw")
include_directories(${LIBELPADIR})


set(OPENMX_C_SOURCES
    openmx.c openmx_common.c Input_std.c Inputtools.c init.c LU_inverse.c ReLU_inverse.c
    truncation.c readfile.c FT_PAO.c FT_NLP.c FT_ProExpn_VNA.c FT_VNA.c FT_ProductPAO.c
    Hamiltonian_Cluster.c Hamiltonian_Cluster_Hs.c Hamiltonian_Cluster_NC_Hs2.c Hamiltonian_Band_NC_Hs2.c
    Overlap_Cluster_NC_Ss2.c Overlap_Band_NC_Ss2.c Overlap_Cluster.c Overlap_Cluster_Ss.c
    Set_ContMat_Cluster_LNO.c Hamiltonian_Band.c Matrix_Band_LNO.c Overlap_Band.c Hamiltonian_Cluster_NC.c
    Hamiltonian_Band_NC.c Hamiltonian_Cluster_SO.c Get_OneD_HS_Col.c SetPara_DFT.c XC_Ceperly_Alder.c XC_CA_LSDA.c
    XC_PW92C.c XC_PBE.c XC_EX.c DFT.c Mixing_DM.c Mixing_H.c Mixing_V.c Force.c Stress.c Poisson.c Poisson_ESM.c
    Cluster_DFT_Col.c Cluster_DFT_NonCol.c Cluster_DFT_Dosout.c Cluster_DFT_ON2.c Cluster_DFT_LNO.c
    Band_DFT_Col.c Band_DFT_NonCol.c Band_DFT_NonCol_GB.c Band_DFT_kpath.c Band_DFT_kpath_LNO.c Band_DFT_MO.c
    Unfolding_Bands.c Band_DFT_Dosout.c Set_Density_Grid.c Set_Orbitals_Grid.c Set_Aden_Grid.c
    Gauss_Legendre.c zero_cfrac.c xyz2spherical.c AngularF.c RadialF.c Dr_RadialF.c PhiF.c
    VNAF.c Dr_VNAF.c VH_AtomF.c Dr_VH_AtomF.c RF_BesselF.c QuickSort.c Nonlocal_RadialF.c
    KumoF.c Dr_KumoF.c Mulliken_Charge.c Occupation_Number_LDA_U.c Eff_Hub_Pot.c Coulomb_Interaction.c
    EulerAngle_Spin.c Smoothing_Func.c Orbital_Moment.c Pot_NeutralAtom.c Simple_Mixing_DM.c
    DIIS_Mixing_DM.c ADIIS_Mixing_DM.c GR_Pulay_DM.c Kerker_Mixing_Rhok.c DIIS_Mixing_Rhok.c
    Total_Energy.c Contract_Hamiltonian.c Contract_iHNL.c Cont_Matrix0.c Cont_Matrix1.c Cont_Matrix2.c
    Cont_Matrix3.c Cont_Matrix4.c Opt_Contraction.c Initial_CntCoes.c Initial_CntCoes2.c
    Set_XC_Grid.c Set_XC_NL1_Grid.c Get_Orbitals.c Get_dOrbitals.c Get_Cnt_Orbitals.c
    Get_Cnt_dOrbitals.c Gaunt.c Find_CGrids.c MD_pac.c RestartFileDFT.c Output_CompTime.c Merge_LogFile.c Make_FracCoord.c
    Make_InputFile_with_FinalCoord.c Output_Energy_Decomposition.c Divide_Conquer.c Divide_Conquer_LNO.c Krylov.c
    Divide_Conquer_Dosout.c EGAC_DFT.c LNO.c Eigen_lapack.c Eigen_lapack2.c Eigen_lapack3.c EigenBand_lapack.c
    Eigen_PReHH.c BroadCast_ReMatrix.c Eigen_PHH.c BroadCast_ComplexMatrix.c lapack_dstedc1.c lapack_dstedc2.c lapack_dstedc3.c
    lapack_dstegr1.c lapack_dstegr2.c lapack_dstegr3.c lapack_dstevx1.c lapack_dstevx2.c lapack_dstevx3.c
    lapack_dstevx4.c lapack_dstevx5.c lapack_dsteqr1.c Nonlocal_Basis.c Set_OLP_Kin.c Set_Nonlocal.c Set_ProExpn_VNA.c
    Set_CoreHoleMatrix.c Set_OLP_p.c Set_Hamiltonian.c Set_Vpot.c Voronoi_Charge.c Voronoi_Orbital_Moment.c Fuzzy_Weight.c
    dampingF.c deri_dampingF.c Spherical_Bessel.c iterout.c iterout_md.c Allocate_Arrays.c Free_Arrays.c
    Init_List_YOUSO.c outputfile1.c malloc_multidimarray.c PrintMemory.c PrintMemory_Fix.c
    dtime.c OutData.c OutData_Binary.c init_alloc_first.c File_CntCoes.c SCF2File.c mimic_sse.c Make_Comm_Worlds.c
    Set_Allocate_Atom2CPU.c Cutoff.c Generating_MP_Special_Kpt.c Maketest.c Runtest.c Memory_Leak_test.c
    Force_test.c Stress_test.c Show_DFT_DATA.c Generate_Wannier.c TRAN_Allocate.c TRAN_DFT.c TRAN_DFT_Dosout.c TRAN_Apply_Bias2e.c
    TRAN_Deallocate_Electrode_Grid.c TRAN_Deallocate_RestartFile.c TRAN_RestartFile.c TRAN_Calc_CentGreen.c TRAN_Input_std.c
    TRAN_Set_CentOverlap.c TRAN_Calc_CentGreenLesser.c TRAN_Input_std_Atoms.c TRAN_Set_Electrode_Grid.c
    TRAN_Calc_GridBound.c TRAN_Set_IntegPath.c TRAN_Output_HKS.c TRAN_Set_MP.c TRAN_Calc_SelfEnergy.c TRAN_Output_Trans_HS.c
    TRAN_Calc_Hopping_G.c TRAN_Calc_SurfGreen.c TRAN_Set_SurfOverlap.c TRAN_Add_Density_Lead.c TRAN_Add_ADensity_Lead.c TRAN_Set_Value.c
    TRAN_Poisson.c TRAN_adjust_Ngrid.c TRAN_Print.c TRAN_Print_Grid.c Lapack_LU_inverse.c TRAN_Distribute_Node.c TRAN_Output_HKS_Write_Grid.c
    TRAN_Credit.c TRAN_Check_Region_Lead.c TRAN_Check_Region.c TRAN_Check_Input.c DFTDvdW_init.c DFTD3vdW_init.c neb.c neb_run.c neb_check.c
    TRAN_Allocate_NC.c TRAN_DFT_NC.c TRAN_Set_CentOverlap_NC.c TRAN_Set_SurfOverlap_NC.c
    TRAN_Calc_OneTransmission.c TRAN_Main_Analysis.c TRAN_Main_Analysis_NC.c MTRAN_EigenChannel.c TRAN_Channel_Functions.c TRAN_Channel_Output.c
    TRAN_Calc_CurrentDensity.c TRAN_CDen_Main.c NBO_Cluster.c NBO_Krylov.c Population_Analysis_Wannier.c Population_Analysis_Wannier2.c
    NabraMatrixElements.c Set_dOrbitals_Grid.c Calc_optical.c Band_DFT_NonCol_Optical.c Cluster_DFT_Optical.c
    Band_DFT_Col_Optical_ScaLAPACK.c Cluster_DFT_Optical_ScaLAPACK.c
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
set(LOCAL_FORTRAN_SOURCES elpa1.f90 solve_evp_real.f90 solve_evp_complex.f90)


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
target_link_libraries(openmx PRIVATE $ENV{LINKER_LIBS})
install(TARGETS openmx RUNTIME DESTINATION ../work)


set(GCC_UTILS bandgnu13 bin2txt cube2xsf intensity_map md2axsf gcube2oned)
foreach(util ${GCC_UTILS})
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${util}
        COMMAND gcc ${CMAKE_CURRENT_SOURCE_DIR}/${util}.c -lm -o ${CMAKE_CURRENT_BINARY_DIR}/${util}
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${util}.c
        COMMENT "Building ${util} with gcc"
    )
    add_custom_target(build_${util} ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${util})
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${util} DESTINATION ../work PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)
endforeach()


function(add_mpi_util name)
    add_executable(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE $ENV{LINKER_LIBS})
    install(TARGETS ${name} RUNTIME DESTINATION ../work)
endfunction()

add_mpi_util(DosMain DosMain.c Inputtools.c malloc_multidimarray.c Tetrahedron_Blochl.c)
add_mpi_util(jx jx.c jx_quicksort.c jx_LNO.c jx_config.c jx_tools.c jx_cluster.c jx_band_psum.c jx_band_indiv.c Inputtools.c read_scfout.c)
add_mpi_util(analysis_example analysis_example.c read_scfout.c)
add_mpi_util(esp esp.c Inputtools.c)
add_mpi_util(polB polB.c read_scfout.c)
add_mpi_util(calB calB.c read_scfout.c)
add_mpi_util(Z2FH Z2FH.c read_scfout.c)
add_mpi_util(tp tp.c read_scfout.c)
add_mpi_util(example_mpi_spawn example_mpi_spawn.c)
