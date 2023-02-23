#!/bin/bash
#################################################################
#             This is an automatically generated patch          #
#                for the Cc4s interface to VASP.6.3.1 .         #
#    For questions and problems visit https://github.com/cc4s   #
#################################################################
echo "Trying to patch VASP version 6.3.1 to include Cc4s interface. " 
 
#######################
# diff of   chi.F
#######################
 
 
cat > .tmp_vaspcc4s_patch<<"EOF"
39c39,40
<        IF (LCHI .OR. LMP2 .OR. LMP2KPAR .OR. LMP2NO .OR. LFCIDUMP .OR. LRPAX .OR. LCCSD .OR. LBRACKETST) THEN
---
>        IF (LCHI .OR. LMP2 .OR. LMP2KPAR .OR. LMP2NO .OR. LFCIDUMP .OR. LRPAX .OR. LCCSD .OR. LBRACKETST &
>           .OR. LCC4S) THEN
59c60
<       & .OR. LCCSD .OR. LBRACKETST) .AND. ENCUTGWSOFT_OLD ==-2 ) THEN
---
>       & .OR. LCCSD .OR. LBRACKETST .OR. LCC4S) .AND. ENCUTGWSOFT_OLD ==-2 ) THEN
87c88,89
<     IF (IU6>=0 .AND. (LCHI.OR. LTIME_EVOLUTION .OR.LMP2.OR.LMP2KPAR.OR.LMP2NO .OR. LFCIDUMP  .OR. LRPAX .OR.  LCCSD .OR. LBRACKETST)) THEN
---
>     IF (IU6>=0 .AND. (LCHI.OR. LTIME_EVOLUTION .OR.LMP2.OR.LMP2KPAR.OR.LMP2NO .OR. LFCIDUMP  .OR. LRPAX &
>        .OR.  LCCSD .OR. LBRACKETST .OR. LCC4S)) THEN
2939c2941,2945
<           T_INFO, DYN, IO,KPOINTS,SYMM,GRID,INFO,AMIX,BMIX)
---
>           T_INFO, DYN, IO,KPOINTS,SYMM,GRID,INFO,AMIX,BMIX, &
>           HAMILTONIAN,LMDIM,CDIJ,CQIJ,SV,E,&
>           GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C,&
>           CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
>           CRHODE, N_MIX_PAW, RHOLM, CHDEN)
2943a2950
>     USE cc4s
2946c2953
< !    USE ccsd
---
>     USE ccsd
2948c2955
< !    USE bracketst
---
>     USE bracketst
2963a2971
>     USE hamil_struct_def
2980a2989,3009
>     TYPE (ham_handle)  HAMILTONIAN
>     INTEGER  LMDIM
>     OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
>     OVERLAP  CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
>     RGRID       SV(DIMREAL(WDES%GRID%MPLWV),WDES%NCDIJ)
>     TYPE (energy)      E
>     TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
>     TYPE (grid_3d)     GRID_SOFT  ! grid for soft chargedensity
>     TYPE (grid_3d)     GRIDUS     ! temporary grid in us.F
>     TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
>     TYPE (transit)     SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
>     COMPLEX(q) CHTOT(GRIDC%MPLWV,WDES%NCDIJ) ! charge density
>     COMPLEX(q) CHTOTL(GRIDC%MPLWV,WDES%NCDIJ)! old charge-density
>     RGRID      DENCOR(GRIDC%RL%NP)
>     COMPLEX(q)  CVTOT(GRIDC%MPLWV,WDES%NCDIJ) ! local potential
>     COMPLEX(q)  CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
>     INTEGER     IRDMAX
>     OVERLAP  CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
>     INTEGER N_MIX_PAW
>     REAL(q)  RHOLM(N_MIX_PAW,WDES%NCDIJ),RHOLM_LAST(N_MIX_PAW,WDES%NCDIJ)
>     COMPLEX(q)  CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
2984a3014
> 
2987c3017
<     CALL READ_CDER_BETWEEN_STATES(WDES, IO%IU0, 55)
---
>     ! CALL READ_CDER_BETWEEN_STATES(WDES, IO%IU0, 55)
3043c3073,3075
<     IF (LMP2) THEN
---
>     IF (LMP2KPAR) THEN
>        CALL CALCULATE_MP2_KPAR( &
>           P,WDES,W,LATT_CUR,T_INFO,IO,KPOINTS,WGW, ENCUTGW, ENCUTGWSOFT, LMAXMP2)
3044a3077
>     ELSEIF (LMP2) THEN
3047,3050d3079
< #endif
<     ELSEIF (LMP2KPAR) THEN
<        CALL CALCULATE_MP2_KPAR( &
<           P,WDES,W,LATT_CUR,T_INFO,IO,KPOINTS,WGW, ENCUTGW, ENCUTGWSOFT, LMAXMP2)
3052d3080
< #ifdef scaLAPACK
3055c3083,3089
< #endif
---
>     ELSEIF (LCC4S) THEN
>        CALL CC4S_INTERFACE(P,WDES,W,LATT_CUR,LATT_INI,T_INFO,IO,KPOINTS,WGW, ENCUTGW, ENCUTGWSOFT, LMAXMP2, INFO, &
>                      HAMILTONIAN,SYMM,GRID,NONLR_S,NONL_S,LMDIM,CDIJ,CQIJ,SV,E,&
>                      GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C,&
>                      CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
>                      CRHODE, N_MIX_PAW, RHOLM, CHDEN)
>  
3057d3090
< #ifdef scaLAPACK
3060,3063c3093,3094
< #endif
< #ifdef vasp6
< !    ELSEIF (LCCSD) THEN
< ! #ifdef gammareal
---
>     ELSEIF (LCCSD) THEN
> #ifdef gammareal
3066c3097
< ! #else
---
> #else
3069,3072c3100,3105
< ! #endif
< !    ELSEIF (LBRACKETST) THEN
< !       CALL CALCULATE_BRACKETST(P,WDES,W,LATT_CUR,T_INFO,IO,KPOINTS,WGW, &
< !         & LMAXMP2)
---
>        CALL CALCULATE_CCSD(SYMM,P,WDES,W,LATT_CUR,T_INFO,INFO,IO,KPOINTS,WGW,ENCUTGW,  & 
>          & ENCUTGWSOFT, LMAXMP2,AMIX,BMIX,NONLR_S,NONL_S,GRID,LATT_INI)
> #endif
>     ELSEIF (LBRACKETST) THEN
>        CALL CALCULATE_BRACKETST(P,WDES,W,LATT_CUR,T_INFO,IO,KPOINTS,WGW, &
>          & LMAXMP2)
EOF

md5sum_orig=5aa367bef1a214f55766b39416686e7e
md5sum_new=`md5sum chi.F  | awk '{print $1}' `
if [ "$md5sum_orig" = "$md5sum_new" ]; then
    echo "Patching file  chi.F " 
    patch chi.F .tmp_vaspcc4s_patch 
else
   echo "You dont have the correct version of file  chi.F  to apply patch. " 
fi
 
#######################
# diff of   chi_glb.F
#######################
 
 
cat > .tmp_vaspcc4s_patch<<"EOF"
427a428
>    LOGICAL :: LCC4S=.FALSE.
1079a1081,1086
>     ELSE IF (TEXT(1:N)=='ccsd') THEN
>        LCHI=.TRUE.
>        LCCSD=.TRUE.
>     ELSE IF (TEXT(1:N)=='cc4s') THEN
>        LCHI=.TRUE.
>        LCC4S=.TRUE.
1168a1176,1177
>     ELSE IF (LCC4S) THEN
>        ALGO_FROM_GW='CC4S'
1745c1754
<     IF (LCHI.OR.LMP2 .OR. LMP2KPAR .OR.L2E4W .OR. LCRPA) THEN
---
>     IF (LCHI.OR.LMP2 .OR. LMP2NO .OR. LMP2KPAR .OR.L2E4W .OR. LCRPA .OR. LCCSD .OR. LCC4S ) THEN
EOF

md5sum_orig=999a3a5d2afcf1df8d2a49771b2c5565
md5sum_new=`md5sum chi_glb.F  | awk '{print $1}' `
if [ "$md5sum_orig" = "$md5sum_new" ]; then
    echo "Patching file  chi_glb.F " 
    patch chi_glb.F .tmp_vaspcc4s_patch 
else
   echo "You dont have the correct version of file  chi_glb.F  to apply patch. " 
fi
 
#######################
# diff of   main.F
#######################
 
 
cat > .tmp_vaspcc4s_patch<<"EOF"
2862a2863,2868
>          
>          !we don't want to reorthonormalize previously obtained orbitals for wavefunction based methods
>          !this causes problems when we work with shifted k-meshes!
>          IF (.NOT.(LMP2NO.OR.LMP2.OR.LCCSD.OR.LBRACKETST.OR.LCC4S)) THEN
>             CALL ORTHCH(WDES,W, INFO%LOVERL, LMDIM,CQIJ)
>          ENDIF
2864,2865d2869
<          CALL ORTHCH(WDES,W, INFO%LOVERL, LMDIM,CQIJ)
< 
2992c2996,2998
<          IF (LMP2 .OR. LMP2KPAR .OR. LMP2NO .OR. LFCIDUMP  .OR. LRPAX .OR.  LCCSD .OR. LBRACKETST) THEN
---
>          IF (LMP2 .OR. LMP2KPAR .OR. LMP2NO .OR. LFCIDUMP  .OR. LCC4S .OR. &
>             & LRPAX .OR.  LCCSD .OR. LBRACKETST) THEN
> 
2994c3000,3005
<              & T_INFO, DYN, IO,KPOINTS,SYMM,GRID,INFO,MIX%AMIX,MIX%BMIX)
---
>                T_INFO, DYN, IO,KPOINTS,SYMM,GRID,INFO,MIX%AMIX,MIX%BMIX, &
>                HAMILTONIAN,LMDIM,CDIJ,CQIJ,SV,E,&
>                GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C,&
>                CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
>                CRHODE, N_MIX_PAW, RHOLM, CHDEN)
> 
3003,3004c3014,3015
<          IF (LMP2 .OR. LMP2KPAR .OR. LMP2NO .OR. LFCIDUMP  .OR. LRPAX .OR.  LCCSD .OR. &
<             & LBRACKETST .OR. FOURORBIT==1 ) THEN
---
>          IF (LMP2 .OR. LMP2KPAR .OR. LMP2NO .OR. LFCIDUMP  .OR. LRPAX .OR.  LCCSD .OR. LCC4S .OR. &
>             LBRACKETST .OR. FOURORBIT==1 ) THEN
EOF

md5sum_orig=71397c4fabdbee57861d2f8baba24fa8
md5sum_new=`md5sum main.F  | awk '{print $1}' `
if [ "$md5sum_orig" = "$md5sum_new" ]; then
    echo "Patching file  main.F " 
    patch main.F .tmp_vaspcc4s_patch 
else
   echo "You dont have the correct version of file  main.F  to apply patch. " 
fi
 
#######################
# diff of   ump2.F
#######################
 
 
cat > .tmp_vaspcc4s_patch<<"EOF"
11a12,13
> ! DESCRIPTION: 
> !
20,24d21
< ! There are however some caveats to this routine. 
< ! One problem is that the difference vectors between any two
< ! k-points must be included in the k-point set. This requires to
< ! use Gamma centered meshes.
< !
34a32,39
> ! ================================================
> ! THE FOLLOWING SETTINGS ARE SUPPORTED:
> ! ================================================
> !
> !   * RESTRICTED (RHF) AND UNRESTRICTED (UHF) CANONICAL ORBITALS ONLY
> !   * GAMMA-CENTERED OR SHIFTED K-MESHES
> !
> !
42,43c47,56
<       COMPLEX(q) , ALLOCATABLE, PRIVATE, SAVE :: FTOD_PW(:,:,:,:,:,:,:)      
<       GDEF      , ALLOCATABLE, PRIVATE, SAVE :: FTOD_OC(:,:,:,:,:,:,:)      
---
>       COMPLEX(q) , ALLOCATABLE, PRIVATE, SAVE :: FTOD_PW_POTFAK(:,:,:)
>       COMPLEX(q) , ALLOCATABLE, PRIVATE, SAVE :: FTOD_PW(:,:,:,:,:,:,:)
>       COMPLEX(q) , ALLOCATABLE, PRIVATE, SAVE :: POTFAK_FULL(:,:)
>       COMPLEX(q) , ALLOCATABLE, PRIVATE, SAVE :: POTFAK_FULL_DIST(:,:)
>       COMPLEX(q) , ALLOCATABLE, PRIVATE, SAVE :: SFACTOR_DIRECT(:,:)
>       COMPLEX(q) , ALLOCATABLE, PRIVATE, SAVE :: SFACTOR_EXCHANGE(:,:)
>       COMPLEX(q) , ALLOCATABLE, PRIVATE, SAVE :: SFACTOR_DIRECT_PAIR(:)
>       COMPLEX(q) , ALLOCATABLE, PRIVATE, SAVE :: SFACTOR_EXCHANGE_PAIR(:)
>       COMPLEX(q) , ALLOCATABLE, PRIVATE, SAVE :: SFACTOR_FULL(:,:)
>       GDEF      , ALLOCATABLE, PRIVATE, SAVE :: FTOD_OC(:,:,:,:,:,:,:)
45a59
>       GDEF , ALLOCATABLE, PRIVATE, SAVE :: TWOE4ORBITAL_FULL(:,:)
71c85
<       COMPLEX(q) :: E_MP2, energy_c, energy_x,energy_cc,EFOCK, tmp2
---
>       COMPLEX(q) :: E_MP2, energy_c, energy_x,energy_cc,EFOCK, energy_ctest
91c105
<       !Restrict band index to active space: NBANDSHIGH >= N > NFREEZE 
---
>       !Restrict band index to active space: NBANDSHIGH_ >= N > NFREEZE 
94c108
<       INTEGER :: NBANDSHIGH
---
>       INTEGER :: NBANDSHIGH_   ! comflicts with NBANDSHIGH from some used module
96,97c110,133
<       LOGICAL :: LDUMPPAIRS
<       COMPLEX(q), ALLOCATABLE :: EMP2_PAIR(:,:,:,:)
---
>       COMPLEX(q), ALLOCATABLE :: EMP2_PAIR_SINGLET_ENCUTS(:,:,:,:,:,:,:)
>       COMPLEX(q), ALLOCATABLE :: EMP2_PAIR_TRIPLET_ENCUTS(:,:,:,:,:,:,:)
>       COMPLEX(q), ALLOCATABLE :: EMP2_PAIR_SINGLET_CBS(:,:,:,:,:,:)
>       COMPLEX(q), ALLOCATABLE :: EMP2_PAIR_TRIPLET_CBS(:,:,:,:,:,:)
>       COMPLEX(q), ALLOCATABLE :: EMP2_SINGLET(:)
>       COMPLEX(q), ALLOCATABLE :: EMP2_TRIPLET(:)
>       COMPLEX(q) :: EMP2_SINGLET_CBS
>       COMPLEX(q) :: EMP2_TRIPLET_CBS
>       !Structure factor calculation and automatic basis set extrapolation variables
>       LOGICAL :: LSFACTOR
>       !energy cutoff extrapolation
>       !different ENCUTs and ENCUTSOFTs for extrapolation
>       INTEGER, PARAMETER, PRIVATE :: N_MP2_ENCUTS = 8        ! number of different ENCUTs
>       REAL(q), PARAMETER, PRIVATE :: ENCUT_DIVISOR = 1.05_q  ! ENCUTS(n+1) = ENCUTS(n) / ENCUT_DIVISOR
>       REAL(q), ALLOCATABLE, PRIVATE :: ENCUTS(:)             ! different ENCUTs
>       REAL(q), ALLOCATABLE, PRIVATE :: ENCUTSOFTS(:)         ! different ENCUTSOFTs
> #ifdef gammareal
>       INTEGER, PARAMETER :: BLOCKSIZE=64 !choose reasonable block size here.
> #else
>       INTEGER, PARAMETER :: BLOCKSIZE=64 !choose reasonable block size here.
>                                           ! intel compilers/libs 2018-2021 show
>                                           ! memory leaks except for certain
>                                           ! "magic" block sizes like 124
> #endif
138a175
>          INTEGER :: FTOD_PW_rows,NG,MB
152c189
<          NBANDSHIGH=WDES%NB_TOT
---
>          NBANDSHIGH_=WDES%NB_TOT
186,187c223,226
<          
<                   
---
>          CALL SETUP_POTFAK_FULL_DIST(WDES)
>          IF (LSFACTOR) CALL SETUP_ENCUTS(ENCUTGW,ENCUTGWSOFT)
>          IF (LSFACTOR) CALL SETUP_SFACTOR(WDES)
> 
188a228
> 
233a274,301
>                      IF (IO%IU0>=0) WRITE(IO%IU0,*) "NBI=",NBI
> 
>                      !APPLY POTFAK here
> #ifdef gammareal
>                      FTOD_PW_POTFAK(:,:,1)=FTOD_PW(:,:,NBI-NFREEZE,KI,KQ,ISP,1)
>                      FTOD_PW_POTFAK(:,:,2)=FTOD_PW(:,:,NBI-NFREEZE,KI,KQ_,ISP,1)
> #else
>                      FTOD_PW_POTFAK(:,:,1)=FTOD_PW(:,:,NBI-NFREEZE,RKIofKI(KI),RKQofKQ(KQ),ISP,1) !FTOD_PW(:,:,NBI-NFREEZE,KI,KQ,ISP,1)
>                      FTOD_PW_POTFAK(:,:,2)=FTOD_PW(:,:,NBI-NFREEZE,RKIofKI(KI),RKQofKQ(KQ_),ISP,1) !FTOD_PW(1,1,NBI-NFREEZE,KI,KQ_,ISP,1)
> #endif
> 
> !                    call BLACS_GRIDINFO(CONTXT_GRID, NPROW, NPCOL, MYROW, MYCOL)
> !         !Blocking size for block-cyclic distribution of FTOD_PW        
>                      FTOD_PW_rows = SIZE(FTOD_PW_POTFAK,1) !numroc(NGVECTOR,MIN(((NGVECTOR)/NPROW),50),MYROW,0,NPROW)
> !                     IF ((NBI==1)) WRITE(*,*)'ME, FTOD rows',ME,FTOD_PW_rows
>  
>                      DO NG=1,FTOD_PW_rows
> #ifdef gammareal
>                         FTOD_PW_POTFAK(NG,:,1)=POTFAK_FULL_DIST(NG,(KQ))*FTOD_PW_POTFAK(NG,:,1)
>                         FTOD_PW_POTFAK(NG,:,2)=POTFAK_FULL_DIST(NG,(KQ_))*FTOD_PW_POTFAK(NG,:,2)
> #else
>                         FTOD_PW_POTFAK(NG,:,1)=POTFAK_FULL_DIST(NG,RKQofKQ(KQ))*FTOD_PW_POTFAK(NG,:,1)
>                         FTOD_PW_POTFAK(NG,:,2)=POTFAK_FULL_DIST(NG,RKQofKQ(KQ_))*FTOD_PW_POTFAK(NG,:,2)
> #endif
>                      ENDDO
> 
> 
> 
249c317
<                           m_ NGVECTOR,-one, FTOD_PW(1,1,NBI-NFREEZE,KI,KQ,ISP,1),1,1,&
---
>                           m_ NGVECTOR,-one, FTOD_PW_POTFAK(1,1,1),1,1,&
254c322
<                           NGVECTOR,-one, FTOD_PW(1,1,NBI-NFREEZE,RKIofKI(KI),RKQofKQ(KQ),ISP,1),1,1,&
---
>                           NGVECTOR,-one, FTOD_PW_POTFAK(1,1,1),1,1,&
297c365
<                              desc_FTOD_PW,FTOD_PW(1,1,NBI-NFREEZE,KI,KQ_,ISP,1),1,1,&
---
>                              desc_FTOD_PW,FTOD_PW_POTFAK(1,1,2),1,1,&
304c372
<                              desc_FTOD_PW,FTOD_PW(1,1,NBI-NFREEZE,RKIofKI(KI),RKQofKQ(KQ_),ISP,1),1,1,&
---
>                              desc_FTOD_PW,FTOD_PW_POTFAK(1,1,2),1,1,&
331a400,444
> 
>                         IF (LSFACTOR) THEN
>                            !direct contribition to structure factor
>                            TWOE4ORBITAL_FULL(:,:)=GCONJG(TWOE4ORBITAL(:,:)) !-TWOE4ORBITAL_X(:,:)
>                            CALL APPLY_DENOM_TWOE4ORBITAL_FULL(W,KI,KJ,KA,KB,NBI,NBJ,ISP,ISP,LATT_CUR,WDES%ISPIN)
>                            ! contract na
> #ifdef gammareal
>                            CALL PDGEMM('n','n',(m_ NGVECTOR),(PROCS*WDES%NBANDS),&
>                              (PROCS*WDES%NBANDS),-one, FTOD_PW(1,1,NBI-NFREEZE,KI,KQ,ISP,1),1,1,&
>                              desc_FTOD_PW,TWOE4ORBITAL_FULL(1,1),1,1,&
>                              desc_TWOE4ORBITAL,zero, FTOD_PW_POTFAK(1,1,3),1,1,desc_FTOD_PW)
> #else
>                            CALL PZGEMM('n','n',(NGVECTOR),(PROCS*WDES%NBANDS),&
>                              (PROCS*WDES%NBANDS),-one, FTOD_PW(1,1,NBI-NFREEZE,RKIofKI(KI),RKQofKQ(KQ),ISP,1),1,1,&
>                              desc_FTOD_PW,TWOE4ORBITAL_FULL(1,1),1,1,&
>                              desc_TWOE4ORBITAL,zero, FTOD_PW_POTFAK(1,1,3),1,1,desc_FTOD_PW)
> #endif
> 
>                            ! contract nb and add to S_ij(G), S_direct,
>                            ! S_exchange abd S_full
>                            call CONTRACT_NB_ADD_DIRECT_TO_SFACTOR(W,WDES,KI,KJ,KA,KB,NBI,NBJ,ISP,ISP,LATT_CUR,WDES%ISPIN)
> 
>                            !exchange contribition to structure factor
>                            TWOE4ORBITAL_FULL(:,:)=GCONJG(TWOE4ORBITAL_X(:,:)) !-TWOE4ORBITAL_X(:,:)
>                            CALL APPLY_DENOM_TWOE4ORBITAL_FULL(W,KI,KJ,KA,KB,NBI,NBJ,ISP,ISP,LATT_CUR,WDES%ISPIN)
>                            ! contract na
> #ifdef gammareal
>                            CALL PDGEMM('n','n',(m_ NGVECTOR),(PROCS*WDES%NBANDS),&
>                              (PROCS*WDES%NBANDS),-one, FTOD_PW(1,1,NBI-NFREEZE,KI,KQ,ISP,1),1,1,&
>                              desc_FTOD_PW,TWOE4ORBITAL_FULL(1,1),1,1,&
>                              desc_TWOE4ORBITAL,zero, FTOD_PW_POTFAK(1,1,3),1,1,desc_FTOD_PW)
> #else
>                            CALL PZGEMM('n','n',(NGVECTOR),(PROCS*WDES%NBANDS),&
>                              (PROCS*WDES%NBANDS),-one, FTOD_PW(1,1,NBI-NFREEZE,RKIofKI(KI),RKQofKQ(KQ),ISP,1),1,1,&
>                              desc_FTOD_PW,TWOE4ORBITAL_FULL(1,1),1,1,&
>                              desc_TWOE4ORBITAL,zero, FTOD_PW_POTFAK(1,1,3),1,1,desc_FTOD_PW)
> #endif
> 
>                            ! contract nb and add to S_ij(G), S_direct,
>                            ! S_exchange abd S_full
> 
>                            CALL CONTRACT_NB_ADD_EXCHANGE_TO_SFACTOR(W,WDES,KI,KJ,KA,KB,NBI,NBJ,ISP,ISP,LATT_CUR,WDES%ISPIN)
> 
>                            CALL GET_PAIR_CBS_ENERGY_FOR_KQ(WDES,WGW,LATT_CUR,KI,KJ,KQ,NBI,NBJ,ISP,ISP,FSG_STORE(1))
>                         ENDIF
372a486,508
> 
>                         IF (LSFACTOR) THEN
>                            !direct contribition to structure factor
>                            TWOE4ORBITAL_FULL(:,:)=GCONJG(TWOE4ORBITAL(:,:)) !-TWOE4ORBITAL_X(:,:)
>                            CALL APPLY_DENOM_TWOE4ORBITAL_FULL(W,KI,KJ,KA,KB,NBI,NBJ,ISP,ISP-1,LATT_CUR,WDES%ISPIN)
>                            ! contract na
> #ifdef gammareal
>                            CALL PDGEMM('n','n',(m_ NGVECTOR),(PROCS*WDES%NBANDS),&
>                              (PROCS*WDES%NBANDS),-one, FTOD_PW(1,1,NBI-NFREEZE,KI,KQ,ISP,1),1,1,&
>                              desc_FTOD_PW,TWOE4ORBITAL_FULL(1,1),1,1,&
>                              desc_TWOE4ORBITAL,zero, FTOD_PW_POTFAK(1,1,3),1,1,desc_FTOD_PW)
> #else
>                            CALL PZGEMM('n','n',(NGVECTOR),(PROCS*WDES%NBANDS),&
>                              (PROCS*WDES%NBANDS),-one, FTOD_PW(1,1,NBI-NFREEZE,RKIofKI(KI),RKQofKQ(KQ),ISP,1),1,1,&
>                              desc_FTOD_PW,TWOE4ORBITAL_FULL(1,1),1,1,&
>                              desc_TWOE4ORBITAL,zero, FTOD_PW_POTFAK(1,1,3),1,1,desc_FTOD_PW)
> #endif
> 
>                            ! contract nb and add to S_ij(G)
>                            call CONTRACT_NB_ADD_DIRECT_TO_SFACTOR(W,WDES,KI,KJ,KA,KB,NBI,NBJ,ISP,ISP-1,LATT_CUR,WDES%ISPIN)
> 
>                         ENDIF
> 
384c520,522
<          
---
>         
> 
> 
386,393c524,531
<          CALL M_sum_z(WGW%COMM_INTER, energy_c, 1)
<          CALL M_sum_z(WGW%COMM_INTER, energy_cc, 1)
<          CALL M_sum_z(WGW%COMM_INTER, energy_x, 1)
<          CALL M_sum_z(WGW%COMM_INTER, EFOCK, 1)
<          E_MP2=energy_c+energy_x
<          IF (LDUMPPAIRS) THEN
<             CALL M_sum_z(WGW%COMM_INTER, EMP2_PAIR(1,1,1,1), VBMAX*VBMAX*WDES%ISPIN*WDES%ISPIN)
<          ENDIF
---
> 
> !         IF (.not. LSFACTOR) THEN
>             CALL M_sum_z(WGW%COMM_INTER, energy_c, 1)
>             CALL M_sum_z(WGW%COMM_INTER, energy_cc, 1)
>             CALL M_sum_z(WGW%COMM_INTER, energy_x, 1)
>             CALL M_sum_z(WGW%COMM_INTER, EFOCK, 1)
>             E_MP2=energy_c+energy_x
> !         ENDIF
428d565
<                ENDIF
430,434c567,572
<                WRITE(IO%IU6,*)
<                WRITE(IO%IU6,*) 'Moeller Plesset 2 correlation:'
<                WRITE(IO%IU6,*) '================================'
<                WRITE(IO%IU6,11) REAL(EFOCK,Kind=q),REAL(energy_c,Kind=q),REAL(energy_x,Kind=q),REAL(E_MP2,Kind=q)
<                WRITE(IO%IU6,*)         
---
>                   WRITE(IO%IU6,*)
>                   WRITE(IO%IU6,*) 'Moeller Plesset 2 correlation:'
>                   WRITE(IO%IU6,*) '================================'
>                   WRITE(IO%IU6,11) REAL(EFOCK,Kind=q),REAL(energy_c,Kind=q),REAL(energy_x,Kind=q),REAL(E_MP2,Kind=q)
>                   WRITE(IO%IU6,*)         
>                ENDIF
436c574
< 11          FORMAT('     Hartree Fock energy: ',F20.8/&
---
> 11          FORMAT('    Fock exchange energy: ',F20.8/&
460a599
>          ENDIF
462,478d600
<             IF (LDUMPPAIRS) THEN
<                WRITE(IO%IU6,*)
<                WRITE(IO%IU6,*)"Writing out MP2 pair energies:"
<                WRITE(IO%IU6,*)"# EMP2  NBI  NBJ  ISPI ISPJ"
< 
<                DO ISPI=1,WDES%ISPIN
<                DO ISPJ=1,WDES%ISPIN
<                DO NBI=1,VBMAX
<                DO NBJ=1,VBMAX
< 
<                WRITE(IO%IU6,*)EMP2_PAIR(NBI,NBJ,ISPI,ISPJ),NBI,NBJ,ISPI,ISPJ
<                                     
<                ENDDO
<                ENDDO
<                ENDDO
<                ENDDO
<             ENDIF
479a602,606
>          IF (LSFACTOR) THEN
>             IF (IO%IU0>0) WRITE(IO%IU0,*)''
>             IF (IO%IU0>0) WRITE(IO%IU0,*)'Automatic extrapolation to the CBS limit will be performed'
>             CALL DUMP_CORR_ENERGIES(IO,WDES,WGW,LATT_CUR,ENCUTGW,ENCUTGWSOFT,FSG_STORE(1))
> !               WRITE(IO%IU6,11) REAL(EFOCK,Kind=q),REAL(energy_c,Kind=q),REAL(energy_x,Kind=q),REAL(E_MP2,Kind=q)
480a608
> 
826c954
<                IF (RNB>NBANDSHIGH) THEN
---
>                IF (RNB>NBANDSHIGH_) THEN
831c959
<                IF (RNA>NBANDSHIGH) THEN
---
>                IF (RNA>NBANDSHIGH_) THEN
880,882d1007
<                      IF (LDUMPPAIRS) THEN
<                         EMP2_PAIR(NI,NJ,ISP1,ISP2)=EMP2_PAIR(NI,NJ,ISP1,ISP2)+OCC*VIRT*FAC1*(FAC2*Inte_c)/DENOM
<                      ENDIF
905,907d1029
<                      IF (LDUMPPAIRS) THEN
<                         EMP2_PAIR(NI,NJ,ISP1,ISP2)=EMP2_PAIR(NI,NJ,ISP1,ISP2)-OCC*VIRT*FAC1*Inte_x/DENOM
<                      ENDIF
936,938d1057
<                      IF (LDUMPPAIRS) THEN
<                         EMP2_PAIR(NI,NJ,ISP1,ISP2)=EMP2_PAIR(NI,NJ,ISP1,ISP2)+OCC*VIRT*Inte_c/DENOM
<                      ENDIF
952c1071
<                IF (RNB>NBANDSHIGH) THEN
---
>                IF (RNB>NBANDSHIGH_) THEN
958c1077
<                IF (RNA>NBANDSHIGH) THEN
---
>                IF (RNA>NBANDSHIGH_) THEN
989,991d1107
<                      IF (LDUMPPAIRS) THEN
<                         EMP2_PAIR(NI,NJ,ISP1,ISP2)=EMP2_PAIR(NI,NJ,ISP1,ISP2)+OCC*VIRT*FAC1*(FAC2*Inte_c)/DENOM-OCC*VIRT*FAC1*(Inte_x)/DENOM
<                      ENDIF
1001,1004d1116
<                      IF (LDUMPPAIRS) THEN
<                         EMP2_PAIR(NI,NJ,ISP1,ISP2)=EMP2_PAIR(NI,NJ,ISP1,ISP2)+OCC*VIRT*(Inte_c)/DENOM
<                      ENDIF
< 
1017a1130,1656
> 
>       SUBROUTINE APPLY_DENOM_TWOE4ORBITAL_FULL(W,KI,KJ,KA,KB,NI,NJ,ISP1,ISP2,LATT_CUR,WISPIN)
>          USE constant
>          USE full_kpoints
>          USE mkpoints
>          USE wave
>          IMPLICIT NONE 
>          TYPE (wavespin) W
>          INTEGER :: KI,KJ,KA,KB,NI,NJ,I,WISPIN
>          TYPE(latt) LATT_CUR
>          INTEGER :: NA,NB,TWOE4ORBITAL_ROWS,TWOE4ORBITAL_COLS
>          INTEGER :: RNA,RNB,ISP1,ISP2,KI_IN_FULL_ORIG,KJ_IN_FULL_ORIG,kq1,kq2
>          REAL(q) :: DENOM, OCC, VIRT
>          COMPLEX(q) :: Inte_c
>          
>          TWOE4ORBITAL_ROWS = numroc(desc_TWOE4ORBITAL(3),desc_TWOE4ORBITAL(5),MYROW,0,NPROW)
>          TWOE4ORBITAL_COLS = numroc(desc_TWOE4ORBITAL(4),desc_TWOE4ORBITAL(6),MYCOL,0,NPCOL)
>          
>          !KI_IN_FULL_ORIG=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,KI),KPOINTS_FULL_ORIG) 
>          KJ_IN_FULL_ORIG=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,KJ),KPOINTS_FULL_ORIG)
>          
>          kq1=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,KB)-W%WDES%VKPT(:,KJ),KPOINTS_FULL)
>          ! k_a = k_i + k_q - G
>          KQ2=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,KI)-W%WDES%VKPT(:,KA),KPOINTS_FULL)
>          
>          if (kq1/=kq2) WRITE(*,*)'error: q-point changed.'
>          !write(*,*)'kq1',W%WDES%VKPT(:,kq1)
>          KQ2=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,KI)-W%WDES%VKPT(:,KB),KPOINTS_FULL)
> 
> 
>          DO NB=1,TWOE4ORBITAL_COLS
>             CALL LOC2GLOB(NB,MYCOL,desc_TWOE4ORBITAL(4),NPCOL,desc_TWOE4ORBITAL(6),RNB)
>                IF (RNB>NBANDSHIGH_) THEN
>                   TWOE4ORBITAL_FULL(NA,NB)=zero
>                   CYCLE
>                ENDIF
> 
>             DO NA=1,TWOE4ORBITAL_ROWS
>                CALL LOC2GLOB(NA,MYROW,desc_TWOE4ORBITAL(3),NPROW,desc_TWOE4ORBITAL(5),RNA)     
>                IF (RNA>NBANDSHIGH_) THEN
>                   TWOE4ORBITAL_FULL(NA,NB)=zero
>                   CYCLE
>                ENDIF
> 
>                OCC =W%FERTOT(NI,KI,ISP1)*W%FERTOT(NJ,KJ,ISP2)*KPOINTS_ORIG%WTKPT(KI)
>                VIRT=(1._q-W%FERTOT(RNA,KA,ISP1))*(1._q-W%FERTOT(RNB,KB,ISP2))
> 
> 
>                IF (FILLED_MP2_ORBITAL(W%FERTOT(RNA,KA,ISP1))) VIRT=0.0_q
>                IF (FILLED_MP2_ORBITAL(W%FERTOT(RNB,KB,ISP2))) VIRT=0.0_q
> 
>                TWOE4ORBITAL_FULL(NA,NB)=OCC*VIRT*(TWOE4ORBITAL_FULL(NA,NB))
> 
>                DENOM=REAL(W%CELTOT(NI,KI,ISP1),KIND=q)+&
>                  REAL(W%CELTOT(NJ,KJ,ISP2),KIND=q)-REAL(W%CELTOT(RNB,KB,ISP2),KIND=q)-&
>                  REAL(W%CELTOT(RNA,KA,ISP1),KIND=q)
> ! not gamma point
> #ifdef gammareal
>          IF ((NJ<NI) .and. (WISPIN==1)) DENOM=DENOM*0.5_q
> #else
> 
> #endif
> 
> 
>              IF (DENOM<0._q) THEN
>                IF (ISP2==ISP1) THEN               
> !                  TWOE4ORBITAL_FULL(NA,NB)=OCC*VIRT*FAC1*(FAC2*TWOE4ORBITAL_FULL(NA,NB))/DENOM
>                   Inte_c=TWOE4ORBITAL(NA,NB)*GCONJG(TWOE4ORBITAL(NA,NB))*OCC*VIRT/DENOM
>                   TWOE4ORBITAL_FULL(NA,NB)=OCC*VIRT*(TWOE4ORBITAL_FULL(NA,NB))/DENOM
>                   Energy_ctest=Energy_c+OCC*VIRT*(Inte_c)/DENOM
> 
> !                     Inte_c=zero
> !                     Inte_c=TWOE4ORBITAL(NA,NB)*GCONJG(TWOE4ORBITAL(NA,NB))             
> !                     Energy_c=Energy_c+OCC*VIRT*FAC1*(FAC2*Inte_c)/DENOM
> !                     Inte_x=zero
> !                     Inte_x=Inte_x+GCONJG(TWOE4ORBITAL(NA,NB))*(TWOE4ORBITAL_X(NA,NB))
> !                     Energy_x=Energy_x-OCC*VIRT*FAC1*(Inte_x)/DENOM                  
> 
>                ENDIF
>                
>                IF (ISP2/=ISP1) THEN
>                   TWOE4ORBITAL_FULL(NA,NB)=OCC*VIRT*(TWOE4ORBITAL_FULL(NA,NB))/DENOM
> !                     Inte_c=zero
> !                     Energy_cc=Energy_cc+OCC*VIRT*(Inte_c)/DENOM
> 
> 
>                ENDIF !isp1/=isp2
> 
>              ENDIF !denom<0 
>                
>            ENDDO
>          ENDDO
> 
> 
>       END SUBROUTINE APPLY_DENOM_TWOE4ORBITAL_FULL
> 
> 
>       SUBROUTINE CONTRACT_NB_ADD_DIRECT_TO_SFACTOR(W,WDES,KI,KJ,KA,KB,NI,NJ,ISP1,ISP2,LATT_CUR,WISPIN)
>          USE constant
>          USE full_kpoints
>          USE mkpoints
>          USE wave
>          IMPLICIT NONE 
>          TYPE (wavespin) W
>          TYPE(wavedes) WDES
>          INTEGER :: KI,KJ,KA,KB,KQ,NI,NJ,ISP1,ISP2,WISPIN
>          TYPE(latt) LATT_CUR
>          INTEGER :: NG,NA,NB,TWOE4ORBITAL_ROWS,TWOE4ORBITAL_COLS
>          INTEGER :: RNG,RNA,RNB,KI_IN_FULL_ORIG,KJ_IN_FULL_ORIG,kq1,kq2
>          REAL(q) :: DENOM, OCC, VIRT
>          INTEGER :: FTOD_PW_rows, FTOD_PW_cols,cc,FTOD_OC_rows, FTOD_OC_cols
> 
>          call BLACS_GRIDINFO(CONTXT_GRID, NPROW, NPCOL, MYROW, MYCOL)
>          !Blocking size for block-cyclic distribution of FTOD_PW        
>          MB=MIN(((NGVECTOR)/NPROW),BLOCKSIZE)   !Row blocking size
>          NB=MIN((NPROW*WDES%NBANDS),MB)   !column blocking size          
>          
>          FTOD_PW_rows = numroc(NGVECTOR,MB,MYROW,0,NPROW)
>          FTOD_PW_cols = numroc(PROCS*WDES%NBANDS,NB,MYCOL,0,NPCOL)
>          
>          !KI_IN_FULL_ORIG=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,KI),KPOINTS_FULL_ORIG) 
>          KJ_IN_FULL_ORIG=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,KJ),KPOINTS_FULL_ORIG)
>          
>          ! k_a = k_i + k_q - G
>          KQ=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,KI)-W%WDES%VKPT(:,KA),KPOINTS_FULL)
>          
>          if (kq1/=kq2) WRITE(*,*)'error: q-point changed.'
>          KQ2=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,KI)-W%WDES%VKPT(:,KB),KPOINTS_FULL)
> 
> 
>          SFACTOR_DIRECT_PAIR=(0.0_q,0.0_q)
> 
>          DO NB=1,FTOD_PW_COLS
>             CALL LOC2GLOB(NB,MYCOL,desc_FTOD_PW(4),NPCOL,desc_FTOD_PW(6),RNB)
>             IF (RNB>NBANDSHIGH_) THEN
>                CYCLE
>             ENDIF
> 
>             OCC =W%FERTOT(NI,KI,ISP1)*W%FERTOT(NJ,KJ,ISP2)*KPOINTS_ORIG%WTKPT(KI)
>             VIRT=(1._q-W%FERTOT(RNB,KB,ISP2))
> 
>             IF (FILLED_MP2_ORBITAL(W%FERTOT(RNB,KB,ISP2))) VIRT=0.0_q
> 
>             IF (ISP2==ISP1) THEN               
>                DO NG=1,FTOD_PW_rows
>                   CALL LOC2GLOB(NG,MYROW,NGVECTOR,NPROW,MB,RNG)
> !                  SFACTOR_FULL(RNG,RKQofKQ(KQ))=SFACTOR_FULL(RNG,RKQofKQ(KQ))+FAC2*CONJG(FTOD_PW_POTFAK(NG,NB,3))*(FTOD_PW(NG,NB,NJ-NFREEZE,RKIofKI(KJ),RKQofKQ(KQ),ISP2,ncc))*VIRT*OCC
>                   SFACTOR_DIRECT_PAIR(RNG)=SFACTOR_DIRECT_PAIR(RNG)+FAC2*CONJG(FTOD_PW_POTFAK(NG,NB,3))*(FTOD_PW(NG,NB,NJ-NFREEZE,RKIofKI(KJ),RKQofKQ(KQ),ISP2,ncc))*VIRT*OCC
> !                  SFACTOR_DIRECT(RNG,RKQofKQ(KQ))=SFACTOR_DIRECT(RNG,RKQofKQ(KQ))+FAC2*CONJG(FTOD_PW_POTFAK(NG,NB,3))*(FTOD_PW(NG,NB,NJ-NFREEZE,RKIofKI(KJ),RKQofKQ(KQ),ISP2,ncc))*VIRT*OCC
>                ENDDO
>             ENDIF
>                
>             IF (ISP2/=ISP1) THEN
> 
>             ENDIF !isp1/=isp2
>          ENDDO
> 
>          DO NG=1,FTOD_PW_rows
>             CALL LOC2GLOB(NG,MYROW,NGVECTOR,NPROW,MB,RNG)
>             SFACTOR_FULL(RNG,RKQofKQ(KQ))=SFACTOR_FULL(RNG,RKQofKQ(KQ))+SFACTOR_DIRECT_PAIR(RNG)
>             SFACTOR_DIRECT(RNG,RKQofKQ(KQ))=SFACTOR_DIRECT(RNG,RKQofKQ(KQ))+SFACTOR_DIRECT_PAIR(RNG)
>          ENDDO
> 
>       END SUBROUTINE CONTRACT_NB_ADD_DIRECT_TO_SFACTOR
> 
>       SUBROUTINE CONTRACT_NB_ADD_EXCHANGE_TO_SFACTOR(W,WDES,KI,KJ,KA,KB,NI,NJ,ISP1,ISP2,LATT_CUR,WISPIN)
>          USE constant
>          USE full_kpoints
>          USE mkpoints
>          USE wave
>          IMPLICIT NONE 
>          TYPE (wavespin) W
>          TYPE(wavedes) WDES
>          INTEGER :: KI,KJ,KA,KB,KQ,NI,NJ,ISP1,ISP2,WISPIN
>          TYPE(latt) LATT_CUR
>          INTEGER :: NG,NA,NB,TWOE4ORBITAL_ROWS,TWOE4ORBITAL_COLS
>          INTEGER :: RNG,RNA,RNB,KI_IN_FULL_ORIG,KJ_IN_FULL_ORIG,kq1,kq2
>          REAL(q) :: DENOM, OCC, VIRT
>          INTEGER :: FTOD_PW_rows, FTOD_PW_cols,cc,FTOD_OC_rows, FTOD_OC_cols
> 
>          call BLACS_GRIDINFO(CONTXT_GRID, NPROW, NPCOL, MYROW, MYCOL)
>          !Blocking size for block-cyclic distribution of FTOD_PW        
>          MB=MIN(((NGVECTOR)/NPROW),BLOCKSIZE)   !Row blocking size
>          NB=MIN((NPROW*WDES%NBANDS),MB)   !column blocking size          
>          
>          FTOD_PW_rows = numroc(NGVECTOR,MB,MYROW,0,NPROW)
>          FTOD_PW_cols = numroc(PROCS*WDES%NBANDS,NB,MYCOL,0,NPCOL)
>          
>          !KI_IN_FULL_ORIG=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,KI),KPOINTS_FULL_ORIG) 
>          KJ_IN_FULL_ORIG=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,KJ),KPOINTS_FULL_ORIG)
>          
>          ! k_a = k_i + k_q - G
>          KQ=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,KI)-W%WDES%VKPT(:,KA),KPOINTS_FULL)
>          
>          if (kq1/=kq2) WRITE(*,*)'error: q-point changed.'
>          !write(*,*)'kq1',W%WDES%VKPT(:,kq1)
>          KQ2=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,KI)-W%WDES%VKPT(:,KB),KPOINTS_FULL)
> 
>          SFACTOR_EXCHANGE_PAIR=(0.0_q,0.0_q)
> 
>          DO NB=1,FTOD_PW_COLS
>             CALL LOC2GLOB(NB,MYCOL,desc_FTOD_PW(4),NPCOL,desc_FTOD_PW(6),RNB)
>             IF (RNB>NBANDSHIGH_) THEN
>                CYCLE
>             ENDIF
> 
>             OCC =W%FERTOT(NI,KI,ISP1)*W%FERTOT(NJ,KJ,ISP2)*KPOINTS_ORIG%WTKPT(KI)
>             VIRT=(1._q-W%FERTOT(RNB,KB,ISP2))
> 
>             IF (FILLED_MP2_ORBITAL(W%FERTOT(RNB,KB,ISP2))) VIRT=0.0_q
> 
>             IF (ISP2==ISP1) THEN               
>                DO NG=1,FTOD_PW_rows
>                   !CALL LOC2GLOB(NG,MYROW,desc_FTOD_PW(3),NPROW,desc_FTOD_PW(5),RNG)
>                   CALL LOC2GLOB(NG,MYROW,NGVECTOR,NPROW,MB,RNG)
> !!                  SFACTOR_FULL(RNG,RKQofKQ(KQ))=SFACTOR_FULL(RNG,RKQofKQ(KQ))-FTOD_PW_POTFAK(NG,NB,3)*CONJG(FTOD_PW(NG,NB,NJ-NFREEZE,RKIofKI(KJ),RKQofKQ(KQ),ISP2,ncc))*VIRT*OCC
> 
>                   SFACTOR_EXCHANGE_PAIR(RNG)=SFACTOR_EXCHANGE_PAIR(RNG)-FTOD_PW_POTFAK(NG,NB,3)*CONJG(FTOD_PW(NG,NB,NJ-NFREEZE,RKIofKI(KJ),RKQofKQ(KQ),ISP2,ncc))*VIRT*OCC
>                ENDDO
>             ENDIF
>                
>             IF (ISP2/=ISP1) THEN
> 
> 
>             ENDIF !isp1/=isp2
> 
>          ENDDO
> 
>          DO NG=1,FTOD_PW_rows
>             !CALL LOC2GLOB(NG,MYROW,desc_FTOD_PW(3),NPROW,desc_FTOD_PW(5),RNG)
>             CALL LOC2GLOB(NG,MYROW,NGVECTOR,NPROW,MB,RNG)
>             SFACTOR_FULL(RNG,RKQofKQ(KQ))=SFACTOR_FULL(RNG,RKQofKQ(KQ))+SFACTOR_EXCHANGE_PAIR(RNG)
>             SFACTOR_EXCHANGE(RNG,RKQofKQ(KQ))=SFACTOR_EXCHANGE(RNG,RKQofKQ(KQ))+SFACTOR_EXCHANGE_PAIR(RNG)
>          ENDDO
> 
>       END SUBROUTINE CONTRACT_NB_ADD_EXCHANGE_TO_SFACTOR
> 
>       SUBROUTINE GET_PAIR_CBS_ENERGY_FOR_KQ(WDES,WGW,LATT_CUR,KI,KJ,KQ,NBI,NBJ,ISP1,ISP2,FSG)
>          USE prec
>          USE poscar
>          USE pseudo
>          USE wave_high
>          USE full_kpoints
>          USE mkpoints
>          USE lattice
>          USE constant
>          USE base
>          IMPLICIT NONE
>          TYPE(wavedes) WDES
>          TYPE(wavedes) WGW
>          TYPE(latt) LATT_CUR
>          INTEGER :: KI,KJ,KQ,NBI,NBJ,ISP1,ISP2
>          REAL(q) :: FSG ! singularity correction
>          ! local variables
>          REAL(q) :: POTFAK(GRIDHF%MPLWV)
>          COMPLEX(q) :: TMPSFAC(GRIDHF%MPLWV)
>          TYPE(wavedes1) WGWQ
>          INTEGER :: NG, I, RKI, RKJ
>          REAL(q) :: E_SINGLET(N_MP2_ENCUTS), E_TRIPLET(N_MP2_ENCUTS)
>          REAL(q) :: E_SINGLET_CBS, E_TRIPLET_CBS
>  
> 
>          CALL M_sum_z(WGW%COMM_INTER, SFACTOR_DIRECT_PAIR(1), NGVECTOR)
>          CALL M_sum_z(WGW%COMM_INTER, SFACTOR_EXCHANGE_PAIR(1), NGVECTOR)
> 
>          E_SINGLET=(0.0_q)
>          E_TRIPLET=(0.0_q)
> 
>          DO I=1,N_MP2_ENCUTS
>             CALL SETWDES(WGW,WGWQ,KQ)
> !            IF (ENCUTS(I) /= ENCUTSOFTS(I) .AND. ENCUTS(I) > 0 .AND. ENCUTSOFTS(I) > 0) THEN
>             POTFAK=(0.0_q)
>             CALL SET_GFAC_WAVEFUN(WGWQ,LATT_CUR,FSG,POTFAK,ENCUTS(I),ENCUTSOFTS(I))
> !            ELSE
> !               CALL SET_GFAC_WAVEFUN(WGWQ,LATT_CUR,FSG,POTFAK)
> !            ENDIF
>             IF (ISP2 == ISP1) THEN
>                IF (MCALPHA/=0) THEN
>                   TMPSFAC=(0.0_q,0.0_q)
>                   TMPSFAC(1:NGVECTOR)=REAL((0.25_q*SFACTOR_DIRECT_PAIR(1:NGVECTOR)-0.5_q*SFACTOR_EXCHANGE_PAIR(1:NGVECTOR)),kind=q)
>                   !  include multipole correction
>                   CALL APPLY_GFAC_MULTIPOLE_WAVEFUN(WGWQ,TMPSFAC(1), &
>                        POTFAK(1))
>                   DO NG=1,NGVECTOR
>                      E_SINGLET(I)=E_SINGLET(I)+REAL(TMPSFAC(NG),kind=q)
>                   ENDDO
> 
>                   TMPSFAC=(0.0_q,0.0_q)
>                   TMPSFAC(1:NGVECTOR)=E_TRIPLET(I)+REAL((0.75_q*SFACTOR_DIRECT_PAIR(NG)+1.5_q*SFACTOR_EXCHANGE_PAIR(NG)),kind=q)
>                   !  include multipole correction
>                   CALL APPLY_GFAC_MULTIPOLE_WAVEFUN(WGWQ,TMPSFAC(1), &
>                        POTFAK(1))
>                   DO NG=1,NGVECTOR
>                      E_TRIPLET(I)=E_TRIPLET(I)+REAL(TMPSFAC(NG),kind=q)
>                   ENDDO
> 
>                ELSE
>                   DO NG=1,NGVECTOR
> !                  WRITE(*,*)'potfal,sfac',SQRT(1.0_q/POTFAK(NG)),REAL(SFACTOR_DIRECT(NG,RKQofKQ(KQ)),kind=q)
>                      E_SINGLET(I)=E_SINGLET(I)+REAL((0.25_q*SFACTOR_DIRECT_PAIR(NG)-0.5_q*SFACTOR_EXCHANGE_PAIR(NG)),kind=q)*POTFAK(NG)
>                      E_TRIPLET(I)=E_TRIPLET(I)+REAL((0.75_q*SFACTOR_DIRECT_PAIR(NG)+1.5_q*SFACTOR_EXCHANGE_PAIR(NG)),kind=q)*POTFAK(NG)
>                   ENDDO
>                ENDIF
>             ENDIF !isp1=isp2
> 
> !            IF (ISP2/=ISP1) THEN
> !            ENDIF !isp1/=isp2
> !            WRITE(*,'(A,2I4,2F14.6)')'ni,nj,encut,emp2',NBI,NBJ,ENCUTS(I),E_SINGLET(I)+E_TRIPLET(I)
> 
>          ENDDO
> 
>          CALL LIN_REG_CONVERGE_TRIPLET(E_TRIPLET, ENCUTS, E_TRIPLET_CBS)
>          CALL LIN_REG_CONVERGE_SINGLET(E_SINGLET, ENCUTS, E_SINGLET_CBS)
> !         WRITE(*,*)'ni,nj,converged,emp2',NBI,NBJ,'inf',E_SINGLET_CBS,E_TRIPLET_CBS
> 
>          RKI=RKIofKI(KI)
>          RKJ=RKIofKI(KJ)
> 
>          EMP2_PAIR_SINGLET_ENCUTS(NBI-NFREEZE,NBJ-NFREEZE,RKI,RKJ,ISP1,ISP2,:)=EMP2_PAIR_SINGLET_ENCUTS(NBI-NFREEZE,NBJ-NFREEZE,RKI,RKJ,ISP1,ISP2,:)+E_SINGLET(:)*REALNKPTS
>          EMP2_PAIR_TRIPLET_ENCUTS(NBI-NFREEZE,NBJ-NFREEZE,RKI,RKJ,ISP1,ISP2,:)=EMP2_PAIR_TRIPLET_ENCUTS(NBI-NFREEZE,NBJ-NFREEZE,RKI,RKJ,ISP1,ISP2,:)+E_TRIPLET(:)*REALNKPTS
>          EMP2_PAIR_SINGLET_CBS(NBI-NFREEZE,NBJ-NFREEZE,RKI,RKJ,ISP1,ISP2)=EMP2_PAIR_SINGLET_CBS(NBI-NFREEZE,NBJ-NFREEZE,RKI,RKJ,ISP1,ISP2)+E_SINGLET_CBS*REALNKPTS
>          EMP2_PAIR_TRIPLET_CBS(NBI-NFREEZE,NBJ-NFREEZE,RKI,RKJ,ISP1,ISP2)=EMP2_PAIR_TRIPLET_CBS(NBI-NFREEZE,NBJ-NFREEZE,RKI,RKJ,ISP1,ISP2)+E_TRIPLET_CBS*REALNKPTS
> 
> !          WRITE(*,*)'ni,nj,converged,emp2',NBI,NBJ,'inf',EMP2_PAIR_TRIPLET_CBS(NBI-NFREEZE,NBJ-NFREEZE,KI,KJ,ISP1,ISP2)
> !         WRITE(*,*)'ni,nj,converged,emp2',NBI,NBJ,'inf',E_SINGLET_CBS+E_TRIPLET_CBS
>          
> 
>       END SUBROUTINE GET_PAIR_CBS_ENERGY_FOR_KQ
> 
> 
> 
>       SUBROUTINE DUMP_CORR_ENERGIES(IO,WDES,WGW,LATT_CUR,ENCUTGW,ENCUTGWSOFT,FSG)
>          USE prec
>          USE poscar
>          USE pseudo
>          USE wave_high
>          USE full_kpoints
>          USE mkpoints
>          USE lattice
>          USE constant
>          USE base
>          IMPLICIT NONE
>          TYPE (in_struct) IO
>          TYPE(wavedes) WDES
>          TYPE(wavedes) WGW
>          TYPE(latt) LATT_CUR
>          REAL(q) :: ENCUTGW,ENCUTGWSOFT
>          REAL(q) :: POTFAK(GRIDHF%MPLWV) 
>          ! local variables
>          TYPE(wavedes1) WGWQ
>          REAL(q) :: FSG ! singularity correction
>          COMPLEX(q) :: ETMP
>          INTEGER :: KQ,NG,KI,KJ,NBI,NBJ,I,ISP2,ISP1, RKI, RKJ
> 
> !=======================================================================================
> !    COMPUTE AND WRITE FINAL MP2 ENERGIES EXTRAPOLATED TO THE CBS LIMIT
> !=======================================================================================
> 
>          CALL M_sum_z(WGW%COMM_INTER, SFACTOR_FULL(1,1), NGVECTOR*REALNKPTS)
>          CALL M_sum_z(WGW%COMM_INTER, SFACTOR_DIRECT(1,1), NGVECTOR*REALNKPTS)
>          CALL M_sum_z(WGW%COMM_INTER, SFACTOR_EXCHANGE(1,1), NGVECTOR*REALNKPTS)
> 
> 
>          energy_c=(0.0_q,0.0_q)
>          energy_x=(0.0_q,0.0_q)
>          E_MP2=(0.0_q,0.0_q)
> 
>          DO KQ=1,WDES%NKPTS
>             IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0)) CYCLE
>             CALL SETWDES(WGW,WGWQ,KQ)
>             IF (ENCUTGW /= ENCUTGWSOFT .AND. ENCUTGW > 0 .AND. ENCUTGWSOFT > 0) THEN
>                CALL SET_GFAC_WAVEFUN(WGWQ,LATT_CUR,FSG,POTFAK,ENCUTGW,ENCUTGWSOFT)
>             ELSE
>                CALL SET_GFAC_WAVEFUN(WGWQ,LATT_CUR,FSG,POTFAK)
>             ENDIF
> 
>             DO NG=1,NGVECTOR
>                energy_c=energy_c+SFACTOR_DIRECT(NG,RKQofKQ(KQ))*POTFAK(NG)
>                energy_x=energy_x+SFACTOR_EXCHANGE(NG,RKQofKQ(KQ))*POTFAK(NG)
>             ENDDO
> 
>          ENDDO
>          
>          E_MP2=energy_c+energy_x
> !         WRITE(*,*)'eMP2=',E_MP2
> 
>          ISP1=1 !fixme for spin-polarized cases
>          ISP2=1
> 
>          DO KI=1,KPOINTS_ORIG%NKPTS
>             IF (WDES%WTKPT(KI)==0) CYCLE
>             IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI)==0)) CYCLE
>             RKI=RKIofKI(KI)
>                
>             DO KJ=1,WDES%NKPTS
>                IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ)==0)) CYCLE
>                RKJ=RKIofKI(KJ)
> 
>                DO NBI=1,VBMAX !loop over valence bands i                     
>                   IF ((LFREEZE) .and. (NBI<=NFREEZE)) CYCLE
>                DO NBJ=1,VBMAX !loop over valence bands i                     
>                   IF ((LFREEZE) .and. (NBJ<=NFREEZE)) CYCLE
> 
> #ifdef gammareal
>          IF ((NBJ>NBI) .and. (WDES%ISPIN==1)) CYCLE
> #else
> 
> #endif
>                   DO I=1,N_MP2_ENCUTS
>                      EMP2_SINGLET(I)=EMP2_SINGLET(I)+EMP2_PAIR_SINGLET_ENCUTS(NBI-NFREEZE,NBJ-NFREEZE,RKI,RKJ,ISP1,ISP2,I)
>                      EMP2_TRIPLET(I)=EMP2_TRIPLET(I)+EMP2_PAIR_TRIPLET_ENCUTS(NBI-NFREEZE,NBJ-NFREEZE,RKI,RKJ,ISP1,ISP2,I)
>                   ENDDO
> 
>                   EMP2_SINGLET_CBS=EMP2_SINGLET_CBS+EMP2_PAIR_SINGLET_CBS(NBI-NFREEZE,NBJ-NFREEZE,RKI,RKJ,ISP1,ISP2)
>                   EMP2_TRIPLET_CBS=EMP2_TRIPLET_CBS+EMP2_PAIR_TRIPLET_CBS(NBI-NFREEZE,NBJ-NFREEZE,RKI,RKJ,ISP1,ISP2)
> 
>                ENDDO
>                ENDDO
>             ENDDO
>          ENDDO
> 
>      ! output
>          IF (IO%IU6>=0) THEN
>             WRITE(IO%IU6,*) " "
>             WRITE(IO%IU6,*) " -------------------------------------------------------------------------"
>             WRITE(IO%IU6,'(5A)') "      ENCUTGW", "  ENCUTGWSOFT", "      singlet MP2", "    triplet MP2", "       total MP2"
>             WRITE(IO%IU6,*) " -------------------------------------------------------------------------"
>             DO I = N_MP2_ENCUTS, 1, -1
>               WRITE(IO%IU6,'(2F13.3,3F16.8)') ENCUTS(I), ENCUTSOFTS(I), REAL(EMP2_SINGLET(I)), REAL(EMP2_TRIPLET(I)), &
>                                         REAL(EMP2_SINGLET(I)+EMP2_TRIPLET(I))
>             ENDDO
>             WRITE(IO%IU6,'(A)') "      linear regression"
>             WRITE(IO%IU6,'(A,3F16.8)') "      converged values    ", REAL(EMP2_SINGLET_CBS), REAL(EMP2_TRIPLET_CBS),REAL(EMP2_SINGLET_CBS+EMP2_TRIPLET_CBS)
>          ENDIF
>          IF (IO%IU0>=0) THEN
>             WRITE(IO%IU0,*) " "
>             WRITE(IO%IU0,*) " -------------------------------------------------------------------------"
>             WRITE(IO%IU0,'(5A)') "      ENCUTGW", "  ENCUTGWSOFT", "      singlet MP2", "    triplet MP2", "       total MP2"
>             WRITE(IO%IU0,*) " -------------------------------------------------------------------------"
>             DO I = N_MP2_ENCUTS, 1, -1
>               WRITE(IO%IU0,'(2F13.3,3F16.8)') ENCUTS(I), ENCUTSOFTS(I), REAL(EMP2_SINGLET(I)), REAL(EMP2_TRIPLET(I)), &
>                                         REAL(EMP2_SINGLET(I)+EMP2_TRIPLET(I))
>             ENDDO
>             WRITE(IO%IU0,'(A)') "      linear regression"
>             WRITE(IO%IU0,'(A,3F16.8)') "      converged values    ", REAL(EMP2_SINGLET_CBS), REAL(EMP2_TRIPLET_CBS), REAL(EMP2_SINGLET_CBS+EMP2_TRIPLET_CBS)
>          ENDIF
> 
> !=======================================================================================
> !    DUMP MP2 PAIR ENERGIES FOR CC4S
> !=======================================================================================
> 
>          IF (ME==0) write(*,*)''
>          IF (ME==0) write(*,*)'Writing Mp2PairEnergies .'
> 
>          IF (ME==0)  THEN
> 
>             OPEN(unit = 7,file = "Mp2PairEnergies.yaml",FORM='FORMATTED',access='stream',STATUS='REPLACE')
> 
>             WRITE(7,'(A)')    'version: 100'
>             WRITE(7,'(A)')    'type: Tensor'
>             WRITE(7,'(A)') 'scalarType: Real64'
>             WRITE(7,'(A)')    'dimensions:'
>             WRITE(7,'(A,I6)') '- length: ',((VBMAX-NFREEZE)*REALNKPTS*WDES%ISPIN)
>             WRITE(7,'(A)')    '  type: State'
>             WRITE(7,'(A,I6)') '- length: ',((VBMAX-NFREEZE)*REALNKPTS*WDES%ISPIN)
>             WRITE(7,'(A)')    '  type: State'
>             WRITE(7,'(A)')    'elements:'
>             WRITE(7,'(A)')    '  type: TextFile'
>             WRITE(7,'(A)') 'unit: 0.03674932217563878       # = (Eh/eV)'
>             CLOSE(7)
> 
>             OPEN(unit = 7,file = "Mp2PairEnergies.elements",FORM='FORMATTED',access='stream',STATUS='REPLACE')
> 
> 
>             ISP1=1
>             ISP2=1
> 
> 
>             DO KI=1,KPOINTS_ORIG%NKPTS
>                IF (WDES%WTKPT(KI)==0) CYCLE
>                IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI)==0)) CYCLE
>                RKI=RKIofKI(KI)
>                
>                DO KJ=1,WDES%NKPTS
>                   IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ)==0)) CYCLE
>                   RKJ=RKIofKI(KJ)
>                
>                   DO NBI=1,VBMAX !loop over valence bands i                     
>                      IF ((LFREEZE) .and. (NBI<=NFREEZE)) CYCLE
>                   DO NBJ=1,VBMAX !loop over valence bands i                     
>                      IF ((LFREEZE) .and. (NBJ<=NFREEZE)) CYCLE
> 
> #ifdef gammareal
>                      IF ((NBJ<NBI) .and. (WDES%ISPIN==1)) THEN 
>                         ETMP=EMP2_PAIR_SINGLET_CBS(NBI-NFREEZE,NBJ-NFREEZE,KI,KJ,ISP1,ISP2)/2.0_q
>                         ETMP=ETMP+EMP2_PAIR_TRIPLET_CBS(NBI-NFREEZE,NBJ-NFREEZE,KI,KJ,ISP1,ISP2)/2.0_q
>                      ENDIF
> 
>                      IF ((NBJ>NBI) .and. (WDES%ISPIN==1)) THEN
>                         ETMP=EMP2_PAIR_SINGLET_CBS(NBJ-NFREEZE,NBI-NFREEZE,KJ,KI,ISP2,ISP1)/2.0_q
>                         ETMP=ETMP+EMP2_PAIR_TRIPLET_CBS(NBJ-NFREEZE,NBI-NFREEZE,KJ,KI,ISP2,ISP1)/2.0_q
>                      ENDIF
> 
>                      IF ((NBJ==NBI) .and. (WDES%ISPIN==1)) THEN
>                         ETMP=EMP2_PAIR_SINGLET_CBS(NBJ-NFREEZE,NBI-NFREEZE,KJ,KI,ISP2,ISP1)
>                         ETMP=ETMP+EMP2_PAIR_TRIPLET_CBS(NBJ-NFREEZE,NBI-NFREEZE,KJ,KI,ISP2,ISP1)
>                      ENDIF
> #else
>                      ETMP=EMP2_PAIR_SINGLET_CBS(NBI-NFREEZE,NBJ-NFREEZE,RKI,RKJ,ISP1,ISP2)
>                      ETMP=ETMP+EMP2_PAIR_TRIPLET_CBS(NBI-NFREEZE,NBJ-NFREEZE,RKI,RKJ,ISP1,ISP2)
> #endif
>                      IF (ME==0) WRITE(7,*) REAL(ETMP,KIND=8)
> 
>                   ENDDO
>                   ENDDO
>                ENDDO
>             ENDDO
> 
> 
> 
>             CLOSE(7)
>          ENDIF
> 
> 
> 
>       END SUBROUTINE DUMP_CORR_ENERGIES
> 
1060a1700
>          REAL(q) :: POTFAK_DUMMY(GRIDHF%MPLWV) 
1092a1733
> 
1093a1735,1737
>          DO KI=1,WDES%NKPTS
>             VBMAX=MAX(LAST_FILLED_XI_NOMOD(W,KI,1),VBMAX)
>          ENDDO
1102,1104c1746,1758
<          IF (LDUMPPAIRS) THEN
<             ALLOCATE(EMP2_PAIR(VBMAX,VBMAX,WDES%ISPIN,WDES%ISPIN))
<             EMP2_PAIR=(0.0_q,0.0_q)
---
>          IF (LSFACTOR) THEN
>             ALLOCATE(EMP2_PAIR_SINGLET_ENCUTS(VBMAX,VBMAX,REALNKPTS,REALNKPTS,WDES%ISPIN,WDES%ISPIN,N_MP2_ENCUTS))
>             ALLOCATE(EMP2_PAIR_TRIPLET_ENCUTS(VBMAX,VBMAX,REALNKPTS,REALNKPTS,WDES%ISPIN,WDES%ISPIN,N_MP2_ENCUTS))
>             ALLOCATE(EMP2_PAIR_SINGLET_CBS(VBMAX,VBMAX,REALNKPTS,REALNKPTS,WDES%ISPIN,WDES%ISPIN))
>             ALLOCATE(EMP2_PAIR_TRIPLET_CBS(VBMAX,VBMAX,REALNKPTS,REALNKPTS,WDES%ISPIN,WDES%ISPIN))
>             ALLOCATE(EMP2_SINGLET(N_MP2_ENCUTS))
>             ALLOCATE(EMP2_TRIPLET(N_MP2_ENCUTS))
>             EMP2_PAIR_SINGLET_ENCUTS=(0.0_q,0.0_q)
>             EMP2_PAIR_TRIPLET_ENCUTS=(0.0_q,0.0_q)
>             EMP2_PAIR_SINGLET_CBS=(0.0_q,0.0_q)
>             EMP2_PAIR_TRIPLET_CBS=(0.0_q,0.0_q)
>             EMP2_SINGLET=(0.0_q,0.0_q)
>             EMP2_TRIPLET=(0.0_q,0.0_q)
1164a1819,1820
>          ALLOCATE(POTFAK_FULL(NGVECTOR,REALNKPTS))
>          POTFAK_FULL=(0.0_q,0.0_q)
1238c1894
<                   POTFAK=SQRT(POTFAK)
---
> !                  POTFAK=SQRT(POTFAK)
1239a1896
>                POTFAK_FULL(1:NP,RKQofKQ(KQ))=POTFAK(1:NP)
1283c1940,1941
<                            !  include multipole correction
---
>                            !  include multipole correction with dummy POTFAK here. We add POTFAK later.
>                            POTFAK_DUMMY=1.0_q
1285,1288c1943
<                             POTFAK(1))
<                         ELSE
<                            CALL APPLY_GFAC_WAVEFUN(WGWQ,GCHGIA(1,NBAA,1), &
<                             POTFAK(1))
---
>                             POTFAK_DUMMY(1))
1305c1960,1962
< #ifdef gammareal
---
> ! #ifdef gammareal
> !                           tmp_FTOD_PW(1:NP,NBAA+NBA-1,NBI-NFREEZE,1)=(GCHGIA(1:NP,NBAA,1))*SQRT(1.0_q/GRIDHF%NPLWV)
> ! #else
1307,1309c1964,1965
< #else
<                            tmp_FTOD_PW(1:NP,NBAA+NBA-1,NBI-NFREEZE,1)=(GCHGIA(1:NP,NBAA,1))
< #endif
---
> 
> ! #endif
1410c2066
<                         tmp_FTOD_PW(1:NP,NBAA+NBA-1,NBI-NFREEZE,2)=(GCHGIA(1:NP,NBAA,2)*(1.0_q/GRIDHF%NPLWV))
---
>                         tmp_FTOD_PW(1:NP,NBAA+NBA-1,NBI-NFREEZE,2)=GCHGIA(1:NP,NBAA,2)*SQRT(1.0_q/GRIDHF%NPLWV)
1561c2217
<          MB=MIN(((NGVECTOR)/NPROW),50)   !Row blocking size
---
>          MB=MIN(((NGVECTOR)/NPROW),BLOCKSIZE)   !Row blocking size
1672a2329,2332
>          IF (ALLOCATED(TWOE4ORBITAL_FULL)) THEN
>             DEALLOCATE(TWOE4ORBITAL_FULL)
>          ENDIF
> 
1674a2335
>          IF (LSFACTOR) allocate(TWOE4ORBITAL_FULL(TWOE4ORBITAL_ROWS,TWOE4ORBITAL_COLS))
1692c2353
<          MB=MIN(((NGVECTOR)/NPROW),50)   !Row blocking size
---
>          MB=MIN(((NGVECTOR)/NPROW),BLOCKSIZE)   !Row blocking size
1710a2372
> !         WRITE(*,*)'ME, FTOD rows',ME,FTOD_PW_rows
1711a2374
>          ALLOCATE(FTOD_PW_POTFAK(FTOD_PW_rows,FTOD_PW_cols,3))
1715c2378
<             MB=MIN(((NHVECTOR)/NPROW),50)   !Row blocking size
---
>             MB=MIN(((NHVECTOR)/NPROW),BLOCKSIZE)   !Row blocking size
1738a2402,2461
>       SUBROUTINE SETUP_POTFAK_FULL_DIST(WDES)
>          IMPLICIT NONE
>          TYPE(wavedes) WDES              
>          INTEGER :: FTOD_PW_rows
>          INTEGER :: NG,RNG,KQ
> 
>          call BLACS_GRIDINFO(CONTXT_GRID, NPROW, NPCOL, MYROW, MYCOL)
>          !Blocking size for block-cyclic distribution of FTOD_PW        
>          MB=MIN(((NGVECTOR)/NPROW),BLOCKSIZE)   !Row blocking size
>          
>          FTOD_PW_rows = numroc(NGVECTOR,mb,MYROW,0,NPROW)
> !         WRITE(*,*)'potfak ME, FTOD rows',ME,FTOD_PW_rows
>          
>          ALLOCATE(POTFAK_FULL_DIST(FTOD_PW_rows,REALNKPTS))
> 
>          DO NG=1,FTOD_PW_rows
>             CALL LOC2GLOB(NG,MYROW,NGVECTOR,NPROW,MB,RNG)     
>             DO KQ=1,WDES%NKPTS
>                IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0)) CYCLE
>                POTFAK_FULL_DIST(NG,RKQofKQ(KQ))=POTFAK_FULL(RNG,RKQofKQ(KQ))
>             ENDDO
>          ENDDO
> 
>       END SUBROUTINE SETUP_POTFAK_FULL_DIST
>  
>       SUBROUTINE SETUP_ENCUTS(ENCUTGW,ENCUTGWSOFT)
>          IMPLICIT NONE
>          REAL(q) :: ENCUTGW,ENCUTGWSOFT
>          INTEGER :: I
> 
>          ! allocate memory for the different ENCUTs
>          ALLOCATE(ENCUTS(N_MP2_ENCUTS))
>          ALLOCATE(ENCUTSOFTS(N_MP2_ENCUTS))
>          ! calculate different ENCUTs
>          DO I = 1, N_MP2_ENCUTS
>             ENCUTS(I) = ENCUTGW / ( ENCUT_DIVISOR**(I-1) )
>             ENCUTSOFTS(I) = ENCUTGWSOFT / ( ENCUT_DIVISOR**(I-1) )
>          ENDDO
>       END SUBROUTINE SETUP_ENCUTS
> 
>       SUBROUTINE SETUP_SFACTOR(WDES)
>          IMPLICIT NONE
>          TYPE(wavedes) WDES              
> 
>          ALLOCATE(SFACTOR_FULL(NGVECTOR,REALNKPTS))
>          ALLOCATE(SFACTOR_DIRECT(NGVECTOR,REALNKPTS))
>          ALLOCATE(SFACTOR_EXCHANGE(NGVECTOR,REALNKPTS))
>          ALLOCATE(SFACTOR_DIRECT_PAIR(NGVECTOR))
>          ALLOCATE(SFACTOR_EXCHANGE_PAIR(NGVECTOR))
> 
>          SFACTOR_FULL=(0.0_q,0.0_q)
>          SFACTOR_DIRECT=(0.0_q,0.0_q)
>          SFACTOR_EXCHANGE=(0.0_q,0.0_q)
>          SFACTOR_DIRECT_PAIR=(0.0_q,0.0_q)
>          SFACTOR_EXCHANGE_PAIR=(0.0_q,0.0_q)
> 
>       END SUBROUTINE SETUP_SFACTOR
>  
> 
> 
1859a2583,2656
>   !****************************** LIN_REG_CONVERGE_LTMP2 *****************************
>   !
>   ! Perform a linear regression of the input data: y = M*x + B, where 
>   ! x = ENCUTS^(-3/2) and y = ENERGIES. 
>   !
>   ! This is used to extrapolate ENCUT -> infinity
>   ! Because: x(ENCUT -> infinity) = 0 => y -> b = converged energy
>   !
>   ! see: http://dx.doi.org/10.1103/PhysRevB.90.075125 for singlet, triplets need
>   ! different extrapolation
>   !
>   !***********************************************************************************
>   SUBROUTINE LIN_REG_CONVERGE_SINGLET(ENERGIES, ENCUTS, CONVERGED_ENERGY)
>     REAL(q) :: ENERGIES(:)         ! calculated energies, each for a given ENCUT
>     REAL(q) :: ENCUTS(:)           ! all ENCUTs
>     REAL(q) :: CONVERGED_ENERGY   ! the resulting converged energy
>     ! local
>     INTEGER :: I
>     INTEGER :: N                            ! number of sample points
>     REAL(q) :: SUM_X, SUM_X2, SUM_Y, SUM_XY ! sums of the samples
>     REAL(q) :: M, B                         ! parameters of y = M*x + B
> 
>     N = SIZE(ENERGIES)
>     SUM_X  = 0
>     SUM_X2 = 0
>     SUM_Y  = 0
>     SUM_XY = 0
> 
>     ! perform linear regression
>     DO I = 1, N
>       SUM_X  = SUM_X  + ENCUTS(I)**(-3.0_q/2)
>       SUM_X2 = SUM_X2 + ENCUTS(I)**(-3.0_q)    ! ENCUTS^3 corresponds to x^2
>       SUM_Y  = SUM_Y  + ENERGIES(I)
>       SUM_XY = SUM_XY + ENCUTS(I)**(-3.0_q/2) * ENERGIES(I)
>     ENDDO
>     M = ( N*SUM_XY - SUM_X*SUM_Y) / ( N*SUM_X2 - SUM_X**2 )
>     B = ( SUM_Y - M*SUM_X ) / N
> 
>     CONVERGED_ENERGY = B
> 
>   END SUBROUTINE LIN_REG_CONVERGE_SINGLET
> 
>   SUBROUTINE LIN_REG_CONVERGE_TRIPLET(ENERGIES, ENCUTS, CONVERGED_ENERGY)
>     REAL(q) :: ENERGIES(:)         ! calculated energies, each for a given ENCUT
>     REAL(q) :: ENCUTS(:)           ! all ENCUTs
>     REAL(q) :: CONVERGED_ENERGY   ! the resulting converged energy
>     ! local
>     INTEGER :: I
>     INTEGER :: N                            ! number of sample points
>     REAL(q) :: SUM_X, SUM_X2, SUM_Y, SUM_XY ! sums of the samples
>     REAL(q) :: M, B                         ! parameters of y = M*x + B
> 
>     N = SIZE(ENERGIES)
>     SUM_X  = 0
>     SUM_X2 = 0
>     SUM_Y  = 0
>     SUM_XY = 0
> 
>     ! perform linear regression
>     DO I = 1, N
>       SUM_X  = SUM_X  + ENCUTS(I)**(-3.0_q/2) !fix me
>       SUM_X2 = SUM_X2 + ENCUTS(I)**(-3.0_q)    ! ENCUTS^3 corresponds to x^2
>       SUM_Y  = SUM_Y  + ENERGIES(I)
>       SUM_XY = SUM_XY + ENCUTS(I)**(-3.0_q/2) * ENERGIES(I)
>     ENDDO
>     M = ( N*SUM_XY - SUM_X*SUM_Y) / ( N*SUM_X2 - SUM_X**2 )
>     B = ( SUM_Y - M*SUM_X ) / N
> 
>     CONVERGED_ENERGY = B
> 
>   END SUBROUTINE LIN_REG_CONVERGE_TRIPLET
> 
> 
> 
1863c2660
<       SUBROUTINE MP2_READER(IU5,IU6,IU0)
---
>    SUBROUTINE MP2_READER(IU5,IU6,IU0)
1880c2677,2683
<       CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBANDSHIGH', NBANDSHIGH, IERR, WRITEXMLINCAR)
---
>       CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBANDSHIGH', NBANDSHIGH_, IERR, WRITEXMLINCAR)
> 
> !      LDUMPPAIRS=.FALSE.
> !      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LDUMPPAIRS', LDUMPPAIRS, IERR, WRITEXMLINCAR)
> 
>       LSFACTOR=.FALSE.
>       CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSFACTOR', LSFACTOR, IERR, WRITEXMLINCAR)
1882,1883d2684
<       LDUMPPAIRS=.FALSE.
<       CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LDUMPPAIRS', LDUMPPAIRS, IERR, WRITEXMLINCAR)
1886c2687
<       END SUBROUTINE MP2_READER
---
>    END SUBROUTINE MP2_READER
EOF

md5sum_orig=55a341fc0bb081834d9fa1f33a9294ec
md5sum_new=`md5sum ump2.F  | awk '{print $1}' `
if [ "$md5sum_orig" = "$md5sum_new" ]; then
    echo "Patching file  ump2.F " 
    patch ump2.F .tmp_vaspcc4s_patch 
else
   echo "You dont have the correct version of file  ump2.F  to apply patch. " 
fi
 
#######################
# diff of   ump2no.F
#######################
 
 
cat > .tmp_vaspcc4s_patch<<"EOF"
11a12,13
> ! AUTHOR :: Andreas Grüneis ( andreas.grueneis@tuwien.ac.at )
> !
25c27,35
< ! aG 2012
---
> ! ================================================
> ! THE FOLLOWING SETTINGS ARE SUPPORTED:
> ! ================================================
> !
> !   * RESTRICTED (RHF) AND UNRESTRICTED (UHF) CANONICAL ORBITALS ONLY
> !   * GAMMA-CENTERED OR SHIFTED K-MESHES
> !
> ! see https://doi.org/10.1021/ct200263g for more details
> !
59a70
>       LOGICAL :: LOVERWRITEVBMAX, LMETAL
68c79
<       LOGICAL :: OEAPPROXIMATE_NOs
---
>       LOGICAL :: LLOWMEM
69a81,83
>       LOGICAL :: LFREEZE
>       INTEGER :: NFREEZE
>       INTEGER :: NBANDSHIGH_
70a85,91
> !!!!!
>       LOGICAL :: SHIFTED_KPOINTS
>       REAL(q) :: WTKPT
>       INTEGER :: REALNKPTS
>       INTEGER, ALLOCATABLE :: RKQofKQ(:), RKIofKI(:)
>       INTEGER, ALLOCATABLE :: KQofRKQ(:), KIofRKI(:)
> 
107c128
<          INTEGER KI,KQ,KQ_,KA,KB,ISP,KJ,NBI,NBJ,ISP2,PINFO,N
---
>          INTEGER KI,KQ,KQ_,KA,KB,ISP,KJ,NBI,NBJ,ISP2,PINFO,N, MNBI,MNBJ
117c138,139
<          OEAPPROXIMATE_NOs=.FALSE.
---
>          IF ((IO%IU0>0)) write(IO%IU0,*)'MP2 NO routine starting.'
>          LLOWMEM=.FALSE.
118a141,146
>          LOVERWRITEVBMAX=.FALSE.
> 
>          NFREEZE=0
>          
>          CALL NFREEZE_READER(IO%IU5, IO%IU6, IO%IU0)
>          IF (NFREEZE/=0) LFREEZE=.TRUE.
129a158,166
>          CALL CHECK_SHIFTED_KPOINTS(WDES,W,KPOINTS) 
>          IF ((IO%IU0>0) .and. (SHIFTED_KPOINTS)) write(IO%IU0,*)'You are using a shifted k-mesh.'
> 
>          NBANDSHIGH_=WDES%NB_TOT
>          CALL NOINCAR_READER(IO%IU5, IO%IU6, IO%IU0)
>          IF ((IO%IU0>0) .and. (APPROXIMATE_NOs)) write(IO%IU0,*)'APPROXIMATE NOs will be calculated.'
>          IF ((APPROXIMATE_NOs) .and. (REALNKPTS==1)) LLOWMEM=.TRUE.
>          IF ((APPROXIMATE_NOs) .and. (REALNKPTS>1)) LLOWMEM=.FALSE.   !low memory variant not tested for more than one k-point
>          
135c172
<          ALLOCATE(CELTOT_NEW(WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN))
---
>          ALLOCATE(CELTOT_NEW(WDES%NB_TOT,REALNKPTS,WDES%ISPIN))
137,138c174
<          
<          !Initialize the 1D process grid
---
>         !Initialize the 1D process grid
146,147c182,187
<          !Calculate the FTOD functions
<          CALL CALC_2ORBITAL_FTOD(WDES,WGW,W,P,T_INFO,LATT_CUR,LMAXMP2,ENCUTGW,ENCUTGWSOFT,IO,FSG_STORE(1))         
---
>          IF (LLOWMEM) THEN
>             CALL CALC_2ORBITAL_FTOD_LOWMEM(WDES,WGW,W,P,T_INFO,INFO,LATT_CUR,LMAXMP2,ENCUTGW,ENCUTGWSOFT,IO,FSG_STORE(1),0,0,0,0)
>          ELSE
>             !Calculate the FTOD functions
>             CALL CALC_2ORBITAL_FTOD(WDES,WGW,W,P,T_INFO,INFO,LATT_CUR,LMAXMP2,ENCUTGW,ENCUTGWSOFT,IO,FSG_STORE(1))         
>          ENDIF
158c198,202
<             WRITE(IO%IU0,*)'Calculating MP2 natural orbitals:'
---
>             IF (LLOWMEM) THEN
>                WRITE(IO%IU0,*)'Calculating approximate MP2 natural orbitals:'
>             ELSE
>                WRITE(IO%IU0,*)'Calculating MP2 natural orbitals:'
>             ENDIF
161a206
>             IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(RKA)==0)) CYCLE
177a223
>                IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0)) CYCLE
179c225
<                CALL GWPROGRESS(IO%IU0, KI,KPOINTS_ORIG%NKPTS,KJ,WDES%NKPTS)
---
>                IF (.not. LLOWMEM) CALL GWPROGRESS(IO%IU0, KI,KPOINTS_ORIG%NKPTS,KJ,WDES%NKPTS)
180a227
>                   IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI)==0)) CYCLE
187a235
>                      IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ)==0)) CYCLE
189c237
<                      IF ((OEAPPROXIMATE_NOs) .and. (KI/=KJ)) CYCLE
---
>                      IF ((APPROXIMATE_NOs) .and. (KI/=KJ)) CYCLE
197,203c245,257
<                      DO NBI=1,VBMAX(ISP) !loop over valence bands i                     
<                         IF ((OEAPPROXIMATE_NOs) .and. (NBI>1)) CYCLE
<                      
<                         DO NBJ=1,VBMAX(ISP) !loop over valence bands j
<                            IF ((OEAPPROXIMATE_NOs) .and. (NBJ>1)) CYCLE
<                            IF ((APPROXIMATE_NOs) .and. (WDES%ISPIN==2)) CYCLE
<                            IF ((APPROXIMATE_NOs) .and. (NBI/=NBJ)) CYCLE
---
>                      DO MNBI=1,VBMAX(ISP) !loop over valence bands i
>                         NBI=MNBI
>                         IF ((LFREEZE) .and. (NFREEZE>=MNBI))  CYCLE
>                         IF (LFREEZE) NBI=MNBI-NFREEZE
>                         
>                     
>                         DO MNBJ=1,VBMAX(ISP) !loop over valence bands j
>                            NBJ=MNBJ
> !                           IF ((APPROXIMATE_NOs) .and. (WDES%ISPIN==2)) CYCLE
>                            IF ((APPROXIMATE_NOs) .and. (MNBI/=MNBJ)) CYCLE
>                            IF ((LFREEZE) .and. (NFREEZE>=MNBJ))  CYCLE
>                            IF (LFREEZE) NBJ=MNBJ-NFREEZE
>                            
207a262,267
> 
>                            IF (LLOWMEM) THEN
> 
>                               CALL CALC_2ORBITAL_FTOD_LOWMEM(WDES,WGW,W,P,T_INFO,INFO,LATT_CUR, &
>                                                  LMAXMP2,ENCUTGW,ENCUTGWSOFT,IO,FSG_STORE(1),MNBI,KI,KQ,ISP)
>  
209,212c269,272
<                            CALL PDGEMM(trans,'n',(PROCS*WDES%NBANDS),(PROCS*WDES%NBANDS),&
<                                   m_ NGVECTOR,-one, FTOD_PW(1,1,NBI,KI,KQ,ISP,1),1,1,&
<                                   desc_FTOD_PW,FTOD_PW(1,1,NBJ,KJ,KQ,ISP,ncc),1,1,&
<                                   desc_FTOD_PW,zero, TWOE4ORBITAL(1,1),1,1,desc_TWOE4ORBITAL)
---
>                               CALL PDGEMM(trans,'n',(PROCS*WDES%NBANDS),(PROCS*WDES%NBANDS),&
>                                      m_ NGVECTOR,-one, FTOD_PW(1,1,1,1,1,1,1),1,1,&
>                                      desc_FTOD_PW,FTOD_PW(1,1,1,1,1,1,ncc),1,1,&
>                                      desc_FTOD_PW,zero, TWOE4ORBITAL(1,1),1,1,desc_TWOE4ORBITAL)
214,217c274,277
<                            CALL PZGEMM(trans,'n',(PROCS*WDES%NBANDS),(PROCS*WDES%NBANDS),&
<                                   NGVECTOR,-one, FTOD_PW(1,1,NBI,KI,KQ,ISP,1),1,1,&
<                                   desc_FTOD_PW,FTOD_PW(1,1,NBJ,KJ,KQ,ISP,ncc),1,1,&
<                                   desc_FTOD_PW,zero, TWOE4ORBITAL(1,1),1,1,desc_TWOE4ORBITAL)
---
>                               CALL PZGEMM(trans,'n',(PROCS*WDES%NBANDS),(PROCS*WDES%NBANDS),&
>                                      NGVECTOR,-one, FTOD_PW(1,1,1,1,1,1,1),1,1,&
>                                      desc_FTOD_PW,FTOD_PW(1,1,1,1,1,1,ncc),1,1,&
>                                      desc_FTOD_PW,zero, TWOE4ORBITAL(1,1),1,1,desc_TWOE4ORBITAL)
219c279,293
<                            IF (ASSOCIATED(H)) THEN
---
>                               IF (ASSOCIATED(H)) THEN
> #ifdef gammareal
>                                  CALL PDGEMM(trans,'n',(PROCS*WDES%NBANDS),(PROCS*WDES%NBANDS),&
>                                         NHVECTOR,-one, FTOD_OC(1,1,1,1,1,1,1),1,1,&
>                                         desc_FTOD_OC,FTOD_OC(1,1,1,1,1,1,2),1,1,&
>                                         desc_FTOD_OC,one, TWOE4ORBITAL(1,1),1,1,desc_TWOE4ORBITAL)
> #else
>                                  CALL PZGEMM(trans,'n',(PROCS*WDES%NBANDS),(PROCS*WDES%NBANDS),&
>                                         NHVECTOR,-one, FTOD_OC(1,1,1,1,1,1,1),1,1,&
>                                         desc_FTOD_OC,FTOD_OC(1,1,1,1,1,ISP,2),1,1,&
>                                         desc_FTOD_OC,one, TWOE4ORBITAL(1,1),1,1,desc_TWOE4ORBITAL)
> #endif
>                               ENDIF
> 
>                            ELSE
222,224c296,298
<                                      NHVECTOR,-one, FTOD_OC(1,1,NBI,KI,KQ,ISP,1),1,1,&
<                                      desc_FTOD_OC,FTOD_OC(1,1,NBJ,KJ,KQ,ISP,2),1,1,&
<                                      desc_FTOD_OC,one, TWOE4ORBITAL(1,1),1,1,desc_TWOE4ORBITAL)
---
>                                      m_ NGVECTOR,-one, FTOD_PW(1,1,NBI,KI,KQ,ISP,1),1,1,&
>                                      desc_FTOD_PW,FTOD_PW(1,1,NBJ,KJ,KQ,ISP,ncc),1,1,&
>                                      desc_FTOD_PW,zero, TWOE4ORBITAL(1,1),1,1,desc_TWOE4ORBITAL)
227,229c301,303
<                                      NHVECTOR,-one, FTOD_OC(1,1,NBI,KI,KQ,ISP,1),1,1,&
<                                      desc_FTOD_OC,FTOD_OC(1,1,NBJ,KJ,KQ,ISP,2),1,1,&
<                                      desc_FTOD_OC,one, TWOE4ORBITAL(1,1),1,1,desc_TWOE4ORBITAL)
---
>                                      NGVECTOR,-one, FTOD_PW(1,1,NBI,RKIofKI(KI),RKQofKQ(KQ),ISP,1),1,1,&
>                                      desc_FTOD_PW,FTOD_PW(1,1,NBJ,RKIofKI(KJ),RKQofKQ(KQ),ISP,ncc),1,1,&
>                                      desc_FTOD_PW,zero, TWOE4ORBITAL(1,1),1,1,desc_TWOE4ORBITAL)
230a305,317
>                               IF (ASSOCIATED(H)) THEN
> #ifdef gammareal
>                                  CALL PDGEMM(trans,'n',(PROCS*WDES%NBANDS),(PROCS*WDES%NBANDS),&
>                                         NHVECTOR,-one, FTOD_OC(1,1,NBI,KI,KQ,ISP,1),1,1,&
>                                         desc_FTOD_OC,FTOD_OC(1,1,NBJ,KJ,KQ,ISP,2),1,1,&
>                                         desc_FTOD_OC,one, TWOE4ORBITAL(1,1),1,1,desc_TWOE4ORBITAL)
> #else
>                                  CALL PZGEMM(trans,'n',(PROCS*WDES%NBANDS),(PROCS*WDES%NBANDS),&
>                                         NHVECTOR,-one, FTOD_OC(1,1,NBI,RKIofKI(KI),RKQofKQ(KQ),ISP,1),1,1,&
>                                         desc_FTOD_OC,FTOD_OC(1,1,NBJ,RKIofKI(KJ),RKQofKQ(KQ),ISP,2),1,1,&
>                                         desc_FTOD_OC,one, TWOE4ORBITAL(1,1),1,1,desc_TWOE4ORBITAL)
> #endif
>                               ENDIF
254a342,377
> 
>                            IF (LLOWMEM) THEN
>                               CYCLE
>                               CALL CALC_2ORBITAL_FTOD_LOWMEM(WDES,WGW,W,P,T_INFO,INFO,LATT_CUR, &
>                                                  LMAXMP2,ENCUTGW,ENCUTGWSOFT,IO,FSG_STORE(1),MNBI,KI,KQ_,ISP)
> 
> #ifdef gammareal
>                               CALL PDGEMM(trans,'n',(PROCS*WDES%NBANDS),&
>                                      (PROCS*WDES%NBANDS),&
>                                      m_ NGVECTOR,-one, FTOD_PW(1,1,1,1,1,1,ncc),1,1, &
>                                      desc_FTOD_PW,FTOD_PW(1,1,1,1,1,1,1),1,1,&
>                                      desc_FTOD_PW,zero, TWOE4ORBITAL_X(1,1),1,1,&
>                                      desc_TWOE4ORBITAL_X)
> #else
>                               CALL PZGEMM(trans,'n',(PROCS*WDES%NBANDS),&
>                                      (PROCS*WDES%NBANDS),&
>                                      NGVECTOR,-one, FTOD_PW(1,1,1,1,1,1,ncc),1,1, &
>                                      desc_FTOD_PW,FTOD_PW(1,1,1,1,1,1,1),1,1,&
>                                      desc_FTOD_PW,zero, TWOE4ORBITAL_X(1,1),1,1,&
>                                      desc_TWOE4ORBITAL_X)
> #endif
>                               IF (ASSOCIATED(H)) THEN
> #ifdef gammareal
>                                  CALL PDGEMM(trans,'n',(PROCS*WDES%NBANDS),(PROCS*WDES%NBANDS),&
>                                         NHVECTOR,-one, FTOD_OC(1,1,1,1,1,1,2),1,1,&
>                                         desc_FTOD_OC,FTOD_OC(1,1,1,1,1,1,1),1,1,&
>                                         desc_FTOD_OC,one, TWOE4ORBITAL_X(1,1),1,1,desc_TWOE4ORBITAL_X)
> #else
>                                  CALL PZGEMM(trans,'n',(PROCS*WDES%NBANDS),(PROCS*WDES%NBANDS),&
>                                         NHVECTOR,-one, FTOD_OC(1,1,1,1,1,1,2),1,1,&
>                                         desc_FTOD_OC,FTOD_OC(1,1,1,1,1,1,1),1,1,&
>                                         desc_FTOD_OC,one, TWOE4ORBITAL_X(1,1),1,1,desc_TWOE4ORBITAL_X)
> #endif
>                               ENDIF
> 
>                            ELSE
265,266c388,389
<                                      NGVECTOR,-one, FTOD_PW(1,1,NBJ,KJ,KQ_,ISP,ncc),1,1, &
<                                      desc_FTOD_PW,FTOD_PW(1,1,NBI,KI,KQ_,ISP,1),1,1,&
---
>                                      NGVECTOR,-one, FTOD_PW(1,1,NBJ,RKIofKI(KJ),RKQofKQ(KQ_),ISP,ncc),1,1, &
>                                      desc_FTOD_PW,FTOD_PW(1,1,NBI,RKIofKI(KI),RKQofKQ(KQ_),ISP,1),1,1,&
278,279c401,402
<                                         NHVECTOR,-one, FTOD_OC(1,1,NBJ,KJ,KQ_,ISP,2),1,1,&
<                                         desc_FTOD_OC,FTOD_OC(1,1,NBI,KI,KQ_,ISP,1),1,1,&
---
>                                         NHVECTOR,-one, FTOD_OC(1,1,NBJ,RKIofKI(KJ),RKQofKQ(KQ_),ISP,2),1,1,&
>                                         desc_FTOD_OC,FTOD_OC(1,1,NBI,RKIofKI(KI),RKQofKQ(KQ_),ISP,1),1,1,&
282a406
>                            ENDIF
289c413
<                            CALL APPLY_DENOM(W,KI,KJ,KA,KB,NBI,NBJ,ISP,ISP,LATT_CUR,NBANDSGW)
---
>                            CALL APPLY_DENOM(W,KI,KJ,KA,KB,MNBI,MNBJ,ISP,ISP,LATT_CUR,NBANDSGW)
307a432
>                               IF ((APPROXIMATE_NOs) .and. (ISP/=ISP2)) CYCLE
310,312c435,439
<                               DO NBJ=1,VBMAX(ISP2) !loop over valence bands j
<                                  IF ((OEAPPROXIMATE_NOs) .and. (NBJ>1)) CYCLE
<                                  IF ((APPROXIMATE_NOs) .and. (NBI/=NBJ)) CYCLE
---
>                               DO MNBJ=1,VBMAX(ISP2) !loop over valence bands j
>                                  NBJ=MNBJ
>                                  IF ((APPROXIMATE_NOs) .and. (MNBI/=MNBJ)) CYCLE
>                                  IF ((LFREEZE) .and. (NFREEZE>=MNBJ))  CYCLE
>                                  IF (LFREEZE) NBJ=MNBJ-NFREEZE
323,324c450,451
<                                         NGVECTOR,-one, FTOD_PW(1,1,NBI,KI,KQ,ISP,1),1,1,&
<                                         desc_FTOD_PW,FTOD_PW(1,1,NBJ,KJ,KQ,ISP2,ncc),1,1,&
---
>                                         NGVECTOR,-one, FTOD_PW(1,1,NBI,RKIofKI(KI),RKQofKQ(KQ),ISP,1),1,1,&
>                                         desc_FTOD_PW,FTOD_PW(1,1,NBJ,RKIofKI(KJ),RKQofKQ(KQ),ISP2,ncc),1,1,&
335,336c462,463
<                                            NHVECTOR,-one, FTOD_OC(1,1,NBI,KI,KQ,ISP,1),1,1,&
<                                            desc_FTOD_OC,FTOD_OC(1,1,NBJ,KJ,KQ,ISP2,2),1,1,&
---
>                                            NHVECTOR,-one, FTOD_OC(1,1,NBI,RKIofKI(KI),RKQofKQ(KQ),ISP,1),1,1,&
>                                            desc_FTOD_OC,FTOD_OC(1,1,NBJ,RKIofKI(KJ),RKQofKQ(KQ),ISP2,2),1,1,&
344c471
<                                  CALL APPLY_DENOM(W,KI,KJ,KA,KB,NBI,NBJ,ISP,ISP2,LATT_CUR,NBANDSGW)
---
>                                  CALL APPLY_DENOM(W,KI,KJ,KA,KB,MNBI,MNBJ,ISP,ISP2,LATT_CUR,NBANDSGW)
375c502
<             D2=D2*(KPOINTS_ORIG%WTKPT(1)*KPOINTS_ORIG%WTKPT(1))*one
---
>             D2=D2*(WTKPT*WTKPT)*one
384,385c511,514
<                      IF (RNA<=VBMAX(ISP)) D2(NA,NB,ISP)=zero
<                      IF (RNB<=VBMAX(ISP)) D2(NA,NB,ISP)=zero
---
>                      IF (.NOT. (EMPTY_MP2_ORBITAL(W%FERTOT(RNA,RKA,ISP)))) D2(NA,NB,ISP)=zero
> !                     IF (RNA<=VBMAX(ISP)) D2(NA,NB,ISP)=zero
>                      IF (.NOT. (EMPTY_MP2_ORBITAL(W%FERTOT(RNB,RKA,ISP)))) D2(NA,NB,ISP)=zero
> !                     IF (RNB<=VBMAX(ISP)) D2(NA,NB,ISP)=zero
389a519,520
>                      IF (RNA>NBANDSHIGH_) D2(NA,NB,ISP)=zero
>                      IF (RNB>NBANDSHIGH_) D2(NA,NB,ISP)=zero
391,392c522,524
<                      IF ((RNA==RNB) .and. (RNA<=VBMAX(ISP)))   D2(NA,NB,ISP)=one*1000
<                      IF ((RNA==RNB) .and. (RNA<=VBMAX(ISP)))   D2(NA,NB,ISP)=one*1000
---
>                      IF ((RNA==RNB) .and. (.NOT. (EMPTY_MP2_ORBITAL(W%FERTOT(RNA,RKA,ISP))))) D2(NA,NB,ISP)=one*1000
> !                     IF ((RNA==RNB) .and. (RNA<=VBMAX(ISP)))   D2(NA,NB,ISP)=one*1000
> !                     IF ((RNA==RNB) .and. (RNA<=VBMAX(ISP)))   D2(NA,NB,ISP)=one*1000
516,517c648,651
<                      IF (RNA<=VBMAX(ISP)) D2(NA,NB,ISP)=zero
<                      IF (RNB<=VBMAX(ISP)) D2(NA,NB,ISP)=zero
---
>                      IF (.NOT. (EMPTY_MP2_ORBITAL(W%FERTOT(RNA,RKA,ISP)))) D2(NA,NB,ISP)=zero
> !                     IF (RNA<=VBMAX(ISP)) D2(NA,NB,ISP)=zero
>                      IF (.NOT. (EMPTY_MP2_ORBITAL(W%FERTOT(RNB,RKA,ISP)))) D2(NA,NB,ISP)=zero
> !                     IF (RNB<=VBMAX(ISP)) D2(NA,NB,ISP)=zero
519c653,654
<                      IF ((RNA==RNB) .and. (RNB<=VBMAX(ISP)))  D2(NA,NB,ISP)=one*1000
---
>                      IF ((RNA==RNB) .and. (.NOT. (EMPTY_MP2_ORBITAL(W%FERTOT(RNA,RKA,ISP))))) D2(NA,NB,ISP)=one*1000
> !                     IF ((RNA==RNB) .and. (RNB<=VBMAX(ISP)))  D2(NA,NB,ISP)=one*1000
541,542c676,679
<                      IF (RNA<=VBMAX(ISP)) D2(NA,NB,ISP)=zero
<                      IF (RNB<=VBMAX(ISP)) D2(NA,NB,ISP)=zero
---
>                      IF (.NOT. (EMPTY_MP2_ORBITAL(W%FERTOT(RNA,RKA,ISP)))) D2(NA,NB,ISP)=zero
>                      !IF (RNA<=VBMAX(ISP)) D2(NA,NB,ISP)=zero
>                      IF (.NOT. (EMPTY_MP2_ORBITAL(W%FERTOT(RNB,RKA,ISP)))) D2(NA,NB,ISP)=zero
>                      !IF (RNB<=VBMAX(ISP)) D2(NA,NB,ISP)=zero
547c684,685
<                      IF ((RNA==RNB) .and. (RNB<=VBMAX(ISP)))  D2(NA,NB,ISP)=one*1000
---
>                      IF ((RNA==RNB) .and. (.NOT. (EMPTY_MP2_ORBITAL(W%FERTOT(RNA,RKA,ISP))))) D2(NA,NB,ISP)=one*1000
>                      !IF ((RNA==RNB) .and. (RNB<=VBMAX(ISP)))  D2(NA,NB,ISP)=one*1000
609c747
<                ENDIF
---
>                   ENDIF
611c749
<             ENDDO
---
>                ENDDO
613c751
<             CALLMPI( M_sum_z(WGW%COMM_INTER, CELTOT_NEW(1,RKA,ISP), WDES%NB_TOT))
---
>                CALLMPI( M_sum_z(WGW%COMM_INTER, CELTOT_NEW(1,RKA,ISP), WDES%NB_TOT))
615,618c753,756
<             CALL FOCKM_OUT(IO,WGW,WDES,RKA,ISP)
<     
<          ENDIF
<          ENDDO
---
>                CALL FOCKM_OUT(IO,WGW,WDES,RKA,ISP)
>   
>                ENDIF
>             ENDDO
674c812
<          CALL EXIT(0)
---
> !         CALL EXIT(0)
677a816,842
> 
> !******************** SUBROUTINE NFREEZE_READER ************************
> !
> !***********************************************************************
>       SUBROUTINE NFREEZE_READER(IU5,IU6,IU0)
>       USE base
>       USE reader_tags
>       INTEGER IU5,IU6,IU0
>       INTEGER NTYP
>       ! local variables
>       INTEGER IDUM,N,IERR
>       REAL(q) RDUM
>       COMPLEX(q) CDUM
>       LOGICAL LOPEN,LDUM,LRPA
>       CHARACTER (1) CHARAC
>       CHARACTER (40) SZNAM
> 
>       LOPEN=.TRUE.
>       NFREEZE=0
>       CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NFREEZE', NFREEZE, IERR, WRITEXMLINCAR)
> 
>       RETURN
>       END SUBROUTINE NFREEZE_READER
> 
> 
> 
> 
720a886,891
>                IF (LOVERWRITEVBMAX) THEN
>                   VIRT=0.0_q
>                   OCC =0.0_q
>                   IF ( (NI<=VBMAX(ISP1)) .and. (NJ<=VBMAX(ISP2)) ) OCC=1._q
>                   IF ( (RNA>VBMAX(ISP1)) .and. (RNB>VBMAX(ISP2)) ) VIRT=1._q
>                ENDIF
727c898,899
<                IF (RNA<=VBMAX(ISP1)) THEN
---
>                IF (.NOT. (EMPTY_MP2_ORBITAL(W%FERTOT(RNA,KA,ISP1)))) THEN
> !               IF (RNA<=VBMAX(ISP1)) THEN
732c904,905
<                IF (RNB<=VBMAX(ISP2)) THEN
---
>                IF (.NOT. (EMPTY_MP2_ORBITAL(W%FERTOT(RNB,KB,ISP2)))) THEN
> !               IF (RNB<=VBMAX(ISP2)) THEN
736a910,915
>                IF (RNA>NBANDSHIGH_) TWOE4ORBITAL(NA,NB)=zero
>                IF (RNB>NBANDSHIGH_) TWOE4ORBITAL(NA,NB)=zero
>                IF (RNA>NBANDSHIGH_) TWOE4ORBITAL_X(NA,NB)=zero
>                IF (RNB>NBANDSHIGH_) TWOE4ORBITAL_X(NA,NB)=zero
> 
> 
770c949
<       SUBROUTINE CALC_2ORBITAL_FTOD(WDES,WGW,W,P,T_INFO,LATT_CUR,LMAXMP2,ENCUTGW,ENCUTGWSOFT,IO,FSG)
---
>       SUBROUTINE CALC_2ORBITAL_FTOD(WDES,WGW,W,P,T_INFO,INFO,LATT_CUR,LMAXMP2,ENCUTGW,ENCUTGWSOFT,IO,FSG)
784a964
>          TYPE (info_struct) INFO
796c976
<          INTEGER NBI,NBA,NBAA, i
---
>          INTEGER NBI,NBA,NBAA, i, MNBI
828c1008,1013
<          VBMAX(1)=LAST_FILLED_XI_NOMOD(W,1,1)
---
> !         VBMAX(1)=LAST_FILLED_XI_NOMOD(W,1,1)
>          VBMAX(1)=0
>          DO KI=1,W%WDES%NKPTS
>             VBMAX(1)=MAX(LAST_FILLED_XI_NOMOD(W,KI,1),VBMAX(1))
>          ENDDO
> 
830c1015,1019
<             VBMAX(2)=LAST_FILLED_XI_NOMOD(W,1,ISP)
---
>             VBMAX(2)=0
>             DO KI=1,W%WDES%NKPTS
>                VBMAX(2)=MAX(LAST_FILLED_XI_NOMOD(W,KI,1),VBMAX(2))
>             ENDDO
> !            VBMAX(2)=LAST_FILLED_XI_NOMOD(W,1,ISP)
831a1021,1029
>          LMETAL=.FALSE.
>          IF (VBMAX(1)>(INFO%NELECT/2)) LMETAL=.TRUE.
>          IF (IO%IU0>0) THEN
>             WRITE(*,*)'VBMAX=',VBMAX(1),'NELECT=',INT(INFO%NELECT),"LMETAL=",LMETAL
>          ENDIF
> 
> 
>          CALL VBMAX_INCAR_READER(IO%IU5, IO%IU6, IO%IU0)
> !         IF (IO%IU0>0) WRITE(IO%IU0,*)'VBMAX is set to : ',VBMAX(1)
869,870c1067,1068
<          mem_req=mem_req+NGVECTOR*PROCS*WDES%NBANDS*MAX(VBMAX(1),VBMAX(WDES%ISPIN))*WDES%NKPTS*WDES%NKPTS*WDES%ISPIN*2.0_q*16.0_q/1024.0_q/1024.0_q/1024.0_q 
<          IF (LMAXMP2>=0) mem_req=mem_req+NHVECTOR*PROCS*WDES%NBANDS*MAX(VBMAX(1),VBMAX(WDES%ISPIN))*WDES%NKPTS*WDES%NKPTS*WDES%ISPIN*2.0_q*16.0_q/1024.0_q/1024.0_q/1024.0_q
---
>          mem_req=mem_req+NGVECTOR*PROCS*WDES%NBANDS*MAX(VBMAX(1),VBMAX(WDES%ISPIN))*REALNKPTS*REALNKPTS*WDES%ISPIN*2.0_q*16.0_q/1024.0_q/1024.0_q/1024.0_q 
>          IF (LMAXMP2>=0) mem_req=mem_req+NHVECTOR*PROCS*WDES%NBANDS*MAX(VBMAX(1),VBMAX(WDES%ISPIN))*REALNKPTS*REALNKPTS*WDES%ISPIN*2.0_q*16.0_q/1024.0_q/1024.0_q/1024.0_q
896c1094,1099
<          ALLOCATE(tmp_FTOD_PW(NGVECTOR,WDES%NBANDS,MAX(VBMAX(1),VBMAX(WDES%ISPIN)),ncc))
---
>          IF (LFREEZE) THEN
>             ALLOCATE(tmp_FTOD_PW(NGVECTOR,WDES%NBANDS,MAX(VBMAX(1),VBMAX(WDES%ISPIN))-NFREEZE,ncc))
>          ELSE
>             ALLOCATE(tmp_FTOD_PW(NGVECTOR,WDES%NBANDS,MAX(VBMAX(1),VBMAX(WDES%ISPIN)),ncc))
>          ENDIF
> 
898c1101,1105
<             ALLOCATE(tmp_FTOD_OC(NHVECTOR,WDES%NBANDS,MAX(VBMAX(1),VBMAX(WDES%ISPIN)),2))
---
>             IF (LFREEZE) THEN
>                ALLOCATE(tmp_FTOD_OC(NHVECTOR,WDES%NBANDS,MAX(VBMAX(1),VBMAX(WDES%ISPIN))-NFREEZE,2))
>             ELSE
>                ALLOCATE(tmp_FTOD_OC(NHVECTOR,WDES%NBANDS,MAX(VBMAX(1),VBMAX(WDES%ISPIN)),2))
>             ENDIF
915a1123
>             IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0)) CYCLE
933a1142
>                IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI)==0)) CYCLE
944,950d1152
<                IF (OEAPPROXIMATE_NOs) THEN
<                   DO NBI=2,VBMAX(ISP)
<                      IF ((EMPTY_MP2_ORBITAL(W%FERTOT(NBI,KI,ISP)))) CYCLE 
<                      WI(1)%CR(:)=ABS(WI(1)%CR(:))+ABS(WI(NBI)%CR(:))
<                      WI(1)%CPROJ(:)=WI(1)%CPROJ(:)+WI(NBI)%CPROJ(:)
<                   ENDDO
<                ENDIF
986,987c1188,1191
<                   DO NBI=1,MAX(VBMAX(1),VBMAX(WDES%ISPIN))
<                      IF ((OEAPPROXIMATE_NOs) .and. (NBI>1)) CYCLE
---
>                   DO MNBI=1,MAX(VBMAX(1),VBMAX(WDES%ISPIN))
>                      NBI=MNBI
>                      IF ((LFREEZE) .and. (NFREEZE>=MNBI))  CYCLE
>                      IF (LFREEZE) NBI=MNBI-NFREEZE
1001c1205
<                            CALL FOCK_CHARGE_ONE_CENTER_NOINT( WI(NBI),WA(NBAA),&
---
>                            CALL FOCK_CHARGE_ONE_CENTER_NOINT( WI(MNBI),WA(NBAA),&
1004c1208
<                            CALL FOCK_CHARGE_NOINT( WI(NBI),WA(NBAA), GWORK(1), &
---
>                            CALL FOCK_CHARGE_NOINT( WI(MNBI),WA(NBAA), GWORK(1), &
1057c1261
<                         CALL FOCK_CHARGE_ONE_CENTER_NOINT( WA(NBAA),WI(NBI),&
---
>                         CALL FOCK_CHARGE_ONE_CENTER_NOINT( WA(NBAA),WI(MNBI),&
1068c1272
<                            tmp_FTOD_OC(1:H%TOTAL_ENTRIES,NBAA+NBA-1,NBI,2)=(CRHOIA(1:H%TOTAL_ENTRIES,NBAA))    !000000000000000000000
---
>                            tmp_FTOD_OC(1:H%TOTAL_ENTRIES,NBAA+NBA-1,MNBI,2)=(CRHOIA(1:H%TOTAL_ENTRIES,NBAA))    !000000000000000000000
1099,1100c1303,1307
<                   DO NBI=1,MAX(VBMAX(1),VBMAX(WDES%ISPIN))
<                      IF ((OEAPPROXIMATE_NOs) .and. (NBI>1)) CYCLE
---
>                   DO MNBI=1,MAX(VBMAX(1),VBMAX(WDES%ISPIN))
>                      NBI=MNBI
>                      IF ((LFREEZE) .and. (NFREEZE>=MNBI))  CYCLE
>                      IF (LFREEZE) NBI=MNBI-NFREEZE
>  
1113c1320
<                            CALL FOCK_CHARGE_ONE_CENTER_NOINT( WB(NBAA),WI(NBI),&
---
>                            CALL FOCK_CHARGE_ONE_CENTER_NOINT( WB(NBAA),WI(MNBI),&
1116c1323
<                            CALL FOCK_CHARGE_NOINT( WB(NBAA),WI(NBI), GWORK(1), &
---
>                            CALL FOCK_CHARGE_NOINT( WB(NBAA),WI(MNBI), GWORK(1), &
1192a1400
>          NPROW_GRID=1
1206,1207c1414,1415
<         a_PRCS=CEILING(SQRT(Real(PROCS)))
<         IF (a_PRCS==SQRT(Real(PROCS))) THEN
---
>         a_PRCS=CEILING(SQRT(Real(PROCS,kind=q)))
>         IF (a_PRCS==SQRT(Real(PROCS,kind=q))) THEN
1211,1213c1419,1421
<         IF (a_PRCS/=SQRT(Real(PROCS))) THEN
<            DO i=1,CEILING(SQRT(Real(PROCS)))
<               NPCOL_TMP=PROCS/Real(i)
---
>         IF (a_PRCS/=SQRT(Real(PROCS,kind=q))) THEN
>            DO i=1,CEILING(SQRT(Real(PROCS,kind=q)))
>               NPCOL_TMP=PROCS/Real(i,kind=q)
1221,1222c1429,1430
<            WRITE(*,*)'The allocated number of CPUs does not allow for the use of an efficient process grid.'
<            WRITE(*,*)'Suggested number of CPUs: 4, 16, 32, ....'
---
>            IF ((MYROW==0) .AND. (MYCOL==0)) WRITE(*,*)'The allocated number of CPUs does not allow for the use of an efficient process grid.'
>            IF ((MYROW==0) .AND. (MYCOL==0)) WRITE(*,*)'Suggested number of CPUs: 4, 16, 32, ....'
1288c1496
<          INTEGER :: NBI,FTOD_PW_rows_br,VBMAX_tmp
---
>          INTEGER :: NBI,FTOD_PW_rows_br,VBMAX_tmp,MNBI
1309,1310c1517,1520
<          IF (OEAPPROXIMATE_NOs) VBMAX_tmp=1
<          DO NBI=1,VBMAX_tmp
---
>          DO MNBI=1,VBMAX_tmp
>             NBI=MNBI
>             IF ((LFREEZE) .and. (NFREEZE>=MNBI))  CYCLE
>             IF (LFREEZE) NBI=MNBI-NFREEZE
1313c1523
<                 desc_FTOD_PW_br,FTOD_PW(1,1,NBI,KI,KQ,ISP,cc),1,1,desc_FTOD_PW,contxt_grid)
---
>                 desc_FTOD_PW_br,FTOD_PW(1,1,NBI,RKIofKI(KI),RKQofKQ(KQ),ISP,cc),1,1,desc_FTOD_PW,contxt_grid)
1340,1341c1550,1554
<             IF (OEAPPROXIMATE_NOs) VBMAX_tmp=1
<             DO NBI=1,VBMAX_tmp
---
>             DO MNBI=1,VBMAX_tmp
>                NBI=MNBI
>                IF ((LFREEZE) .and. (NFREEZE>=MNBI))  CYCLE
>                IF (LFREEZE) NBI=MNBI-NFREEZE
>  
1348c1561
<                    desc_FTOD_OC_br,FTOD_OC(1,1,NBI,KI,KQ,ISP,cc),1,1,desc_FTOD_OC,contxt_grid)
---
>                    desc_FTOD_OC_br,FTOD_OC(1,1,NBI,RKIofKI(KI),RKQofKQ(KQ),ISP,cc),1,1,desc_FTOD_OC,contxt_grid)
1550,1551c1763,1766
<          IF (OEAPPROXIMATE_NOs) THEN
<             ALLOCATE(FTOD_PW(FTOD_PW_rows,FTOD_PW_cols,1,WDES%NKPTS,WDES%NKPTS,WDES%ISPIN,ncc))
---
>          IF ((LFREEZE) .and. (.not. LLOWMEM)) THEN
>             ALLOCATE(FTOD_PW(FTOD_PW_rows,FTOD_PW_cols,MAX(VBMAX(1),VBMAX(WDES%ISPIN))-NFREEZE,REALNKPTS,REALNKPTS,WDES%ISPIN,ncc))
>          ELSE IF (LLOWMEM) THEN
>             ALLOCATE(FTOD_PW(FTOD_PW_rows,FTOD_PW_cols,1,1,1,1,ncc))
1553c1768
<             ALLOCATE(FTOD_PW(FTOD_PW_rows,FTOD_PW_cols,MAX(VBMAX(1),VBMAX(WDES%ISPIN)),WDES%NKPTS,WDES%NKPTS,WDES%ISPIN,ncc))
---
>             ALLOCATE(FTOD_PW(FTOD_PW_rows,FTOD_PW_cols,MAX(VBMAX(1),VBMAX(WDES%ISPIN)),REALNKPTS,REALNKPTS,WDES%ISPIN,ncc))
1576,1577c1791,1794
<             IF (OEAPPROXIMATE_NOs) THEN            
<                ALLOCATE(FTOD_OC(FTOD_OC_rows,FTOD_OC_cols,1,WDES%NKPTS,WDES%NKPTS,WDES%ISPIN,2))
---
>             IF ((LFREEZE) .and. (.not. LLOWMEM)) THEN
>                ALLOCATE(FTOD_OC(FTOD_OC_rows,FTOD_OC_cols,MAX(VBMAX(1),VBMAX(WDES%ISPIN))-NFREEZE,REALNKPTS,REALNKPTS,WDES%ISPIN,2))
>             ELSE IF (LLOWMEM) THEN
>                ALLOCATE(FTOD_OC(FTOD_OC_rows,FTOD_OC_cols,1,1,1,1,2))
1579c1796
<                ALLOCATE(FTOD_OC(FTOD_OC_rows,FTOD_OC_cols,MAX(VBMAX(1),VBMAX(WDES%ISPIN)),WDES%NKPTS,WDES%NKPTS,WDES%ISPIN,2))
---
>                ALLOCATE(FTOD_OC(FTOD_OC_rows,FTOD_OC_cols,MAX(VBMAX(1),VBMAX(WDES%ISPIN)),REALNKPTS,REALNKPTS,WDES%ISPIN,2))
1693a1911,2478
>     SUBROUTINE CHECK_SHIFTED_KPOINTS(WDES,W,KPOINTS)
>          USE constant
>          USE full_kpoints
>          USE mkpoints
>          USE wave
>          implicit NONE
>          TYPE(wavedes) WDES
>          TYPE(wavespin) W
>          TYPE (kpoints_struct) KPOINTS
>          integer :: MKI,KI
>          LOGICAL :: GAMMA_FOUND
>          LOGICAL :: SECOND_GAMMA
>          REAL(q) :: SECOND_GAMMA_WEIGHT
>          REAL(q) :: GAMMA_WEIGHT
>          INTEGER :: NRKI, NRKQ
>          REAL(q) :: TINY=1.D-8 
>          
>          GAMMA_FOUND=.FALSE.
>          !write(*,*)LSHIFT_KPOINTS
>          SHIFTED_KPOINTS=.TRUE.
>          SECOND_GAMMA=.FALSE.
>          do MKI=1,WDES%NKPTS
>             if ( (ABS(W%WDES%VKPT(1,MKI))<(TINY)) .and. (ABS(W%WDES%VKPT(2,MKI))<(TINY)) .and. (ABS(W%WDES%VKPT(3,MKI))<(TINY))) then
>                IF (GAMMA_FOUND) SECOND_GAMMA=.TRUE.
>                IF (SECOND_GAMMA) SECOND_GAMMA_WEIGHT=W%WDES%WTKPT(MKI)
>                IF (.NOT. GAMMA_FOUND) THEN
>                   GAMMA_WEIGHT=W%WDES%WTKPT(MKI)
>                   GAMMA_FOUND=.true.
>                ENDIF
>             endif
>          enddo
>          !write(*,*)gamma_found,second_gamma,gamma_weight
>          if ((GAMMA_FOUND) .and. (GAMMA_WEIGHT/=0.0_q)) THEN
>             SHIFTED_KPOINTS=.FALSE.
>          endif
>          IF (SECOND_GAMMA) SHIFTED_KPOINTS=.TRUE.
>          !CALL SETUP_KPOINTS_STATIC(KPOINTS)
>          !CALL CHECK_FULL_KPOINTS ! all set up properly ?
> 
>          REALNKPTS=WDES%NKPTS
>          IF (SHIFTED_KPOINTS) REALNKPTS=WDES%NKPTS/2
> 
>          ALLOCATE(RKQofKQ(WDES%NKPTS))
>          ALLOCATE(RKIofKI(WDES%NKPTS))
>          ALLOCATE(KQofRKQ(REALNKPTS))
>          ALLOCATE(KIofRKI(REALNKPTS))
>          NRKI=0
>          NRKQ=0
>          RKQofKQ=0
>          RKIofKI=0
>          KQofRKQ=0
>          KIofRKI=0
>          IF (SHIFTED_KPOINTS) THEN
>             DO KI=1,WDES%NKPTS
>                IF ((WDES%WTKPT(KI)/=0.00_q))  THEN
>                   NRKI=NRKI+1
>                   RKIofKI(KI)=NRKI
>                   KIofRKI(NRKI)=KI
>                ENDIF
>                IF ((WDES%WTKPT(KI)==0.00_q))  THEN
>                   NRKQ=NRKQ+1
>                   RKQofKQ(KI)=NRKQ
>                   KQofRKQ(NRKQ)=KI
>                ENDIF
> 
>             !write(*,*)LSHIFT_KPOINTS,W%WDES%WTKPT(MKI)
>             !write(*,*) mki,',',W%WDES%VKPT(:,MKI),W%WDES%WTKPT(MKI)
>             ENDDO
>          ELSE
>             DO KI=1,WDES%NKPTS
>                RKIofKI(KI)=KI
>                KIofRKI(KI)=KI
>                RKQofKQ(KI)=KI
>                KQofRKQ(KI)=KI
>             ENDDO
>          ENDIF
> 
>         WTKPT=1.0_q/REALNKPTS 
>          
>     END SUBROUTINE CHECK_SHIFTED_KPOINTS
> 
> !******************** SUBROUTINE NOINCAR_READER ************************
> !
> !***********************************************************************
>     SUBROUTINE NOINCAR_READER(IU5,IU6,IU0)
>       USE base
>       USE reader_tags
>       INTEGER IU5,IU6,IU0
>       INTEGER NTYP
>       ! local variables
>       INTEGER IDUM,N,IERR
>       REAL(q) RDUM
>       COMPLEX(q) CDUM
>       LOGICAL LOPEN,LDUM,LRPA
>       CHARACTER (1) CHARAC
>       CHARACTER (40) SZNAM
>       CHARACTER (40)   STRING
>       INTEGER, EXTERNAL :: LENGTH
> 
>       LOPEN=.TRUE.
>       APPROXIMATE_NOs=.FALSE.
>       CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LAPPROX', APPROXIMATE_NOs, IERR, WRITEXMLINCAR)
>       CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBANDSHIGH', NBANDSHIGH_, IERR, WRITEXMLINCAR)
>       RETURN
> 
>     END SUBROUTINE NOINCAR_READER
> 
> !******************** SUBROUTINE NOINCAR_READER ************************
> !
> !***********************************************************************
>     SUBROUTINE VBMAX_INCAR_READER(IU5,IU6,IU0)
>       USE base
>       USE reader_tags
>       INTEGER IU5,IU6,IU0
>       INTEGER NTYP
>       ! local variables
>       INTEGER IDUM,N,IERR
>       REAL(q) RDUM
>       COMPLEX(q) CDUM
>       LOGICAL LOPEN,LDUM,LRPA
>       CHARACTER (1) CHARAC
>       CHARACTER (40) SZNAM
>       CHARACTER (40)   STRING
>       INTEGER, EXTERNAL :: LENGTH
>       INTEGER :: MVBMAX
> 
> !      LOPEN=.TRUE.
> !      MVBMAX=0
> !      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'VBMAX', MVBMAX, IERR, WRITEXMLINCAR)
> !      IF (MVBMAX/=0) THEN
> !         LOVERWRITEVBMAX=.TRUE.
> !         VBMAX(:)=MVBMAX
> !      ENDIF
>       RETURN
> 
>     END SUBROUTINE VBMAX_INCAR_READER
> 
> 
>       !***********************************************************************
>       !This routine calculates the fourier-transformed overlap integrals <i|-G|a>
>       !and <j|G|b>*4*pi*e^2/(G+q)^2 and also calls the redistribution routine.
>       !This is done for the plane-wave part, as well as for the one-center terms.
>       !Note, that in the gamma-only version it is enough to store <i|-G|a>*SQRT(4*pi*e^2/(G+q)^2)
>       !for the plane-wave part, because (<i|-G|a>)*=<i|G|a>.
>       !However, this approximation cannot be made for the one-center terms in the gamma-only
>       !version, only for technical reasons.
>       !The <i|G|a> quantities are most important for the construction of the two electron
>       !four orbital integrals <ij|ab>.  <ij|ab>=4pi e^2 \sum_G <i|-G|a> <j|G|b> /(G+q)^2
>       !***********************************************************************
> 
>       SUBROUTINE CALC_2ORBITAL_FTOD_LOWMEM(WDES,WGW,W,P,T_INFO,INFO,LATT_CUR,LMAXMP2,ENCUTGW,ENCUTGWSOFT,IO,FSG,INBI,IKI,IKQ,IISP)
>          USE prec
>          USE poscar
>          USE pseudo
>          USE wave_high
>          USE full_kpoints
>          USE mkpoints
>          USE lattice
>          USE constant
>          USE base
>          IMPLICIT NONE
>          TYPE(wavedes) WDES
>          TYPE(wavedes) WGW
>          TYPE(wavespin) W
>          TYPE(type_info) T_INFO
>          TYPE (info_struct) INFO
>          TYPE(potcar) P(T_INFO%NTYP)
>          TYPE(latt) LATT_CUR
>          INTEGER LMAXMP2
>          TYPE (in_struct) IO
>          INTEGER :: INBI, IKI, IKQ, IISP
>          ! local variables
>          TYPE(wavespin), SAVE :: WHF
>          TYPE(wavedes1), SAVE :: WGWQ
>          TYPE(wavedes1), SAVE ::  WDESKI,WDESKA,WDESKB
>          TYPE(wavefun1), ALLOCATABLE, SAVE :: WI(:), WA(:), WB(:)
>          INTEGER KA,KB,KI_IN_FULL_ORIG,KQ_, NBI, KI, KQ
>          INTEGER, SAVE :: NSTRIP, NSTRIPA, ISP, NFFT
>          INTEGER NBA,NBAA, i, MNBI
>          INTEGER NP ! number of plane waves for the overlap density i*(r)a(r) GCHGIA(r)
>          REAL(q) :: FSG ! singularity correction
>          REAL(q) :: ENCUTGW,ENCUTGWSOFT 
>          REAL(q) :: POTFAK(GRIDHF%MPLWV) 
>          COMPLEX(q) :: CPHASE(GRIDHF%MPLWV)
>          COMPLEX(q) :: CPHASE2(GRIDHF%MPLWV)
>          LOGICAL :: LPHASE  
>          LOGICAL :: LPHASE2
>          COMPLEX(q) :: GWORK(MAX(GRIDHF%MPLWV,WGW%GRID%MPLWV)) !work array for calculating GCHGIA
>          COMPLEX(q), ALLOCATABLE, SAVE :: GCHGIA(:,:,:)  ! charge
>          GDEF      , ALLOCATABLE, SAVE :: CRHOIA(:,:)    ! one-center charge
>          GDEF      , ALLOCATABLE, SAVE :: CRHOIB(:,:)    ! one-center charge
>          GDEF      , ALLOCATABLE, SAVE :: CRHOLM(:)      ! augmentation occupancy matrix
>          COMPLEX(q), ALLOCATABLE, SAVE :: tmp_FTOD_PW(:,:,:,:)
>          GDEF      , ALLOCATABLE, SAVE :: tmp_FTOD_OC(:,:,:,:)
>          Real(q) :: mem_req
>  
>          NBI=INBI
>          KI=IKI
>          KQ=IKQ
>          ISP=IISP
> 
>          IF (INBI==0) THEN
> 
>             CALL CHECK_FULL_KPOINTS
>          
>             IF (LMAXMP2>=0) THEN
>               CALL SET_UP_ONE_CENTER_H(WDES,P,T_INFO,LMAXMP2,H)
>             ENDIF
>  
>             WHF=W
>             WHF%WDES => WDES_FOCK
>             NSTRIP=30
>             CALL SETWDES(WHF%WDES,WDESKI,0)
>             CALL SETWDES(WHF%WDES,WDESKA,0)
>             CALL SETWDES(WHF%WDES,WDESKB,0)
>          
> !            VBMAX(1)=LAST_FILLED_XI_NOMOD(W,1,1)
>             VBMAX(1)=0
>             DO KI=1,W%WDES%NKPTS
>                VBMAX(1)=MAX(LAST_FILLED_XI_NOMOD(W,KI,1),VBMAX(1))
>             ENDDO
>             DO ISP=2,WDES%ISPIN
>                VBMAX(2)=0
>                DO KI=1,W%WDES%NKPTS
>                   VBMAX(2)=MAX(LAST_FILLED_XI_NOMOD(W,KI,1),VBMAX(2))
>                ENDDO
> !               VBMAX(2)=LAST_FILLED_XI_NOMOD(W,1,ISP)
>             ENDDO
>             LMETAL=.FALSE.
>             IF (VBMAX(1)>(INFO%NELECT/2)) LMETAL=.TRUE.
>             IF (IO%IU0>0) THEN
>                WRITE(*,*)'VBMAX=',VBMAX(1),'NELECT=',INFO%NELECT,"LMETAL=",LMETAL
>             ENDIF
> 
> 
>             CALL VBMAX_INCAR_READER(IO%IU5, IO%IU6, IO%IU0)
> !            WRITE(*,*)'VBMAX is set to : ',VBMAX(1)
> 
>             ALLOCATE(WI(MAX(VBMAX(1),VBMAX(WDES%ISPIN))),WA(NSTRIP),WB(NSTRIP))
>             DO NBI=1,MAX(VBMAX(1),VBMAX(WDES%ISPIN))
>                CALL NEWWAV(WI(NBI),WDESKI,.TRUE.)
>             ENDDO
>             DO NBA=1,NSTRIP
>                CALL NEWWAV(WA(NBA),WDESKA,.TRUE.)
>             ENDDO
>             DO NBA=1,NSTRIP
>                CALL NEWWAV(WB(NBA),WDESKB,.TRUE.)
>             ENDDO
>           
>             NGVECTOR=MAXVAL(WGW%NGVECTOR(:))
>             NHVECTOR=0
>             IF (ASSOCIATED(H)) THEN
>                NHVECTOR=H%TOTAL_ENTRIES
>             ENDIF
>             ! Make sure that NGVECTOR is suited for the BLACS process grid
>             IF (MOD(NGVECTOR,NPROW_GRID)/=0) THEN
>                NGVECTOR=NGVECTOR+(NPROW_GRID-MOD(NGVECTOR,NPROW_GRID))
>             ENDIF
>             ! Make sure that NHVECTOR is suited for the BLACS process grid
>             IF (MOD(NHVECTOR,NPROW_GRID)/=0) THEN
>                NHVECTOR=NHVECTOR+(NPROW_GRID-MOD(NHVECTOR,NPROW_GRID))
>             ENDIF
>    
>             IF (ALLOCATED(FTOD_PW)) THEN
>                DEALLOCATE(FTOD_PW)
>             ENDIF
>          
>             !Initialize the 2D process grid
>             CALL INIT_BLACS_GRID(WDES)
>    
>             CALL SETUP_FTOD(WDES)
>  
> 
>             !write out an estimate for the required memory
>             mem_req=2.0_q*(PROCS*PROCS*WDES%NBANDS*WDES%NBANDS)*16.0_q/1024.0_q/1024.0_q/1024.0_q !TWOE4ORBITAL arrays
>             mem_req=mem_req+NGVECTOR*PROCS*WDES%NBANDS*MAX(VBMAX(1),VBMAX(WDES%ISPIN))*REALNKPTS*REALNKPTS*WDES%ISPIN*2.0_q*16.0_q/1024.0_q/1024.0_q/1024.0_q 
>             IF (LMAXMP2>=0) mem_req=mem_req+NHVECTOR*PROCS*WDES%NBANDS*MAX(VBMAX(1),VBMAX(WDES%ISPIN))*REALNKPTS*REALNKPTS*WDES%ISPIN*2.0_q*16.0_q/1024.0_q/1024.0_q/1024.0_q
> !            IF (IO%IU0>0) THEN
> !               WRITE(*,"(A,E10.3,A)") 'For MP2 calculations approximately',mem_req,'GB RAM will be required.' 
> !            ENDIF
>             IF (IO%IU0>0) THEN
> !               WRITE(IO%IU0,*) 'Allocating memory...'
>             ENDIF
>             !Setup the descriptors for the distributed TWOE4ORBITAL matrix   
>             !and allocate the TWOE4ORBITAL and TWOE4ORBITAL_X matrices
>             CALL SETUP_TWOE4ORBITAL(WDES)
>             IF (ALLOCATED(TWOE4ORBITAL)) THEN
>                DEALLOCATE(TWOE4ORBITAL)
>             ENDIF
>             IF (ALLOCATED(D2)) THEN
>                DEALLOCATE(D2)
>             ENDIF
>             IF (ALLOCATED(D2EV)) THEN
>                DEALLOCATE(D2EV)
>             ENDIF
>             IF (ALLOCATED(TWOE4ORBITAL_X)) THEN
>                DEALLOCATE(TWOE4ORBITAL_X)
>             ENDIF
>             IF (IO%IU0>0) THEN
> !               WRITE(IO%IU0,*) 'succeeded'
>             ENDIF
>          
>             IF (LFREEZE) THEN
>                ALLOCATE(tmp_FTOD_PW(NGVECTOR,WDES%NBANDS,1,ncc))
>             ELSE
>                ALLOCATE(tmp_FTOD_PW(NGVECTOR,WDES%NBANDS,1,ncc))
>             ENDIF
> 
>             IF (ASSOCIATED(H)) THEN
>                IF (LFREEZE) THEN
>                   ALLOCATE(tmp_FTOD_OC(NHVECTOR,WDES%NBANDS,1,2))
>                ELSE
>                   ALLOCATE(tmp_FTOD_OC(NHVECTOR,WDES%NBANDS,1,2))
>                ENDIF
>             ENDIF 
> 
>          
>             IF (IO%IU0>0) THEN
>                WRITE(IO%IU0,*)
>                WRITE(IO%IU0,*)'Calculating fourier transformed overlap densities:'
>             ENDIF
>          
>          
>             call BLACS_GRIDINFO(CONTXT, NPROW, NPCOL, MYROW, MYCOL)
>       
>             IF (ASSOCIATED(H)) THEN
>                ALLOCATE(CRHOIA(NHVECTOR, NSTRIP),CRHOIB(NHVECTOR, NSTRIP))
>                FTOD_OC=zero           
>             ENDIF
> 
>             RETURN
>          ENDIF
> 
>          NBI=INBI
>          KI=IKI
>          KQ=IKQ
>          ISP=IISP
> 
> 
> !         spin: DO ISP=1,WDES%ISPIN
> !         kqloop: DO KQ=1,WDES%NKPTS
>             IF (IO%IU0>=0) THEN
>                WRITE(IO%IU6,'("NI,KI=",2I4)') NBI,KI
>                IF (WDES%ISPIN==1) THEN
> !                 WRITE(IO%IU0,'("NQ=",I4,3F10.4,", ")') KQ,WDES%VKPT(:,KQ)
>                ELSE
> !                 WRITE(IO%IU0,'("NQ=",I4,3F10.4,", ",A1,A1,", ")') KQ,WDES%VKPT(:,KQ),ISP
>                ENDIF
>             ENDIF
>             CALL SETWDES(WGW,WGWQ,KQ)
>          
>             NP=WGWQ%NGVECTOR            
>             IF (NP>NGVECTOR) THEN
>                WRITE(*,*)'Internal error in "Calc_2orbital_FTOD": NP larger than NGVECTOR'
>                STOP
>             ENDIF
>             ALLOCATE(GCHGIA(NP,NSTRIP,2),CRHOLM(AUG_DES%NPRO*WDES%NRSPINORS))
>          
>  !           kiloop: DO KI=1,WDES%NKPTS
>                tmp_FTOD_PW=(0._q,0._q)
>                IF (ASSOCIATED(H)) tmp_FTOD_OC=zero
>                
> !               CALL GWPROGRESS(IO%IU0, KI,WDES%NKPTS,KQ,WDES%NKPTS)
>                KI_IN_FULL_ORIG=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,KI),KPOINTS_FULL_ORIG)
>                ! collect all valence bands at k_i
>                CALL SETWDES(WHF%WDES,WDESKI,KI)
>                CALL W1_GATHER_GLB(WHF,1,VBMAX(ISP),ISP,WI)
>                ! k_b = k_i - k_q - G
>                KB=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KQ)+WDES%VKPT(:,KI),KPOINTS_FULL)
> 
>                ! k_a = k_i + k_q - G
>                KA=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KI)-WDES%VKPT(:,KQ),KPOINTS_FULL)
>                
>                CALL SETWDES(WHF%WDES,WDESKA,KA)
>                
>                ! CPHASE(r) = e^iGr, where G = k_i - k_q - k_b
>                CALL SETPHASE(WDES%VKPT(:,KI)-WDES%VKPT(:,KQ)-WDES%VKPT(:,KA),GRIDHF,CPHASE,LPHASE)
>                
>                CALL SET_GFAC_WITHOUT_WEIGHT(GRIDHF,LATT_CUR,KI,KA,FSG,POTFAK)
>                ! 1/(G+q)**2
>                
>                IF (ENCUTGW /= ENCUTGWSOFT .AND. ENCUTGW > 0 .AND. ENCUTGWSOFT > 0) THEN
>                   CALL SET_GFAC_WAVEFUN(WGWQ,LATT_CUR,FSG,POTFAK,ENCUTGW,ENCUTGWSOFT)
>                ELSE
>                   CALL SET_GFAC_WAVEFUN(WGWQ,LATT_CUR,FSG,POTFAK)
>                ENDIF
>                
> #ifdef gammareal
>                   POTFAK=SQRT(POTFAK)
> #endif
>                
>                ! loop over all bands
>                DO NBA=1,WDES%NBANDS,NSTRIP
>                   NSTRIPA=MIN(WDES%NBANDS+1-NBA,NSTRIP)
>                   ! FFT{psi_a} to real space
>                   DO NBAA=1,NSTRIPA !copy and fourier transform NSTRIP wave functions
>                      CALL W1_COPY( ELEMENT(WHF,WDESKA,NBA+NBAA-1,ISP),WA(NBAA))
>                      CALL FFTWAV_W1(WA(NBAA))
>                   ENDDO
>                   ! loop over valence bands only
>  !                 DO MNBI=1,MAX(VBMAX(1),VBMAX(WDES%ISPIN))
>                      MNBI=NBI
> !                     IF (LFREEZE) NBI=MNBI-NFREEZE
>                      GCHGIA=0 
>                      
>                      IF (ASSOCIATED(H)) THEN
>                         CRHOIA=zero
>                         CRHOIB=zero
>                         CRHOLM=zero
>                      ENDIF
>                      !loop over all bands in NSTRIP 
>                      DO NBAA=1,NSTRIPA
>                      
>                      ! calculate rho(r)=psi_i(r)* psi_a(r) for,
>                      ! one center terms and, on the plane wave grid.
>                         IF (ASSOCIATED(H)) THEN
>                            CALL FOCK_CHARGE_ONE_CENTER_NOINT( WI(MNBI),WA(NBAA),&
>                              GWORK(1),H,CRHOIA(1,NBAA), CRHOLM,  SIZE(CRHOLM))
>                         ELSE
>                            CALL FOCK_CHARGE_NOINT( WI(MNBI),WA(NBAA), GWORK(1), &
>                              CRHOLM, SIZE(CRHOLM))
>                         ENDIF
>                         ! Set phase e^iGr, where G = k_i - k_q - k_b
>                         IF (LPHASE) THEN
>                            CALL APPLY_PHASE( GRIDHF, CPHASE(1), GWORK(1), &
>                              GWORK(1))
>                         ENDIF
>                    
>                      ! FFT{rho} to reciprocal space
>                         CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), &
>                           GWORK(1),GCHGIA(1,NBAA,1),WGWQ%GRID,.FALSE.)
>                         NFFT=NFFT+1
>                      ! multiply with potential factor
>                      
>                         CALL APPLY_GFAC_WAVEFUN(WGWQ,GCHGIA(1,NBAA,1), &
>                          POTFAK(1))
>                          
>                      ENDDO !NBAA (loop over all bands in NSTRIP)
>                      
>                      IF (ASSOCIATED(H)) THEN
>                         CALL APPLY_PHASE_ONE_CENTER(WHF%WDES, H, CRHOIA(:,:), &
>                           WHF%WDES%VKPT(:,KI)-WHF%WDES%VKPT(:,KA))
>                      ENDIF
>                      
>                      IF (ASSOCIATED(H)) THEN
>                         ! use CRHOIA as temporary work array
>                         CALL APPLY_ONE_CENTER_H( WHF%WDES, H, CRHOIA(:,:), CRHOIB(:,:), NSTRIPA)
>                      ENDIF
>                      
>                      DO NBAA=1,NSTRIPA
>                         !copy FTOD functions to FTOD_PW and FTOD_OC
> #ifdef gammareal
>                            tmp_FTOD_PW(1:NP,NBAA+NBA-1,1,1)=(GCHGIA(1:NP,NBAA,1))*SQRT(1.0_q/GRIDHF%NPLWV)
> #else
>                            tmp_FTOD_PW(1:NP,NBAA+NBA-1,1,1)=(GCHGIA(1:NP,NBAA,1))
> #endif
>                         IF (ASSOCIATED(H)) THEN
>                            tmp_FTOD_OC(1:H%TOTAL_ENTRIES,NBAA+NBA-1,1,1)=(CRHOIA(1:H%TOTAL_ENTRIES,NBAA))
>                         ENDIF
> 
>                      ENDDO
>                      
> !!!!!!!!!!!!!!!!! GAMMA-only version !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
> #ifdef gammareal
>                                      
>                      ! calculate rho(r)=psi_i(r)* psi_a(r) for,
>                      ! (1._q,0._q) center terms and, on the plane wave grid.
>                 IF (ASSOCIATED(H)) THEN
>                      CRHOIB=0
>                      CRHOLM=0
>                      !loop over all bands in NSTRIP 
>                      DO NBAA=1,NSTRIPA
>                         CALL FOCK_CHARGE_ONE_CENTER_NOINT( WA(NBAA),WI(MNBI),&
>                           GWORK(1),H,CRHOIA(1,NBAA), CRHOLM,  SIZE(CRHOLM))                
>                      ! FFT{rho} to reciprocal space
>                         CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), &
>                           GWORK(1),GCHGIA(1,NBAA,1),WGWQ%GRID,.FALSE.)
>                         NFFT=NFFT+1
>                      ! multiply with potential factor
>                      ENDDO !NBAA (loop over all bands in NSTRIP)
>                      DO NBAA=1,NSTRIPA
>                         !copy FTOD functions to FTOD_PW and FTOD_OC
>                         IF (ASSOCIATED(H)) THEN
>                            tmp_FTOD_OC(1:H%TOTAL_ENTRIES,NBAA+NBA-1,1,2)=(CRHOIA(1:H%TOTAL_ENTRIES,NBAA))    !000000000000000000000
>                         ENDIF
>                      ENDDO
>                  ENDIF
> #endif
> !!!!!!!!!!!!!!!!!!!!!!!!!!GAMMA only version !!!!!!!!!!!!!!!!!!!!!!!!
> 
> !                  ENDDO !NBI (loop over valence bands only)
>                ENDDO !NBA (loop over all bands)
>                
> #ifdef gammareal
>                
> #else
>                
>                CALL SETWDES(WHF%WDES,WDESKB,KB)        
>                
>                ! k_i - k_a = k_q + G
>                CALL PHASER_HF(GRIDHF,LATT_CUR,FAST_AUG_FOCK,WDES%VKPT(:,KB)-WDES%VKPT(:,KI))
>                
>                ! CPHASE(r) = e^iGr, where G = k_i - k_q - k_a
>                CALL SETPHASE(WDES%VKPT(:,KB)-WDES%VKPT(:,KQ)-WDES%VKPT(:,KI),GRIDHF,CPHASE,LPHASE)
>                
>                ! loop over all bands
>                DO NBA=1,WDES%NBANDS,NSTRIP
>                   NSTRIPA=MIN(WDES%NBANDS+1-NBA,NSTRIP)
>                   ! FFT{psi_a} to real space
>                   DO NBAA=1,NSTRIPA !copy and fourier transform NSTRIP wave functions
>                      CALL W1_COPY( ELEMENT(WHF,WDESKB,NBA+NBAA-1,ISP),WB(NBAA))
>                      CALL FFTWAV_W1(WB(NBAA))
>                   ENDDO
>                   ! loop over valence bands only
> !                  DO MNBI=1,MAX(VBMAX(1),VBMAX(WDES%ISPIN))
>                      MNBI=NBI
> !                     IF ((LFREEZE) .and. (NFREEZE>=MNBI))  CYCLE
>  
>                      GCHGIA=0 
>                      
>                      IF (ASSOCIATED(H)) THEN
>                         CRHOIB=0
>                         CRHOLM=0
>                      ENDIF
>                      !loop over all bands in NSTRIP
>                      DO NBAA=1,NSTRIPA
>                      ! GCHGIA number 2 for X-changed waves
>                      ! calculate rho(r)=psi_i(r)* psi_a(r) for,
>                      ! one center terms and, on the plane wave grid.
>                         IF (ASSOCIATED(H)) THEN
>                            CALL FOCK_CHARGE_ONE_CENTER_NOINT( WB(NBAA),WI(MNBI),&
>                              GWORK(1),H,CRHOIB(1,NBAA), CRHOLM,  SIZE(CRHOLM))
>                         ELSE
>                            CALL FOCK_CHARGE_NOINT( WB(NBAA),WI(MNBI), GWORK(1), &
>                             CRHOLM, SIZE(CRHOLM))
>                         ENDIF
>                         ! Set phase e^iGr, where G = k_i - k_q - k_a
>                         IF (LPHASE) THEN
>                         CALL APPLY_PHASE( GRIDHF, CPHASE(1), GWORK(1), &
>                           GWORK(1) )
>                            !IF (KQ==1) WRITE(*,*)'error: no apply_phase need for kq=1'
>                         ENDIF  
> 
>                      ! FFT{rho} to reciprocal space
>                         CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), &
>                           GWORK(1),GCHGIA(1,NBAA,2),WGWQ%GRID,.FALSE.)
>                         NFFT=NFFT+1
>                      ! multiply with potential factor                        
>                      ENDDO !NBAA (loop over all bands in NSTRIP)
>                      
>                      IF (ASSOCIATED(H)) THEN
>                         CALL APPLY_PHASE_ONE_CENTER(WHF%WDES, H, CRHOIB(:,:), &
>                           WHF%WDES%VKPT(:,KB)-WHF%WDES%VKPT(:,KI))
>                      ENDIF
>                      
>                      DO NBAA=1,NSTRIPA
>                         !copy FTOD functions to FTOD_PW and FTOD_OC
>                         
>                         tmp_FTOD_PW(1:NP,NBAA+NBA-1,1,2)=(GCHGIA(1:NP,NBAA,2)*(1.0_q/GRIDHF%NPLWV))
>                         IF (ASSOCIATED(H)) THEN
>                            tmp_FTOD_OC(1:H%TOTAL_ENTRIES,NBAA+NBA-1,1,2)=CRHOIB(1:H%TOTAL_ENTRIES,NBAA)
>                         ENDIF   
> 
>                         
>                      ENDDO
> !                  ENDDO !NBI (loop over valence bands only)
>                ENDDO !NBA (loop over all bands) 
>                
1694a2480,2595
>                CALL REDISTRIBUTE_FTOD_GRID_LOWMEM(WDES,KI,KQ,ISP,tmp_FTOD_PW,tmp_FTOD_OC)
> 
> !            ENDDO kiloop
>          
>             DEALLOCATE(GCHGIA,CRHOLM)
>          
> !         ENDDO kqloop
> !         ENDDO spin
> 
>          IF ((NBI==-1)) THEN
> 
>             IF (ALLOCATED(CRHOIA)) DEALLOCATE(CRHOIA)
>             IF (ALLOCATED(CRHOIB)) DEALLOCATE(CRHOIB)
>             IF (ALLOCATED(tmp_FTOD_OC)) DEALLOCATE(tmp_FTOD_OC)
>             IF (ALLOCATED(tmp_FTOD_PW)) DEALLOCATE(tmp_FTOD_PW)
>          
>             DO NBI=1,MAX(VBMAX(1),VBMAX(WDES%ISPIN))
>                CALL DELWAV(WI(NBI),.TRUE.)
>             ENDDO
>             DO NBA=1,NSTRIP
>                CALL DELWAV(WA(NBA),.TRUE.)
>             ENDDO
>             DO NBA=1,NSTRIP
>                CALL DELWAV(WB(NBA),.TRUE.)
>             ENDDO
>             DEALLOCATE(WI,WA,WB)
> 
>             RETURN         
> 
>          ELSE
> 
>             RETURN
> 
>          ENDIF
>          
>       END SUBROUTINE CALC_2ORBITAL_FTOD_LOWMEM
>  
>       
>       !***********************************************************************
>       !This subroutine redistributes the fourier-transformed overlap integrals <i|G|a>
>       !from the column-only process grid (CONTXT) to the quadratic 
>       !process grid (CONTXT_GRID)
>       !***********************************************************************
> 
>       SUBROUTINE REDISTRIBUTE_FTOD_GRID_LOWMEM(WDES,KI,KQ,ISP,tmp_FTOD_PW,tmp_FTOD_OC)
>          IMPLICIT NONE
>          TYPE(wavedes) WDES
>          INTEGER :: KI,KQ,ISP
>          GDEF, TARGET :: tmp_FTOD_OC(:,:,:,:)
>          COMPLEX(q), TARGET :: tmp_FTOD_PW(:,:,:,:)
>          INTEGER :: FTOD_PW_rows, FTOD_PW_cols,cc,FTOD_OC_rows, FTOD_OC_cols
>          INTEGER :: NBI,FTOD_PW_rows_br,VBMAX_tmp,MNBI
>          
>          ! Prepare array descriptors for ScaLAPACK
>          call BLACS_GRIDINFO(CONTXT_GRID, NPROW, NPCOL, MYROW, MYCOL)
>          MB=MIN(((NGVECTOR)/NPROW),50)   !Row blocking size
>          FTOD_PW_rows = numroc(NGVECTOR,mb,MYROW,0,NPROW)
>          desc_FTOD_PW(3) = NGVECTOR       ! global number of rows
>          desc_FTOD_PW(5) = mb       ! row block size
>          desc_FTOD_PW(9) = MAX(1,FTOD_PW_rows) ! leading dimension of local array      
> 
>          desc_FTOD_PW_br(1) = 1              ! descriptor type
>          desc_FTOD_PW_br(2) = contxt         ! blacs context
>          desc_FTOD_PW_br(3) = NGVECTOR    ! global number of rows
>          desc_FTOD_PW_br(4) = (PROCS*WDES%NBANDS) ! global number of cols
>          desc_FTOD_PW_br(5) = NGVECTOR     ! row block size
>          desc_FTOD_PW_br(6) = 1             ! col block size
>          desc_FTOD_PW_br(7) = 0              ! initial process row
>          desc_FTOD_PW_br(8) = 0              ! initial process col
>          desc_FTOD_PW_br(9) = MAX(1,(NGVECTOR)) ! leading dimension of local array
>          !Distribute the fourier-transformed overlap integrals
>          VBMAX_tmp=MAX(VBMAX(1),VBMAX(WDES%ISPIN))
>          DO cc=1,ncc
>             CALL PZGEMR2D(NGVECTOR,(PROCS*WDES%NBANDS),tmp_FTOD_PW(1,1,1,cc),1,1,&
>              desc_FTOD_PW_br,FTOD_PW(1,1,1,1,1,1,cc),1,1,desc_FTOD_PW,contxt_grid)
>          ENDDO
>          
> #ifdef gammareal         
>          desc_FTOD_PW(9) = MAX(1,m_ desc_FTOD_PW(9)) ! leading dimension of local array
>          desc_FTOD_PW(5) = m_ desc_FTOD_PW(5)
>          desc_FTOD_PW(3) = m_ desc_FTOD_PW(3)
> #else
> 
> #endif
>          !ONE-center part FTOD_OC
>            
>          IF (ASSOCIATED(H)) THEN
>          
>             ! Prepare array descriptors for ScaLAPACK
>          
>             desc_FTOD_OC_br(1) = 1            ! descriptor type
>             desc_FTOD_OC_br(2) = contxt       ! blacs context
>             desc_FTOD_OC_br(3) = NHVECTOR     ! global number of rows
>             desc_FTOD_OC_br(4) = (PROCS*WDES%NBANDS) ! global number of cols
>             desc_FTOD_OC_br(5) = NHVECTOR     ! row block size
>             desc_FTOD_OC_br(6) = 1            ! col block size
>             desc_FTOD_OC_br(7) = 0            ! initial process row
>             desc_FTOD_OC_br(8) = 0            ! initial process col
>             desc_FTOD_OC_br(9) = MAX(1,(NHVECTOR)) ! leading dimension of local array
>      
>  
>             DO cc=1,2
> #ifdef gammareal
>                CALL PDGEMR2D(NHVECTOR,(PROCS*WDES%NBANDS),tmp_FTOD_OC(1,1,1,cc),1,1,&
>                 desc_FTOD_OC_br,FTOD_OC(1,1,1,1,1,1,cc),1,1,desc_FTOD_OC,contxt_grid)
> #else
>                CALL PZGEMR2D(NHVECTOR,(PROCS*WDES%NBANDS),tmp_FTOD_OC(1,1,NBI,cc),1,1,&
>                 desc_FTOD_OC_br,FTOD_OC(1,1,1,1,1,1,cc),1,1,desc_FTOD_OC,contxt_grid)
> #endif
>             ENDDO           
> 
>          ENDIF
>          
>       END SUBROUTINE REDISTRIBUTE_FTOD_GRID_LOWMEM
>       
> #endif // scaLAPACK
EOF

md5sum_orig=6854ef37fdadf989aa4ade1d4be01cd5
md5sum_new=`md5sum ump2no.F  | awk '{print $1}' `
if [ "$md5sum_orig" = "$md5sum_new" ]; then
    echo "Patching file  ump2no.F " 
    patch ump2no.F .tmp_vaspcc4s_patch 
else
   echo "You dont have the correct version of file  ump2no.F  to apply patch. " 
fi
 
#######################
# diff of   .objects
#######################
 
 
cat > .tmp_vaspcc4s_patch<<"EOF"
217a218,220
> 	bracketst.o \
> 	ccsd.o \
> 	cc4s.o \
EOF

md5sum_orig="acee6f6052e79e4c14ac7bf131f7d410"
md5sum_new=`md5sum .objects  | awk '{print $1}' `
if [ "$md5sum_orig" = "$md5sum_new" ]; then
    echo "Patching file  .objects " 
    patch .objects .tmp_vaspcc4s_patch 
else
   echo "You dont have the correct version of file  .objects  to apply patch. " 
fi
 
#######################
# new file   bracketst.F
#######################
cat > .tmp_vaspcc4s_patch<<"EOF"
#include "symbol.inc"
!*********************************************************************************
!
! MODULE :: (T)
! AUTHOR :: Andreas Grüneis ( andreas.grueneis@tuwien.ac.at )
!
! DESCRIPTION: 
!
! THIS MODULE PERFORMS PERTURBATIVE TRIPLES COUPLED CLUSTER THEORY CALCULATIONS AND IS PART OF VASP
! 
! ================================================
! THE FOLLOWING SETTINGS ARE SUPPORTED:
! ================================================
!
!   * RESTRICTED (WDES%ISPIN=1) CANONICAL ORBITALS ONLY
!   * GAMMA-CENTERED OR SHIFTED K-MESHES
!   * FOR THE SPECIAL CASE OF GAMMA-CENTERED 1X1X1 AND 2X2X2 K-MESHES REAL ORBITALS (LORBITALREAL=.TRUE.)
!      CAN BE EMPLOYED. THIS REDUCES THE MEMORY FOOTPRINT SIGNIFICANTLY.
!
!
!*********************************************************************************

MODULE bracketst
      USE prec
      USE fock
      USE chi_base      
      USE lattice
      USE wpot
      IMPLICIT NONE

#ifdef scaLAPACK

!**********************************************************************
!
!**********************************************************************
      COMPLEX(q) :: E_MP2
      ! important intermediates for real orbitals

      INTEGER :: NGVECTOR, NHVECTOR
      integer, ALLOCATABLE :: KPTS_MKPTS(:)
      integer, ALLOCATABLE :: PROCS_KPTS(:)
      integer, ALLOCATABLE :: MKPTS_KPTS(:)
      integer, PUBLIC :: MY_NKPTS, MAX_MY_NKPTS
      INTEGER, PUBLIC  :: CONTXT_COLS,CONTXT_GRID
      INTEGER, PUBLIC :: NPROW_GRID, ncc
      INTEGER, PUBLIC :: VBMAX, NUNOCC
      INTEGER, PUBLIC :: PROCS, ME
      INTEGER, PRIVATE :: NPROW,NPCOL,MYROW,MYCOL,NBLOCKAB
      INTEGER :: TETS, TET2S, TETSORT, iteration
      !BLACS related variables and function
      !PROCS.. number of processors, ME... processor number, NPROW... number of rows in process grid
      !NPCOL... number of columns in process grid, myrow,mycol... my coordinates in process grid
      !mb... blocking size of rows for block cyclic distribution
      !nb... blocking size of columns for block cyclic distribution      
      
      REAL(q) :: NFLOAT4O,tetfpo,tetcom,tetr2d,tetvv,tetcc,tetvc,tetcontv
      INTEGER, EXTERNAL :: NUMROC, BLACS_PNUM
      INTEGER :: NTHREADS

      TYPE(one_center_handle), POINTER, PUBLIC, SAVE :: H
      LOGICAL :: LORBREAL
      LOGICAL :: LUSESHARED

!!!!! FOR TRIPLES
      COMPLEX(qs), ALLOCATABLE :: W_ijk_abc(:,:,:,:,:,:)
      REAL(qs), ALLOCATABLE :: W_ijk_abc_R(:,:,:,:,:,:)
      COMPLEX(qs), ALLOCATABLE :: Wint_ijk_abc(:,:,:,:,:,:)
      REAL(qs), ALLOCATABLE :: Wint_ijk_abc_R(:,:,:,:,:,:)
      COMPLEX(qs), ALLOCATABLE :: Wtmp_ijk_abc(:,:,:,:,:,:)
      REAL(qs), ALLOCATABLE :: Wtmp_ijk_abc_R(:,:,:,:,:,:)
      COMPLEX(qs), ALLOCATABLE :: T2_MTMP(:,:,:,:)
      REAL(qs), ALLOCATABLE :: T2_MTMP_R(:,:,:,:)
     ! T2_TMP and PW_TMP
      COMPLEX(qs), ALLOCATABLE :: T2_TMP(:,:,:,:,:,:,:)
      REAL(qs), ALLOCATABLE :: T2_TMP_R(:,:,:,:,:,:,:)
      COMPLEX(qs), POINTER :: T2(:,:,:,:,:,:)
      INTEGER :: T2_SHMID
      COMPLEX(qs), ALLOCATABLE :: T1(:,:,:)
      REAL(qs), POINTER :: T2_R(:,:,:,:,:,:)
      INTEGER :: T2_R_SHMID
      REAL(qs), ALLOCATABLE :: T1_R(:,:,:)
      COMPLEX(q), POINTER :: FTOD_PW_IJ(:,:,:,:,:)
      COMPLEX(q), POINTER :: FTOD_OC_IJ(:,:,:,:,:)
      COMPLEX(q), POINTER :: FTOD_PW_AB(:,:,:,:), FTOD_PW_AI(:,:,:,:,:)
      COMPLEX(q), POINTER :: FTOD_OC_AB(:,:,:,:), FTOD_OC_AI(:,:,:,:,:)
      INTEGER :: FTOD_PW_AB_SHMID
      INTEGER :: FTOD_OC_AB_SHMID
      INTEGER :: FTOD_PW_AI_SHMID
      INTEGER :: FTOD_OC_AI_SHMID
      COMPLEX(q), ALLOCATABLE :: PW_AB_TMP(:,:,:,:,:)
      COMPLEX(q), ALLOCATABLE :: OC_AB_TMP(:,:,:,:,:)
      COMPLEX(q), ALLOCATABLE :: PW_AI_TMP(:,:,:,:,:,:)
      COMPLEX(q), ALLOCATABLE :: OC_AI_TMP(:,:,:,:,:,:)
      ! indexing for T2_KA_KJ
      INTEGER, ALLOCATABLE :: KJ_INDEX(:), KA_INDEX(:), PROCS_KA_KJ(:,:)
      INTEGER, ALLOCATABLE :: INDEX_KA_KJ(:,:)
      INTEGER, ALLOCATABLE :: INDEX_KQ_KA(:,:)
      ! indexing for FTOD_AI/AB and 
      INTEGER, ALLOCATABLE :: PROCS_KQ_KA(:,:)
      ! indexing for T2_TMP
      INTEGER, ALLOCATABLE :: KP_TMP(:), NTMP_KP(:)
      INTEGER :: NKINDEX
      INTEGER :: EMSCALCWINT
      COMPLEX(q) :: ETRIPLES, ETRIPLESZ
!!!!!!!!!!
      LOGICAL :: SHIFTED_KPOINTS
      REAL(q) :: WTKPT
      INTEGER :: REALNKPTS
      INTEGER, ALLOCATABLE :: RKQofKQ(:), RKIofKI(:)
      INTEGER, ALLOCATABLE :: KQofRKQ(:), KIofRKI(:)
!!!!!!!!!!
      LOGICAL :: LSFACTOR
      COMPLEX(q), ALLOCATABLE :: SF_ETRIPLES(:,:), SF_ETRIPLESZ(:,:)
      COMPLEX(qs), ALLOCATABLE :: Wint_ijk_abc_kq(:,:,:,:,:,:,:)
      COMPLEX(qs), ALLOCATABLE :: W_ijk_abc_kq(:,:,:,:,:,:,:)
      REAL(q), ALLOCATABLE :: POTFAK_FULL(:,:), GVECLEN(:,:), GVECX(:,:), GVECY(:,:), GVECZ(:,:)

      
      CONTAINS

!***********************************************************************
!***********************************************************************
#ifdef gammareal
      SUBROUTINE CALCULATE_bracketsT(P,WDES,W,LATT_CUR,T_INFO,IO,KPOINTS,WGW, LMAXMP2)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE base
         USE mkpoints
         USE mpimy
         USE ini
         USE choleski
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavedes) WGW
         TYPE(wavespin) W
         TYPE(type_info) T_INFO
         TYPE(potcar) P(T_INFO%NTYP)
         TYPE(latt) LATT_CUR
         INTEGER LMAXMP2,NBANDSGW
         REAL(q) :: ENCUTGW, ENCUTGWSOFT,change,denom,AMIX,BMIX
         TYPE (in_struct) IO
         TYPE (kpoints_struct) KPOINTS
         INTEGER :: NI,NJ,NR,NS,KI,KJ,KR,KS,KQ,KQ_,direction,KR_,NA,MKA,MKJ, &
                    GKQ, MKB, KC
         integer :: TWOE4ORBITAL_COLS,TWOE4ORBITAL_ROWS,niteration,nwfup
         integer :: time_array1(8),time_array2(8),ems,KN,NN,omp_get_num_threads

      END SUBROUTINE CALCULATE_bracketsT
#else
      SUBROUTINE CALCULATE_bracketsT(P,WDES,W,LATT_CUR,T_INFO,IO,KPOINTS,WGW, LMAXMP2)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE base
         USE mkpoints
         USE mpimy
         USE ini
         USE choleski
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavedes) WGW
         TYPE(wavespin) W        
         TYPE(type_info) T_INFO
         TYPE(potcar) P(T_INFO%NTYP)
         TYPE(latt) LATT_CUR
         INTEGER LMAXMP2,NBANDSGW
         TYPE (in_struct) IO
         TYPE (kpoints_struct) KPOINTS
         INTEGER :: NI,NJ,NR,NS,KI,KJ,KR,KS,KQ,KQ_,direction,KR_,NA,MKA,MKJ, &
                    GKQ, MKB, KC
         integer :: time_array1(8),time_array2(8),ems,KN,NN,omp_get_num_threads
         TYPE(wavedes1) WGWQ

         ncc=2

         LUSESHARED=.FALSE.
         LORBREAL=.FALSE.
         IF (WDES%LORBITALREAL) LORBREAL=.TRUE.

         CALL CHECK_FULL_KPOINTS ! all set up properly ?
         CALL CHECK_SHIFTED_KPOINTS(WDES,W,KPOINTS) 
         IF ((IO%IU0>0) .and. (SHIFTED_KPOINTS)) write(IO%IU0,*)'You are using a shifted k-mesh.'
 
         CALL INIT_BLACS_COLS()

         CALL SETUP_T2_KINDEX(WDES)
         CALL IN_T2(IO, WDES, W, WGW)
         IF (IO%IU0>=0) write(IO%IU0,*)'T2CAR file read in'
         CALL IN_T1(IO, WDES, W)
         IF (IO%IU0>=0) write(IO%IU0,*)'T1CAR file read in'
         CALL IN_FTOD_PW(IO, WDES, W, WGW)
         IF (IO%IU0>=0) write(IO%IU0,*)'FTODCAR file read in'


         LSFACTOR=.FALSE.
         IF (LSFACTOR) THEN
            IF (.not. ALLOCATED(GVECLEN)) THEN
               ALLOCATE(GVECLEN(NGVECTOR,REALNKPTS))
               ALLOCATE(GVECX(NGVECTOR,REALNKPTS))
               ALLOCATE(GVECY(NGVECTOR,REALNKPTS))
               ALLOCATE(GVECZ(NGVECTOR,REALNKPTS))
               GVECLEN=0.0_q
               GVECX=0.0_q
               GVECY=0.0_q
               GVECZ=0.0_q
            ENDIF
            DO KQ=1,WDES%NKPTS
               CALL SETWDES(WGW,WGWQ,KQ)
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0.00_q)) CYCLE
               CALL SET_GVEC(WGWQ, LATT_CUR, GVECLEN(1,RKQofKQ(KQ)), GVECX(1,RKQofKQ(KQ)), GVECY(1,RKQofKQ(KQ)), GVECZ(1,RKQofKQ(KQ)))
            ENDDO
            ALLOCATE(SF_ETRIPLES(NGVECTOR,REALNKPTS),SF_ETRIPLESZ(NGVECTOR,REALNKPTS))
            SF_ETRIPLES=zero
            SF_ETRIPLESZ=zero
         ENDIF

         CALL CALC_TRIPLES(WDES,WGW,W,ITERATION,KPOINTS,IO)
         IF (IO%IU0>=0) WRITE(IO%IU6,*) '[T] contribution to the correlation energy is',ETRIPLES*REAL(1.0_q/REALNKPTS/REALNKPTS/REALNKPTS/REALNKPTS,kind=q)
         IF (IO%IU0>=0) WRITE(IO%IU6,*) '(T) contribution to the correlation energy is',(ETRIPLES+ETRIPLESZ)*REAL(1.0_q/REALNKPTS/REALNKPTS/REALNKPTS/REALNKPTS,kind=q)
         
         call BLACS_GRIDEXIT(contxt_cols)
         
      END SUBROUTINE CALCULATE_bracketsT


! this routine reads in the T2CAR file

     SUBROUTINE IN_T2(IO, WDES, W, WGW)
      USE base
      USE wave
      USE lattice
      USE main_mpi
      IMPLICIT NONE
      TYPE (in_struct)   IO
      TYPE (wavedes)  WDES
      TYPE (wavespin) W
      TYPE(wavedes) WGW
    ! local
      INTEGER NPL_TOT, IU0, IRECLW, MRECL
      INTEGER*8 IREC
      COMPLEX(qs) :: CNUNOCC, CVBMAX, CNKPTS
      REAL(qs) :: RNUNOCC, RVBMAX, RNKPTS
      REAL(q) :: MVKPT
      INTEGER :: NI,NJ,NA,NB,KI,KA,KJ, START_IT
      LOGICAL :: file_exists
      CHARACTER ch
      
      INQUIRE(FILE="T2CAR", EXIST=file_exists)
      IF (file_exists) THEN
      IF (ME==0) WRITE(*,*)'Found T2CAR file. Reading in T2 amplitudes.'
      START_IT=2

      NUNOCC=0
      VBMAX=0


      CNUNOCC=NUNOCC
      RNUNOCC=NUNOCC
      CVBMAX=VBMAX
      RVBMAX=VBMAX
      CNKPTS=REALNKPTS
      RNKPTS=REALNKPTS
      call BLACS_PINFO(ME,PROCS)
      MRECL=8
      IF (ME==0) THEN
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'T2CAR',ACCESS='SEQUENTIAL',STATUS='UNKNOWN',FORM='UNFORMATTED')

         !read info about number of orbitals and k-points from header of T2CAR
          READ(12)CNUNOCC
          NUNOCC=CNUNOCC
          READ(12)CVBMAX
          VBMAX=CVBMAX
          READ(12)CNKPTS
          IREC=4
          IF (REALNKPTS/=CNKPTS) THEN
             WRITE(*,*)'Number of k-points in T2CAR and KPOINTS files differ.'
             CALL EXIT
          ENDIF
          DO KA=1,WDES%NKPTS
             IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0)) CYCLE
             READ(12)MVKPT
             WDES%VKPT(1,KA)=MVKPT
             IREC=IREC+1
             READ(12)MVKPT
             WDES%VKPT(2,KA)=MVKPT
             IREC=IREC+1
             READ(12)MVKPT
             WDES%VKPT(3,KA)=MVKPT
             IREC=IREC+1

             DO NA=1,VBMAX+NUNOCC
               READ(12)MVKPT
               IREC=IREC+1
             ENDDO
          ENDDO
      ENDIF
    
      CALLMPI( M_sum_i(WGW%COMM_INTER, VBMAX, 1)) 
      CALLMPI( M_sum_i(WGW%COMM_INTER, NUNOCC, 1)) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE(T2_MTMP(NUNOCC,NUNOCC,VBMAX,VBMAX)) 
      IF (LORBREAL) THEN
         ALLOCATE(T2_R(NUNOCC,NUNOCC,VBMAX,VBMAX,REALNKPTS,NKINDEX))
      ELSE
         ALLOCATE(T2(NUNOCC,NUNOCC,VBMAX,VBMAX,REALNKPTS,NKINDEX)) 
      ENDIF
      IF (LORBREAL) ALLOCATE(T2_MTMP_R(NUNOCC,NUNOCC,VBMAX,VBMAX)) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !now communicate and write the T2 amplitudes
      DO KA=1,WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0)) CYCLE
      DO KJ=1,WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ)==0)) CYCLE
      DO KI=1,WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI)==0)) CYCLE

         IF (.not. LUSESHARED) THEN

         IF (ME==0) THEN
            DO NJ=1,VBMAX
            DO NI=1,VBMAX
            DO NB=1,NUNOCC
            DO NA=1,NUNOCC
               READ(12) T2_MTMP(NA,NB,NI,NJ)
               IREC=IREC+1
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDIF

         IF (LORBREAL) T2_MTMP_R=REAL(T2_MTMP,kind=qs)

         IF (LORBREAL) THEN
            IF (ME==0) THEN
               CALL SGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC),(VBMAX*VBMAX), &
                  T2_MTMP_R(1,1,1,1), (NUNOCC*NUNOCC))
               IF (PROCS_KA_KJ(KA,KJ)==ME) T2_R(:,:,:,:,RKIofKI(KI),INDEX_KA_KJ(KA,KJ))=T2_MTMP_R
            ELSE
               CALL SGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC),(VBMAX*VBMAX), &
                  T2_MTMP_R(1,1,1,1), (NUNOCC*NUNOCC),0,0)
               IF (PROCS_KA_KJ(KA,KJ)==ME) T2_R(:,:,:,:,RKIofKI(KI),INDEX_KA_KJ(KA,KJ))=T2_MTMP_R
            ENDIF
         ELSE 
            IF (ME==0) THEN
               CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC),(VBMAX*VBMAX), &
                  T2_MTMP(1,1,1,1), (NUNOCC*NUNOCC))
               IF (PROCS_KA_KJ(KA,KJ)==ME) T2(:,:,:,:,RKIofKI(KI),INDEX_KA_KJ(KA,KJ))=T2_MTMP(:,:,:,:)
            ELSE
               CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC),(VBMAX*VBMAX), &
                  T2_MTMP(1,1,1,1), (NUNOCC*NUNOCC),0,0)
               IF (PROCS_KA_KJ(KA,KJ)==ME) T2(:,:,:,:,RKIofKI(KI),INDEX_KA_KJ(KA,KJ))=T2_MTMP(:,:,:,:)
            ENDIF
         ENDIF

         ELSE ! if shared memory is used
           !write(*,*)'reading T2 amplitudes',me
           IF (ME==0) THEN
               DO NJ=1,VBMAX
               DO NI=1,VBMAX
               DO NB=1,NUNOCC
               DO NA=1,NUNOCC
                  READ(12) T2_MTMP(NA,NB,NI,NJ)
                  IREC=IREC+1
               ENDDO
               ENDDO
               ENDDO
               ENDDO
            ENDIF

            IF (LORBREAL) T2_MTMP_R=REAL(T2_MTMP,kind=qs)

            IF (LORBREAL) THEN
               IF (ME==0) THEN
                  CALL SGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC),(VBMAX*VBMAX), &
                     T2_MTMP_R(1,1,1,1), (NUNOCC*NUNOCC))
                  T2_R(:,:,:,:,RKIofKI(KI),INDEX_KA_KJ(KA,KJ))=T2_MTMP_R
               ELSE
                  CALL SGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC),(VBMAX*VBMAX), &
                     T2_MTMP_R(1,1,1,1), (NUNOCC*NUNOCC),0,0)
                  T2_R(:,:,:,:,RKIofKI(KI),INDEX_KA_KJ(KA,KJ))=T2_MTMP_R
               ENDIF
            ELSE 
               IF (ME==0) THEN
                  !write(*,*)'sending T2 amplitudes',me
                  CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC),(VBMAX*VBMAX), &
                     T2_MTMP(1,1,1,1), (NUNOCC*NUNOCC))
                  T2(:,:,:,:,RKIofKI(KI),INDEX_KA_KJ(KA,KJ))=T2_MTMP(:,:,:,:)
               ELSE
                  !write(*,*)'receiving T2 amplitudes',me
                  CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC),(VBMAX*VBMAX), &
                     T2_MTMP(1,1,1,1), (NUNOCC*NUNOCC),0,0)
                  T2(:,:,:,:,RKIofKI(KI),INDEX_KA_KJ(KA,KJ))=T2_MTMP(:,:,:,:)
               ENDIF
            ENDIF


         ENDIF
      ENDDO
      ENDDO
      ENDDO

      IF (ME==0) CLOSE(12)

      ENDIF

      DEALLOCATE(T2_MTMP)
      IF (ALLOCATED(T2_MTMP_R)) DEALLOCATE(T2_MTMP_R)

      RETURN

      write(*,*)' T2 done',me
     END SUBROUTINE IN_T2

     SUBROUTINE CLOSET2
        CLOSE(12)
     END SUBROUTINE CLOSET2

! this routine reads in the T1CAR file

     SUBROUTINE IN_T1(IO, WDES, W)
      USE base
      USE wave
      USE lattice
      USE main_mpi
      IMPLICIT NONE
      TYPE (in_struct)   IO
      TYPE (wavedes)  WDES
      TYPE (wavespin) W
    ! local
      INTEGER NPL_TOT, IU0, IRECLW, MRECL
      INTEGER*8 IREC
      COMPLEX(qs) :: CNUNOCC, CVBMAX, CNKPTS
      REAL(qs) :: RNUNOCC, RVBMAX, RNKPTS
      REAL(q) :: MVKPT
      REAL(qs), ALLOCATABLE :: T1_MTMP_R(:,:,:)
      INTEGER :: NI,NJ,NA,NB,KI,KA,KJ, START_IT
      LOGICAL :: file_exists
      
      INQUIRE(FILE="T1CAR", EXIST=file_exists)
      IF (file_exists) THEN

      IF (ME==0) WRITE(*,*)'Found T1CAR file. Reading in T1 amplitudes.'
     
      CNUNOCC=NUNOCC
      RNUNOCC=NUNOCC
      CVBMAX=VBMAX
      RVBMAX=VBMAX
      CNKPTS=REALNKPTS
      RNKPTS=REALNKPTS

      call BLACS_PINFO(ME,PROCS)
      ! RECL : the T2 amplitudes are stored in single precision complex
      MRECL=8
      IF (ME==0) THEN
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'T1CAR',ACCESS='DIRECT', &
             FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=MRECL)

         !write info about number of orbitals and k-points to header of T2CAR
          READ(12,REC=1)CNUNOCC
          NUNOCC=CNUNOCC
          READ(12,REC=2)CVBMAX
          VBMAX=CVBMAX
          READ(12,REC=3)CNKPTS
          IREC=4
          IF (REALNKPTS/=CNKPTS) THEN
             WRITE(*,*)'Number of k-points in T1CAR and KPOINTS files differ.'
             CALL EXIT
          ENDIF
          DO KA=1,WDES%NKPTS
             IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0)) CYCLE
             READ(12,REC=IREC)MVKPT
             WDES%VKPT(1,KA)=MVKPT
             IREC=IREC+1
             READ(12,REC=IREC)MVKPT
             WDES%VKPT(2,KA)=MVKPT
             IREC=IREC+1
             READ(12,REC=IREC)MVKPT
             WDES%VKPT(3,KA)=MVKPT
             IREC=IREC+1
          ENDDO
      ENDIF

    
      ALLOCATE(T1(NUNOCC,VBMAX,REALNKPTS))
      IF (LORBREAL) THEN
         ALLOCATE(T1_R(NUNOCC,VBMAX,REALNKPTS)) 
      ENDIF

      !now communicate and write the T2 amplitudes

      IF (ME==0) THEN
         DO KI=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI)==0)) CYCLE
         DO NI=1,VBMAX
         DO NA=1,NUNOCC
            READ(12,REC=IREC) T1(NA,NI,RKIofKI(KI))
            IREC=IREC+1
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (LORBREAL) T1_R=REAL(T1,kind=qs)

      IF (LORBREAL) THEN
         IF (ME==0) THEN
            CALL SGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC),(VBMAX*REALNKPTS), &
               T1_R(1,1,1), (NUNOCC))
         ELSE
            CALL SGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC),(VBMAX*REALNKPTS), &
               T1_R(1,1,1), (NUNOCC),0,0)
         ENDIF
      ELSE
         IF (ME==0) THEN
            CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC),(VBMAX*REALNKPTS), &
               T1(1,1,1), (NUNOCC))
         ELSE
            CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC),(VBMAX*REALNKPTS), &
               T1(1,1,1), (NUNOCC),0,0)
         ENDIF
      ENDIF


      IF (ME==0) CLOSE(12)

      
      ENDIF

     END SUBROUTINE IN_T1


!! this routine writes out the overlapdensities that are needed in the calculations of (T)-energies to FTODCAR

     SUBROUTINE IN_FTOD_PW(IO, WDES, W, WGW)
      USE base
      USE wave
      USE lattice
      USE main_mpi
      IMPLICIT NONE
      TYPE (in_struct)   IO
      TYPE (wavedes)  WDES
      TYPE (wavespin) W
      TYPE(wavedes) WGW
    ! local
      INTEGER NPL_TOT, IU0, IRECLW, MRECL
      INTEGER*8 IREC
      COMPLEX(qs) :: CNUNOCC, CVBMAX, CNKPTS, CNGVECTOR
      REAL(qs) :: RNUNOCC, RVBMAX, RNKPTS
      COMPLEX(q), ALLOCATABLE :: PW_IJ_TMP(:,:,:,:)
      INTEGER :: NI,NJ,NA,NB,KI,KA,KJ,NG,KQ,mncc,ierror
      LOGICAL :: file_exists
      
      INQUIRE(FILE="FTODCAR", EXIST=file_exists)

      IF ((ME==0) .and. (file_exists)) WRITE(*,*)'Found FTODCAR file. Reading in FTODs.'
      NGVECTOR=0 
      VBMAX=0
      NUNOCC=0
      IF (.not. file_exists) THEN
         WRITE(*,*)'FTODCAR file not found. Exiting.'
         CALL EXIT()
      ENDIF

      call BLACS_PINFO(ME,PROCS)
      ! RECL : the FTODs are stored in double precision complex
      MRECL=16
      IF (ME==0) THEN
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'FTODCAR',ACCESS='DIRECT', &
             FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=MRECL)

         !write info about number of orbitals and k-points to header of T2CAR
          READ(12,REC=1)CNUNOCC
          READ(12,REC=2)CVBMAX
          READ(12,REC=3)CNKPTS
          VBMAX=CVBMAX
          NUNOCC=CNUNOCC

          READ(12,REC=4)CNGVECTOR

          IREC=5
          DO KA=1,WDES%NKPTS
             IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0)) CYCLE
             READ(12,REC=IREC)WDES%VKPT(1,KA)
             READ(12,REC=IREC+1)WDES%VKPT(2,KA)
             READ(12,REC=IREC+2)WDES%VKPT(3,KA)
             IREC=IREC+3
          ENDDO
          NGVECTOR=CNGVECTOR
         IF (VBMAX/=CVBMAX) WRITE(*,*)'SEVERE ERROR! Check if the FTODCAR, T1CAR and T2CAR files correspond to the same system.'
         IF (NUNOCC/=CNUNOCC) WRITE(*,*)'SEVERE ERROR! Check if the FTODCAR, T1CAR and T2CAR files correspond to the same system.'
      ENDIF

      CALLMPI( M_sum_i(WGW%COMM_INTER, NGVECTOR, 1)) 
      CALLMPI( M_sum_i(WGW%COMM_INTER, VBMAX, 1)) 
      CALLMPI( M_sum_i(WGW%COMM_INTER, NUNOCC, 1)) 

      CALLMPI( MPI_barrier(WGW%COMM_INTER,ierror) ) 
      ALLOCATE(PW_AB_TMP(NGVECTOR,NUNOCC,NUNOCC,1,1))
      ALLOCATE(PW_AI_TMP(NGVECTOR,VBMAX,NUNOCC,1,1,1))
      ALLOCATE(PW_IJ_TMP(NGVECTOR,VBMAX,VBMAX,REALNKPTS))

      ALLOCATE(FTOD_PW_AB(NGVECTOR,NUNOCC,NUNOCC,NKINDEX))
      ALLOCATE(FTOD_PW_AI(NGVECTOR,VBMAX,NUNOCC,NKINDEX,2))
      ALLOCATE(FTOD_PW_IJ(NGVECTOR,VBMAX,VBMAX,REALNKPTS,REALNKPTS))

      !now read and communicate the FTOD_PW_AB amplitudes
      DO KA=1,WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0)) CYCLE
      DO KQ=1,WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0.0_q)) CYCLE

         IF (ME==0) THEN
            DO NB=1,NUNOCC
            DO NA=1,NUNOCC
            DO NG=1,NGVECTOR
               READ(12,REC=IREC) PW_AB_TMP(NG,NA,NB,1,1)
               IREC=IREC+1
            ENDDO
            ENDDO
            ENDDO
         ENDIF

         IF (.not. LUSESHARED) THEN

            IF (ME==0) THEN
               CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR)*(NUNOCC),(NUNOCC), &
                  PW_AB_TMP(1,1,1,1,1), (NGVECTOR*NUNOCC))
               IF (PROCS_KQ_KA(KQ,KA)==ME) FTOD_PW_AB(:,:,:,INDEX_KQ_KA(KQ,KA))=PW_AB_TMP(:,:,:,1,1)
            ELSE
               CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR)*(NUNOCC),(NUNOCC), &
                  PW_AB_TMP(1,1,1,1,1), (NGVECTOR*NUNOCC),0,0)
               IF (PROCS_KQ_KA(KQ,KA)==ME) FTOD_PW_AB(:,:,:,INDEX_KQ_KA(KQ,KA))=PW_AB_TMP(:,:,:,1,1)
            ENDIF

         ELSE

            IF (ME==0) THEN
               CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR)*(NUNOCC),(NUNOCC), &
                  PW_AB_TMP(1,1,1,1,1), (NGVECTOR*NUNOCC))
               FTOD_PW_AB(:,:,:,INDEX_KQ_KA(KQ,KA))=PW_AB_TMP(:,:,:,1,1)
            ELSE
               CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR)*(NUNOCC),(NUNOCC), &
                  PW_AB_TMP(1,1,1,1,1), (NGVECTOR*NUNOCC),0,0)
               FTOD_PW_AB(:,:,:,INDEX_KQ_KA(KQ,KA))=PW_AB_TMP(:,:,:,1,1)
            ENDIF

         ENDIF
      ENDDO
      ENDDO


      !now read and communicate the FTOD_PW_AI
      DO mncc=1,2
      DO KA=1,WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0)) CYCLE
      DO KQ=1,WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0.0_q)) CYCLE

         IF (ME==0) THEN
            DO NA=1,NUNOCC
            DO NI=1,VBMAX
            DO NG=1,NGVECTOR
               READ(12,REC=IREC) PW_AI_TMP(NG,NI,NA,1,1,1)
               IREC=IREC+1
            ENDDO
            ENDDO
            ENDDO
         ENDIF

         IF (.not. LUSESHARED) THEN
            IF (ME==0) THEN
               CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR)*(VBMAX),(NUNOCC), &
                  PW_AI_TMP(1,1,1,1,1,1), (NGVECTOR*VBMAX))
               IF (PROCS_KQ_KA(KQ,KA)==ME) FTOD_PW_AI(:,:,:,INDEX_KQ_KA(KQ,KA),mncc)=PW_AI_TMP(:,:,:,1,1,1)
            ELSE
               CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR)*(VBMAX),(NUNOCC), &
                  PW_AI_TMP(1,1,1,1,1,1), (NGVECTOR*VBMAX),0,0)
               IF (PROCS_KQ_KA(KQ,KA)==ME) FTOD_PW_AI(:,:,:,INDEX_KQ_KA(KQ,KA),mncc)=PW_AI_TMP(:,:,:,1,1,1)
            ENDIF
         ELSE
            IF (ME==0) THEN
               CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR)*(VBMAX),(NUNOCC), &
                  PW_AI_TMP(1,1,1,1,1,1), (NGVECTOR*VBMAX))
               FTOD_PW_AI(:,:,:,INDEX_KQ_KA(KQ,KA),mncc)=PW_AI_TMP(:,:,:,1,1,1)
            ELSE
               CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR)*(VBMAX),(NUNOCC), &
                  PW_AI_TMP(1,1,1,1,1,1), (NGVECTOR*VBMAX),0,0)
               FTOD_PW_AI(:,:,:,INDEX_KQ_KA(KQ,KA),mncc)=PW_AI_TMP(:,:,:,1,1,1)
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      !now read and communicate the FTOD_PW_IJ
      DO KI=1,WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI)==0)) CYCLE

         IF (ME==0) THEN
            DO KQ=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0.0_q)) CYCLE
            DO NI=1,VBMAX
            DO NJ=1,VBMAX
            DO NG=1,NGVECTOR
               READ(12,REC=IREC) PW_IJ_TMP(NG,NJ,NI,RKQofKQ(KQ))
               IREC=IREC+1
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDIF

         IF (ME==0) THEN
            CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR)*(VBMAX),(VBMAX*REALNKPTS), &
               PW_IJ_TMP(1,1,1,1), (NGVECTOR*VBMAX))
            FTOD_PW_IJ(:,:,:,:,RKIofKI(KI))=PW_IJ_TMP(:,:,:,:)
         ELSE
            CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR)*(VBMAX),(VBMAX*REALNKPTS), &
               PW_IJ_TMP(1,1,1,1), (NGVECTOR*VBMAX),0,0)
            FTOD_PW_IJ(:,:,:,:,RKIofKI(KI))=PW_IJ_TMP(:,:,:,:)
         ENDIF

      ENDDO

      IF (ME==0) CLOSE(12)

      IF (ALLOCATED(PW_IJ_TMP)) DEALLOCATE(PW_IJ_TMP)

     END SUBROUTINE IN_FTOD_PW


      !***********************************************************************
      !This routine initializes a process grid which is used for the 
      !pre-calculation of the fourier-transformed overlap integrals <i|G|a>
      !This process grid has NPROCS(number of processors used) columns and 1 row.
      !The assigned context variable is called CONTXT_COLS.
      !*********************************************************************** 

      SUBROUTINE INIT_BLACS_COLS()
         implicit none           
         INTEGER :: a_PRCS, i
         REAL :: NPCOL_TMP
         
         !first we create a column-only process grid in order to
         !calculate the <i|-G|a> and <j|G|b> quantities
         call BLACS_PINFO(ME,PROCS)
         nprow=1
         npcol=PROCS
         call BLACS_PINFO(ME,PROCS)
         call BLACS_GET     (0, 0, CONTXT_COLS)
         
         call BLACS_GRIDINIT(CONTXT_COLS, 'R', NPROW, NPCOL)
         call BLACS_GRIDINFO(CONTXT_COLS, NPROW, NPCOL, MYROW, MYCOL)
        
      END SUBROUTINE INIT_BLACS_COLS

      !***********************************************************************
      !This routine initializes the process grid, which is used for the
      !matrix-matrix multiplications and their
      !block-cyclic data distribution. Note that the context variable
      !CONTXT_GRID is used for this grid.
      !The routine tries to create the most quadratic grid possible for the given
      !number of processors.
      !*********************************************************************** 

      SUBROUTINE INIT_BLACS_GRID(WDES)
         implicit none           
         TYPE(wavedes) WDES
         
         call BLACS_PINFO(ME,PROCS)
         
         nprow=NPROW_GRID
         npcol=PROCS/NPROW
         IF (MOD(PROCS/NPROW,1)/=0) THEN
            WRITE(*,*)'internal error in INIT_BLACS_GRID: Bad process grid'
         ENDIF
         
         IF (NGVECTOR>(PROCS*WDES%NBANDS)) THEN
            IF (NPROW_GRID<(PROCS/NPROW)) THEN
               NPCOL=NPROW_GRID
               NPROW=PROCS/NPROW_GRID
               NPROW_GRID=NPROW
            ENDIF
         ENDIF         
         IF (NGVECTOR<(PROCS*WDES%NBANDS)) THEN
            IF (NPROW_GRID>(PROCS/NPROW)) THEN
               NPCOL=NPROW_GRID
               NPROW=PROCS/NPROW_GRID
               NPROW_GRID=NPROW
            ENDIF
         ENDIF
         IF ((MYROW==0) .AND. (MYCOL==0)) THEN
           WRITE(*,'(A,I3,A,I3,A)')'The allocated processors form a',NPROW_GRID,'x',NPCOL,' grid.'
         ENDIF

         call BLACS_PINFO(ME,PROCS)
         call BLACS_GET     (0, 0, CONTXT_GRID)
         call BLACS_GRIDINIT(CONTXT_GRID, 'R', NPROW, NPCOL)
         call BLACS_GRIDINFO(CONTXT_GRID, NPROW, NPCOL, MYROW, MYCOL)

      END SUBROUTINE INIT_BLACS_GRID

      SUBROUTINE SETUP_FTOD(WDES)
         IMPLICIT NONE
         TYPE(wavedes) WDES      

         call BLACS_GRIDINFO(CONTXT_GRID, NPROW, NPCOL, MYROW, MYCOL)
         !Blocking size for block-cyclic distribution of RESPF_PW          
         
      END SUBROUTINE SETUP_FTOD
      
      SUBROUTINE SETUP_T2_KINDEX(WDES)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE mkpoints
         USE lattice
         USE constant
         USE base
         IMPLICIT NONE
         INTEGER KA, KJ, KQ
         TYPE(wavedes) WDES
         INTEGER MKB(WDES%NKPTS*WDES%NKPTS), MKC(WDES%NKPTS*WDES%NKPTS), ATPROCS

         call BLACS_PINFO(ME,PROCS)
         IF (.not. LUSESHARED) THEN

         ATPROCS=0
         NKINDEX=0
         DO KJ=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ)==0)) CYCLE
         DO KA=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0)) CYCLE
            IF (ME==ATPROCS) THEN
               NKINDEX=NKINDEX+1
            ENDIF
            ATPROCS=ATPROCS+1
            IF (ATPROCS>(PROCS-1)) ATPROCS=0
         ENDDO
         ENDDO

         ALLOCATE(KJ_INDEX(NKINDEX))
         ALLOCATE(KA_INDEX(NKINDEX))
         ALLOCATE(PROCS_KA_KJ(WDES%NKPTS,WDES%NKPTS))
         ALLOCATE(PROCS_KQ_KA(WDES%NKPTS,WDES%NKPTS))
         ALLOCATE(INDEX_KA_KJ(WDES%NKPTS,WDES%NKPTS))
         ALLOCATE(INDEX_KQ_KA(WDES%NKPTS,WDES%NKPTS))
         INDEX_KA_KJ=0
         PROCS_KA_KJ=-1
         ATPROCS=0
         NKINDEX=0
         DO KJ=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ)==0)) CYCLE
         DO KA=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0)) CYCLE
            IF (ME==ATPROCS) THEN
               NKINDEX=NKINDEX+1
               KJ_INDEX(NKINDEX)=KJ
               KA_INDEX(NKINDEX)=KA
               PROCS_KA_KJ(KA,KJ)=ME
               INDEX_KA_KJ(KA,KJ)=NKINDEX
            ENDIF
            PROCS_KA_KJ(KA,KJ)=ATPROCS
            ATPROCS=ATPROCS+1
            IF (ATPROCS>PROCS-1) ATPROCS=0
         ENDDO
         ENDDO

         INDEX_KQ_KA=0
         PROCS_KQ_KA=-1
         ATPROCS=0
         NKINDEX=0
         DO KA=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0)) CYCLE
         DO KQ=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0.0_q)) CYCLE
            IF (ME==ATPROCS) THEN
               NKINDEX=NKINDEX+1
               KJ_INDEX(NKINDEX)=KA
               KA_INDEX(NKINDEX)=KQ
               PROCS_KQ_KA(KQ,KA)=ME
               INDEX_KQ_KA(KQ,KA)=NKINDEX
            ENDIF
            PROCS_KQ_KA(KQ,KA)=ATPROCS
            ATPROCS=ATPROCS+1
            IF (ATPROCS>PROCS-1) ATPROCS=0
         ENDDO
         ENDDO


         ELSE

         NKINDEX=0
         DO KJ=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ)==0)) CYCLE
         DO KA=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0)) CYCLE
               NKINDEX=NKINDEX+1
         ENDDO
         ENDDO

         ALLOCATE(KJ_INDEX(NKINDEX))
         ALLOCATE(KA_INDEX(NKINDEX))
         ALLOCATE(PROCS_KA_KJ(WDES%NKPTS,WDES%NKPTS))
         ALLOCATE(PROCS_KQ_KA(WDES%NKPTS,WDES%NKPTS))
         ALLOCATE(INDEX_KA_KJ(WDES%NKPTS,WDES%NKPTS))
         ALLOCATE(INDEX_KQ_KA(WDES%NKPTS,WDES%NKPTS))
         INDEX_KA_KJ=0
         PROCS_KA_KJ=-1
         NKINDEX=0
         DO KJ=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ)==0)) CYCLE
         DO KA=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0)) CYCLE
               NKINDEX=NKINDEX+1
               KJ_INDEX(NKINDEX)=KJ
               KA_INDEX(NKINDEX)=KA
               INDEX_KA_KJ(KA,KJ)=NKINDEX
         ENDDO
         ENDDO
         INDEX_KQ_KA=INDEX_KA_KJ
         PROCS_KQ_KA=PROCS_KA_KJ

         INDEX_KQ_KA=0
         PROCS_KQ_KA=-1
         NKINDEX=0
         DO KA=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0)) CYCLE
         DO KQ=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0.0_q)) CYCLE
               NKINDEX=NKINDEX+1
               KJ_INDEX(NKINDEX)=KA
               KA_INDEX(NKINDEX)=KQ
               INDEX_KQ_KA(KQ,KA)=NKINDEX
         ENDDO
         ENDDO


         ENDIF
     END SUBROUTINE SETUP_T2_KINDEX


      SUBROUTINE SETUP_BLOCKED_KINDEX(WDES,KB,KC,KBLOCK,MKB,MKC)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE mkpoints
         USE lattice
         USE constant
         USE base
         IMPLICIT NONE
         INTEGER COUNTER, RKB,RKC, KB, KC, KBLOCK
         TYPE(wavedes) WDES
         INTEGER MKB(WDES%NKPTS*WDES%NKPTS), MKC(WDES%NKPTS*WDES%NKPTS), ATPROCS

         MY_NKPTS=0
         ATPROCS=0
         DO RKB=KB,MIN(WDES%NKPTS,KB+KBLOCK-1)
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(RKB)==0)) CYCLE
         DO RKC=KC,MIN(WDES%NKPTS,KC+KBLOCK-1)
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(RKC)==0)) CYCLE
            IF (ME==ATPROCS) THEN
               MY_NKPTS=MY_NKPTS+1
               MKB(MY_NKPTS)=RKB
               MKC(MY_NKPTS)=RKC
               KP_TMP(MY_NKPTS*2)=RKB
               KP_TMP(MY_NKPTS*2+1)=RKC
               NTMP_KP(RKB)=MY_NKPTS*2
               NTMP_KP(RKC)=MY_NKPTS*2+1
            ENDIF
            ATPROCS=ATPROCS+1
            IF (ATPROCS>PROCS-1) ATPROCS=0
         ENDDO
         ENDDO

     END SUBROUTINE SETUP_BLOCKED_KINDEX

     SUBROUTINE SETUP_MAX_MY_NKPTS(WDES,KBLOCK)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE mkpoints
         USE lattice
         USE constant
         USE base
         IMPLICIT NONE
         INTEGER COUNTER, KBLOCK
         TYPE(wavedes) WDES
         INTEGER ATPROCS,RKB,RKC

         MAX_MY_NKPTS=0
         ATPROCS=0
         DO RKB=1,MIN(WDES%NKPTS,KBLOCK)
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(RKB)==0)) CYCLE
         DO RKC=1,MIN(WDES%NKPTS,KBLOCK)
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(RKC)==0)) CYCLE
            IF (ME==ATPROCS) THEN
               MAX_MY_NKPTS=MAX_MY_NKPTS+2
            ENDIF
            ATPROCS=ATPROCS+1
            IF (ATPROCS>PROCS-1) ATPROCS=0
         ENDDO
         ENDDO
         MAX_MY_NKPTS=MAX_MY_NKPTS+1 ! the one additional slot is needed to store T2_KA (the outermost k-point-loop)
         ALLOCATE(KP_TMP(MAX_MY_NKPTS))
         KP_TMP=0

     END SUBROUTINE SETUP_MAX_MY_NKPTS


      SUBROUTINE CALC_TRIPLES(WDES,WGW,W,ITERATION,KPOINTS,IO)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE mkpoints
         USE lattice
         USE constant
         USE base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavedes) WGW
         TYPE(wavespin) W
         TYPE (kpoints_struct) KPOINTS
         TYPE (in_struct) IO
         INTEGER ITERATION
         INTEGER KI,KI_,KJ,KK,KA,KA_,KB,KC,NI,NJ,NK,NA,NB,NC,KQ,KQ_,KMAX,KTMP,BKB,BKC,NG
         INTEGER RNBLOCKA, RNBLOCKB, RNBLOCKC, RNA, RNB, RNC
         integer :: time_array1(8),time_array2(8),emstriples
         INTEGER :: NTMP_T2,COMPK,KBLOCK,MKB(WDES%NKPTS*WDES%NKPTS),MKC(WDES%NKPTS*WDES%NKPTS)
         GDEF :: OVOV(VBMAX,NUNOCC,VBMAX,NUNOCC)
         REAL(q) :: OVOV_R(VBMAX,NUNOCC,VBMAX,NUNOCC), fac
         GDEF :: ECCSD, ETRIPLESofKA(REALNKPTS), ETRIPLESZofKA(REALNKPTS)
         INTEGER :: IERROR

         call date_and_time (values=time_array1)


         IF (IO%IU0>0) WRITE(IO%IU0,*) 
         IF (IO%IU0>0) WRITE(IO%IU0,*)'Calculating (T) contribution to correlation energy'
         ETRIPLES=(0.0_q,0.0_q)
         ETRIPLESofKA=(0.0_q,0.0_q)
         ETRIPLESZofKA=(0.0_q,0.0_q)
         ETRIPLESZ=(0.0_q,0.0_q)

         !an optimal NBLOCK should also increase with PROCS
         KBLOCK=MIN(64,WDES%NKPTS)

         CALL SETUP_MAX_MY_NKPTS(WDES,KBLOCK)

         !Here we need to allocate T2_TMP(1:MAX_MY_NKPTS) and CHI_AB_TMP(1:MAX_MY_NKPTS)
         IF (ALLOCATED(T2_TMP)) DEALLOCATE(T2_TMP)
         IF (ALLOCATED(T2_TMP_R)) DEALLOCATE(T2_TMP_R)
         IF (.not. LUSESHARED) THEN
            IF (LORBREAL) THEN
               ALLOCATE(T2_TMP_R(NUNOCC,NUNOCC,VBMAX,VBMAX,REALNKPTS,REALNKPTS,MAX_MY_NKPTS))
            ELSE
               ALLOCATE(T2_TMP(NUNOCC,NUNOCC,VBMAX,VBMAX,REALNKPTS,REALNKPTS,MAX_MY_NKPTS))
            ENDIF
         ENDIF
         IF (ALLOCATED(PW_AB_TMP)) DEALLOCATE(PW_AB_TMP)
         IF (.not. LUSESHARED) THEN
            ALLOCATE(PW_AB_TMP(NGVECTOR,NUNOCC,NUNOCC,REALNKPTS,MAX_MY_NKPTS))
         ENDIF
         IF (ALLOCATED(PW_AI_TMP)) DEALLOCATE(PW_AI_TMP)
         IF (.not. LUSESHARED) THEN
            ALLOCATE(PW_AI_TMP(NGVECTOR,VBMAX,NUNOCC,REALNKPTS,MAX_MY_NKPTS,2),stat=ierror)
         ENDIF
         IF (LORBREAL) THEN
            ALLOCATE(Wint_ijk_abc_R(NUNOCC,NUNOCC,NUNOCC,VBMAX,VBMAX,VBMAX))
            ALLOCATE(W_ijk_abc_R(NUNOCC,NUNOCC,NUNOCC,VBMAX,VBMAX,VBMAX))
            ALLOCATE(Wtmp_ijk_abc_R(NUNOCC,NUNOCC,NUNOCC,VBMAX,VBMAX,VBMAX))
         ELSE
            ALLOCATE(Wint_ijk_abc(NUNOCC,NUNOCC,NUNOCC,VBMAX,VBMAX,VBMAX),stat=ierror)
            ALLOCATE(W_ijk_abc(NUNOCC,NUNOCC,NUNOCC,VBMAX,VBMAX,VBMAX))
            ALLOCATE(Wtmp_ijk_abc(NUNOCC,NUNOCC,NUNOCC,VBMAX,VBMAX,VBMAX))
         ENDIF
         IF (LSFACTOR) ALLOCATE(Wint_ijk_abc_kq(NUNOCC,NUNOCC,NUNOCC,VBMAX,VBMAX,VBMAX,REALNKPTS),stat=ierror)
         IF (LSFACTOR) ALLOCATE(W_ijk_abc_kq(NUNOCC,NUNOCC,NUNOCC,VBMAX,VBMAX,VBMAX,REALNKPTS),stat=ierror)

         IF (ALLOCATED(NTMP_KP)) DEALLOCATE(NTMP_KP)
         ALLOCATE(NTMP_KP(WDES%NKPTS))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         IF (.not. LUSESHARED) THEN
            CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
              FTOD_PW_AI(1,1,1,1,1),FTOD_OC_AI(1,1,1,1,1),&
              (NUNOCC)*(VBMAX),&
              FTOD_PW_AI(1,1,1,1,2),FTOD_OC_AI(1,1,1,1,2), &
              VBMAX*NUNOCC,&
              OVOV(1,1,1,1),OVOV_R(1,1,1,1),zero)
         ENDIF

         ECCSD=zero
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
         DO NA=1,NUNOCC
         DO NB=1,NUNOCC
            IF (LORBREAL) THEN
               ECCSD=ECCSD+(OVOV_R(NI,NA,NJ,NB))*(2.0_q*T2_R(NA,NB,NI,NJ,1,1)-T2_R(NB,NA,NI,NJ,1,1))
            ELSE
               ECCSD=ECCSD+CONJG(OVOV(NI,NA,NJ,NB))*(2.0_q*T2(NA,NB,NI,NJ,1,1)-T2(NB,NA,NI,NJ,1,1))
            ENDIF
!            ECCSD=ECCSD+CONJG(OVOV(NI,NA,NJ,NB))*OVOV(NI,NA,NJ,NB)/ &
!                 REAL(W%CELTOT(NI,1,1)+W%CELTOT(NJ,1,1)-W%CELTOT(NA+VBMAX,1,1)-W%CELTOT(NB+VBMAX,1,1),kind=q)
         ENDDO
         ENDDO
         ENDDO
         ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


         DO KA_=1,KPOINTS_ORIG%NKPTS
         IF (IO%IU0>0) WRITE(IO%IU0,*)
         IF (IO%IU0>0) WRITE(IO%IU0,'("KA=",I4,3F10.4,", ")') KA_,KPOINTS_ORIG%VKPT(:,KA_)

!!!!!!!!!!!!!! for IBZ wedge find corresponding k-point in full BZ 
         KA=KPOINT_IN_FULL_GRID(KPOINTS_ORIG%VKPT(:,KA_),KPOINTS_FULL)
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0)) CYCLE

         !Here we need to (read in and) bcast T2_KA to all T2_tmp_1
         KP_TMP=0
         KP_TMP(1)=KA
         NTMP_KP(KA)=1
         IF (.not. LUSESHARED) THEN
            CALL BCAST2ALL_T2_KA(KA,WDES)
         !Here we need to (read in and) bcast CHI_AB_KA and CHI_AI_KA to all CHI_AB/AI_tmp_1
            CALL BCAST2ALL_FTOD_AB_KA(KA,WDES)
            CALL BCAST2ALL_FTOD_AI_KA(KA,WDES)
         ENDIF
         DO BKB=1,WDES%NKPTS,KBLOCK
         DO BKC=1,WDES%NKPTS,KBLOCK
            CALL SETUP_BLOCKED_KINDEX(WDES,BKB,BKC,KBLOCK,MKB,MKC)

            !Here we have to broadcast all T2_K and store the corresponding ones in T2_TMP(2:1+MY_NKPTS*2) for MKB and MKC
            IF (.not. LUSESHARED) THEN
               CALL BCAST2ALL_T2(WDES)
               CALL BCAST2ALL_FTOD_AB(WDES)
               CALL BCAST2ALL_FTOD_AI(WDES)
            ENDIF
            !We also need to broadcast CHI_AB/AI_(1:MYNKPTS*2), CHI_AB/AI_tmp_(2:1+MY_NKPTS*2) for MKB and MKC

            DO COMPK=1,MY_NKPTS
               KB=MKB(COMPK)
               KC=MKC(COMPK)
               CALL GWPROGRESS(IO%IU0, KB,WDES%NKPTS,KC,WDES%NKPTS)
               DO KI=1,WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI)==0.0_q)) CYCLE
               DO KJ=1,WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ)==0.0_q)) CYCLE
               IF ((KA_==1) .and. (BKB==1) .and. (BKC==1) .and. (COMPK==1)) CALL DATE_AND_TIME(values=time_array1)
               !calculate KK
               KK=KPOINT_IN_FULL_GRID(-WDES%VKPT(:,(KI))- &
                  WDES%VKPT(:,KJ)+WDES%VKPT(:,KA)+WDES%VKPT(:,KB)+WDES%VKPT(:,KC),KPOINTS_FULL)

               ! call routine that calculates and adds
               ! W_ijk_abc KA, MKB(COMPK), MKC(COMPK)
               ! W_kji_cba MKC(COMPK), MKB(COMPK), KA
               ! W_ikj_acb KA, MKC(COMPK), MKB(COMPK)
               ! .........

               ! evaluate Wtmp_ijk_abc=4/3*W_ijk_abc-2*W_ijk_acb+2/3*W_ijk_bca
               IF (LORBREAL) THEN

                  fac=1.0_q*KPOINTS_ORIG%WTKPT(KA_)
                  Wtmp_ijk_abc_R=(0.0_qs)
                  CALL CALC_Wint_ijk_abc(KI,KJ,KK,KA,KC,KB,WDES,KPOINTS)
                  DO NB=1,NUNOCC
                  DO NC=1,NUNOCC
                     Wtmp_ijk_abc_R(:,NB,NC,:,:,:)=Wtmp_ijk_abc_R(:,NB,NC,:,:,:)-2.0_qs*Wint_ijk_abc_R(:,NC,NB,:,:,:)
                  ENDDO
                  ENDDO
                  CALL CALC_Wint_ijk_abc(KI,KJ,KK,KB,KC,KA,WDES,KPOINTS)
                  DO NB=1,NUNOCC
                  DO NC=1,NUNOCC
                    Wtmp_ijk_abc_R(:,NB,NC,:,:,:)=Wtmp_ijk_abc_R(:,NB,NC,:,:,:)+2.0_qs/3.0_qs*Wint_ijk_abc_R(NB,NC,:,:,:,:)
                  ENDDO
                  ENDDO
                  CALL CALC_Wint_ijk_abc(KI,KJ,KK,KA,KB,KC,WDES,KPOINTS)
                  Wtmp_ijk_abc_R(:,:,:,:,:,:)=Wtmp_ijk_abc_R(:,:,:,:,:,:)+4.0_qs/3.0_qs*Wint_ijk_abc_R(:,:,:,:,:,:)

                  ! contract Wbar_ijk_abc W_ijk_abc
                  DO NK=1,VBMAX
                  DO NJ=1,VBMAX
                  DO NI=1,VBMAX
                  DO NC=1,NUNOCC
                  DO NB=1,NUNOCC
                  DO NA=1,NUNOCC
                     ETRIPLES=ETRIPLES+ &
                           fac*(Wtmp_ijk_abc_R(NA,NB,NC,NI,NJ,NK))*(Wint_ijk_abc_R(NA,NB,NC,NI,NJ,NK))/ &
                              REAL(W%CELTOT(NI,KI,1)+W%CELTOT(NJ,KJ,1)+W%CELTOT(NK,(KK),1)-W%CELTOT(NA+VBMAX,KA,1)-W%CELTOT(NB+VBMAX,KB,1)-W%CELTOT(NC+VBMAX,KC,1),kind=q)
                     ETRIPLESofKA(RKIofKI(KA))=ETRIPLESofKA(RKIofKI(KA))+ &
                           fac*(Wtmp_ijk_abc_R(NA,NB,NC,NI,NJ,NK))*(Wint_ijk_abc_R(NA,NB,NC,NI,NJ,NK))/ &
                              REAL(W%CELTOT(NI,KI,1)+W%CELTOT(NJ,KJ,1)+W%CELTOT(NK,(KK),1)-W%CELTOT(NA+VBMAX,KA,1)-W%CELTOT(NB+VBMAX,KB,1)-W%CELTOT(NC+VBMAX,KC,1),kind=q)
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                 
               ! evaluate Vbar_ijk_abc=4/3*V_ijk_abc-2*V_ijk_acb+2/3*V_ijk_bca and store in Wtmp_ijk_abc

                  CALL CALC_Ztmp_ijk_abc(KI,KJ,KK,KA,KB,KC,WDES,KPOINTS)
  
                  ! contract Zbar_ijk_abc W_ijk_abc
                  DO NK=1,VBMAX
                  DO NJ=1,VBMAX
                  DO NI=1,VBMAX
                  DO NC=1,NUNOCC
                  DO NB=1,NUNOCC
                  DO NA=1,NUNOCC
                     ETRIPLESZ=ETRIPLESZ+ &
                           fac*(Wtmp_ijk_abc_R(NA,NB,NC,NI,NJ,NK))*(Wint_ijk_abc_R(NA,NB,NC,NI,NJ,NK))/ &
                              REAL(W%CELTOT(NI,KI,1)+W%CELTOT(NJ,KJ,1)+W%CELTOT(NK,(KK),1)-W%CELTOT(NA+VBMAX,KA,1)-W%CELTOT(NB+VBMAX,KB,1)-W%CELTOT(NC+VBMAX,KC,1),kind=q)
                     ETRIPLESZofKA(RKIofKI(KA))=ETRIPLESZofKA(RKIofKI(KA))+ &
                           fac*(Wtmp_ijk_abc_R(NA,NB,NC,NI,NJ,NK))*(Wint_ijk_abc_R(NA,NB,NC,NI,NJ,NK))/ &
                              REAL(W%CELTOT(NI,KI,1)+W%CELTOT(NJ,KJ,1)+W%CELTOT(NK,(KK),1)-W%CELTOT(NA+VBMAX,KA,1)-W%CELTOT(NB+VBMAX,KB,1)-W%CELTOT(NC+VBMAX,KC,1),kind=q)

                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
 

               ELSE

                  fac=1.0_q*KPOINTS_ORIG%WTKPT(KA_)
                  Wtmp_ijk_abc=(0.0_qs,0.0_qs)
                  CALL CALC_Wint_ijk_abc(KI,KJ,KK,KA,KC,KB,WDES,KPOINTS)
                  DO NB=1,NUNOCC
                  DO NC=1,NUNOCC
                     Wtmp_ijk_abc(:,NB,NC,:,:,:)=Wtmp_ijk_abc(:,NB,NC,:,:,:)-2.0_qs*Wint_ijk_abc(:,NC,NB,:,:,:)
                  ENDDO
                  ENDDO
                  CALL CALC_Wint_ijk_abc(KI,KJ,KK,KB,KC,KA,WDES,KPOINTS)
                  DO NB=1,NUNOCC
                  DO NC=1,NUNOCC
                     Wtmp_ijk_abc(:,NB,NC,:,:,:)=Wtmp_ijk_abc(:,NB,NC,:,:,:)+2.0_qs/3.0_qs*Wint_ijk_abc(NB,NC,:,:,:,:)
                  ENDDO
                  ENDDO
                  CALL CALC_Wint_ijk_abc(KI,KJ,KK,KA,KB,KC,WDES,KPOINTS)
                  Wtmp_ijk_abc(:,:,:,:,:,:)=Wtmp_ijk_abc(:,:,:,:,:,:)+4.0_qs/3.0_qs*Wint_ijk_abc(:,:,:,:,:,:)

                  IF (LSFACTOR) THEN

                  DO NG=1,NGVECTOR

                  CALL CALC_Wint_ijk_abc_kq(KI,KJ,KK,KA,KB,KC,NG,WDES,KPOINTS)

                  DO KQ=1,WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0.0_q)) CYCLE
               ! contract Wbar_ijk_abc W_ijk_abc
                  DO NK=1,VBMAX
                  DO NJ=1,VBMAX
                  DO NI=1,VBMAX
                  DO NC=1,NUNOCC
                  DO NB=1,NUNOCC
                  DO NA=1,NUNOCC
                     SF_ETRIPLES(NG,RKQofKQ(KQ))=SF_ETRIPLES(NG,RKQofKQ(KQ))+ &
                           fac*(Wtmp_ijk_abc(NA,NB,NC,NI,NJ,NK))*GCONJG(Wint_ijk_abc_kq(NA,NB,NC,NI,NJ,NK,RKQofKQ(KQ)))/ &
                              REAL(W%CELTOT(NI,KI,1)+W%CELTOT(NJ,KJ,1)+W%CELTOT(NK,(KK),1)-W%CELTOT(NA+VBMAX,KA,1)-W%CELTOT(NB+VBMAX,KB,1)-W%CELTOT(NC+VBMAX,KC,1),kind=q)
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
 
                  ENDDO ! KQ

                  ENDDO ! NG

                  ENDIF
                         

               ! contract Wbar_ijk_abc W_ijk_abc
                  DO NK=1,VBMAX
                  DO NJ=1,VBMAX
                  DO NI=1,VBMAX
                  DO NC=1,NUNOCC
                  DO NB=1,NUNOCC
                  DO NA=1,NUNOCC
                     ETRIPLES=ETRIPLES+ &
                           fac*(Wtmp_ijk_abc(NA,NB,NC,NI,NJ,NK))*GCONJG(Wint_ijk_abc(NA,NB,NC,NI,NJ,NK))/ &
                              REAL(W%CELTOT(NI,KI,1)+W%CELTOT(NJ,KJ,1)+W%CELTOT(NK,(KK),1)-W%CELTOT(NA+VBMAX,KA,1)-W%CELTOT(NB+VBMAX,KB,1)-W%CELTOT(NC+VBMAX,KC,1),kind=q)
                      ETRIPLESofKA(RKIofKI(KA))=ETRIPLESofKA(RKIofKI(KA))+ &
                           fac*(Wtmp_ijk_abc(NA,NB,NC,NI,NJ,NK))*GCONJG(Wint_ijk_abc(NA,NB,NC,NI,NJ,NK))/ &
                              REAL(W%CELTOT(NI,KI,1)+W%CELTOT(NJ,KJ,1)+W%CELTOT(NK,(KK),1)-W%CELTOT(NA+VBMAX,KA,1)-W%CELTOT(NB+VBMAX,KB,1)-W%CELTOT(NC+VBMAX,KC,1),kind=q)
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
              
               ! evaluate Vbar_ijk_abc=4/3*V_ijk_abc-2*V_ijk_acb+2/3*V_ijk_bca and store in Wtmp_ijk_abc

                  CALL CALC_Ztmp_ijk_abc(KI,KJ,KK,KA,KB,KC,WDES,KPOINTS)

                  IF (LSFACTOR) THEN

                  DO NG=1,NGVECTOR
                  CALL CALC_Wint_ijk_abc_kq(KI,KJ,KK,KA,KB,KC,NG,WDES,KPOINTS)

                  DO KQ=1,WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0.0_q)) CYCLE
               ! contract Wbar_ijk_abc W_ijk_abc
                  DO NK=1,VBMAX
                  DO NJ=1,VBMAX
                  DO NI=1,VBMAX
                  DO NC=1,NUNOCC
                  DO NB=1,NUNOCC
                  DO NA=1,NUNOCC
                     SF_ETRIPLESZ(NG,RKQofKQ(KQ))=SF_ETRIPLESZ(NG,RKQofKQ(KQ))+ &
                           fac*(Wtmp_ijk_abc(NA,NB,NC,NI,NJ,NK))*(Wint_ijk_abc_kq(NA,NB,NC,NI,NJ,NK,RKQofKQ(KQ)))/ &
                              REAL(W%CELTOT(NI,KI,1)+W%CELTOT(NJ,KJ,1)+W%CELTOT(NK,(KK),1)-W%CELTOT(NA+VBMAX,KA,1)-W%CELTOT(NB+VBMAX,KB,1)-W%CELTOT(NC+VBMAX,KC,1),kind=q)

                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
 
                  ENDDO ! KQ

                  ENDDO ! NG

                  ENDIF
 

               ! contract Zbar_ijk_abc W_ijk_abc
                  DO NK=1,VBMAX
                  DO NJ=1,VBMAX
                  DO NI=1,VBMAX
                  DO NC=1,NUNOCC
                  DO NB=1,NUNOCC
                  DO NA=1,NUNOCC
                     ETRIPLESZ=ETRIPLESZ+ &
                           fac*(Wtmp_ijk_abc(NA,NB,NC,NI,NJ,NK))*(Wint_ijk_abc(NA,NB,NC,NI,NJ,NK))/ &
                              REAL(W%CELTOT(NI,KI,1)+W%CELTOT(NJ,KJ,1)+W%CELTOT(NK,(KK),1)-W%CELTOT(NA+VBMAX,KA,1)-W%CELTOT(NB+VBMAX,KB,1)-W%CELTOT(NC+VBMAX,KC,1),kind=q)
                     ETRIPLESZofKA(KA)=ETRIPLESZofKA(KA)+ &
                           fac*(Wtmp_ijk_abc(NA,NB,NC,NI,NJ,NK))*(Wint_ijk_abc(NA,NB,NC,NI,NJ,NK))/ &
                              REAL(W%CELTOT(NI,KI,1)+W%CELTOT(NJ,KJ,1)+W%CELTOT(NK,(KK),1)-W%CELTOT(NA+VBMAX,KA,1)-W%CELTOT(NB+VBMAX,KB,1)-W%CELTOT(NC+VBMAX,KC,1),kind=q)

                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO


               ENDIF

                 IF ((KI==1) .and. (KJ==1) .and. (KA_==1) .and. (BKB==1) .and. (BKC==1) .and. (COMPK==1)) THEN
                    CALL DATE_AND_TIME(values=time_array2)
                    EMSTRIPLES=((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
                      time_array2(8)-time_array1(8))
                  IF (IO%IU0>0) WRITE(IO%IU6,*)'calculating (T) will take approximately ',emstriples*KPOINTS_ORIG%NKPTS*REALNKPTS*REALNKPTS/KBLOCK*REALNKPTS/KBLOCK*REALNKPTS*MY_NKPTS/1000.0_q,' seconds'
                  IF (IO%IU0>0) WRITE(*,*)'calculating (T) will take approximately ',emstriples*KPOINTS_ORIG%NKPTS*REALNKPTS*REALNKPTS/KBLOCK*REALNKPTS/KBLOCK*REALNKPTS*MY_NKPTS/1000.0_q,' seconds'
                  IF (IO%IU0>0) CALL flush(IO%IU6)
                 ENDIF

               ENDDO
               ENDDO

            ENDDO

         ENDDO
         ENDDO
         
         CALLMPI( M_sum_d(WGW%COMM_INTER, ETRIPLESofKA(RKIofKI(KA)), 2)) 
         CALLMPI( M_sum_d(WGW%COMM_INTER, ETRIPLESZofKA(RKIofKI(KA)), 2)) 
         IF (IO%IU0>=0) WRITE(IO%IU6,*) '[T] correlation energy contribution of k-point '
         IF (IO%IU0>=0) WRITE(IO%IU6,*) WDES%VKPT(:,KA),' is'
         IF (IO%IU0>=0) WRITE(IO%IU6,*) ETRIPLESofKA(RKIofKI(KA))*REAL(1.0_q/REALNKPTS/REALNKPTS/REALNKPTS/REALNKPTS,kind=q),' eV'
         IF (IO%IU0>=0) WRITE(IO%IU6,*) '(T) correlation energy contribution of k-point '
         IF (IO%IU0>=0) WRITE(IO%IU6,*) WDES%VKPT(:,KA),' is '
         IF (IO%IU0>=0) WRITE(IO%IU6,*) ETRIPLESZofKA(RKIofKI(KA))*REAL(1.0_q/REALNKPTS/REALNKPTS/REALNKPTS/REALNKPTS,kind=q),'eV'

         ENDDO

         CALLMPI( M_sum_d(WGW%COMM_INTER, ETRIPLES, 2))
         CALLMPI( M_sum_d(WGW%COMM_INTER, ETRIPLESZ, 2))
         IF (IO%IU0>0) WRITE(IO%IU0,*)
         IF (IO%IU0>0) WRITE(*,*)'Total [T] contribution to the correlation energy is',ETRIPLES*REAL(1.0_q/REALNKPTS/REALNKPTS/REALNKPTS/REALNKPTS,kind=q) !,T2
         IF (IO%IU0>0) WRITE(*,*)'Total (T) contribution to the correlation energy is',(ETRIPLES+ETRIPLESZ)*REAL(1.0_q/REALNKPTS/REALNKPTS/REALNKPTS/REALNKPTS,kind=q) !,T2
 
         IF (LSFACTOR) THEN

         WRITE(*,*)'Writing (T) contribution to S-Factor'
         CALLMPI( M_sum_d(WGW%COMM_INTER, SF_ETRIPLES(1,1), NGVECTOR*REALNKPTS*2)) 
         CALLMPI( M_sum_d(WGW%COMM_INTER, SF_ETRIPLESZ(1,1), NGVECTOR*REALNKPTS*2)) 
         SF_ETRIPLES=SF_ETRIPLES*REAL(1.0_q/REALNKPTS/REALNKPTS/REALNKPTS/REALNKPTS,kind=q)
         SF_ETRIPLESZ=SF_ETRIPLESZ*REAL(1.0_q/REALNKPTS/REALNKPTS/REALNKPTS/REALNKPTS,kind=q)

         IF (IO%IU0>=0) THEN
            WRITE(*,*)'Writing to pTCORRofG ...'
            OPEN(unit = 7,file = "pTCORRofG")
            WRITE(7,'(5G20.12)')  0.0, 0.0, 0.0, 0.0, 0.0
            DO KQ=1,WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0.0_q)) CYCLE
            DO NG=1,NGVECTOR
               IF (GVECLEN(NG,RKQofKQ(KQ))==0.0_q) CYCLE
               WRITE(7,'(5G20.12)') GVECX(NG,RKQofKQ(KQ)),GVECY(NG,RKQofKQ(KQ)),GVECZ(NG,RKQofKQ(KQ)),GVECLEN(NG,RKQofKQ(KQ)),REAL(SF_ETRIPLES(NG,RKQofKQ(KQ))+SF_ETRIPLESZ(NG,RKQofKQ(KQ)),kind=q)
            ENDDO
            ENDDO
            CLOSE(7)
         ENDIF


         ENDIF

      END Subroutine CALC_TRIPLES

      SUBROUTINE CONTR_FTOD(NGV,NHV,PW_XY,OC_XY,LNXY,PW_XY_,OC_XY_, &
                            LNXY_,XY2E4ORB,XY2E4ORB_R,beta)
      IMPLICIT NONE
      INTEGER :: NGV,NHV,LNXY,LNXY_
      GDEF :: XY2E4ORB,PW_XY,OC_XY, &
              PW_XY_,OC_XY_, beta
      REAL(q) :: XY2E4ORB_R
              
      IF (.not. LORBREAL) THEN
            CALL ZGEMM(trans,'n',(LNXY),(LNXY_),&
                 (NGV),(1._q,0._q),PW_XY,(NGV),&
                 PW_XY_,(NGV),&
                 beta, XY2E4ORB,(LNXY))

            IF (ASSOCIATED(H)) THEN

            CALL ZGEMM(trans,'n',(LNXY),(LNXY_),&
                 (NHV),(1._q,0._q),OC_XY,(NHV),&
                 OC_XY_,(NHV),&
                 (1._q,0._q), XY2E4ORB,(LNXY))
            ENDIF

      ELSE
     
            CALL DGEMM(trans,'n',(LNXY),(LNXY_),&
                 (NGV*2),(1._q),PW_XY,(NGV*2),&
                 PW_XY_,(NGV*2),&
                 REAL(beta,kind=q), XY2E4ORB_R,(LNXY))

            IF (ASSOCIATED(H)) THEN
            CALL DGEMM(trans,'n',(LNXY),(LNXY_),&
                 (NHV*2),(1._q),OC_XY,(NHV*2),&
                 OC_XY_,(NHV*2),&
                 (1._q), XY2E4ORB_R,(LNXY))
            ENDIF

      ENDIF
      
      END SUBROUTINE

       SUBROUTINE CALC_Ztmp_ijk_abc(KI,KJ,KK,KA,KB,KC,WDES,KPOINTS)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE mkpoints
         USE lattice
         USE constant
         USE base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE (kpoints_struct) KPOINTS
         INTEGER ITERATION,NBLOCK,MKIMAX,MKIMIN,NK,NC
         INTEGER KI,KJ,KK,KA,KB,KC,NI,NJ,NA,NB
         COMPLEX(q) :: VVOO(NUNOCC,NUNOCC,VBMAX,VBMAX)

         IF (LORBREAL) THEN
            Wtmp_ijk_abc_R=(0.0_qs)
         ELSE
            Wtmp_ijk_abc=(0.0_qs,0.0_qs)
         ENDIF

         CALL CALC_Z_ijk_abc(KI,KJ,KK,KA,KB,KC,WDES,KPOINTS)
         IF (LORBREAL) THEN
            Wtmp_ijk_abc_R(:,:,:,:,:,:)=Wtmp_ijk_abc_R(:,:,:,:,:,:)+4.0_qs/3.0_qs*W_ijk_abc_R(:,:,:,:,:,:)
         ELSE
            Wtmp_ijk_abc(:,:,:,:,:,:)=Wtmp_ijk_abc(:,:,:,:,:,:)+4.0_qs/3.0_qs*W_ijk_abc(:,:,:,:,:,:)
         ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         CALL CALC_Z_ijk_abc(KI,KJ,KK,KB,KC,KA,WDES,KPOINTS)
         IF (LORBREAL) THEN

         DO NI=1,VBMAX
         DO NJ=1,VBMAX
         DO NK=1,VBMAX
            DO NC=1,NUNOCC
            DO NB=1,NUNOCC
            DO NA=1,NUNOCC
             Wtmp_ijk_abc_R(NA,NB,NC,NI,NJ,NK)=Wtmp_ijk_abc_R(NA,NB,NC,NI,NJ,NK)+2.0_qs/3.0_qs*W_ijk_abc_R(NB,NC,NA,NI,NJ,NK)
            ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDDO

         ELSE

         DO NI=1,VBMAX
         DO NJ=1,VBMAX
         DO NK=1,VBMAX
            DO NC=1,NUNOCC
            DO NB=1,NUNOCC
            DO NA=1,NUNOCC
             Wtmp_ijk_abc(NA,NB,NC,NI,NJ,NK)=Wtmp_ijk_abc(NA,NB,NC,NI,NJ,NK)+2.0_qs/3.0_qs*W_ijk_abc(NB,NC,NA,NI,NJ,NK)
            ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDDO

         ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         CALL CALC_Z_ijk_abc(KI,KJ,KK,KA,KC,KB,WDES,KPOINTS)
         IF (LORBREAL) THEN

         DO NI=1,VBMAX
         DO NJ=1,VBMAX
         DO NK=1,VBMAX
            DO NC=1,NUNOCC
            DO NB=1,NUNOCC
            DO NA=1,NUNOCC
             Wtmp_ijk_abc_R(NA,NB,NC,NI,NJ,NK)=Wtmp_ijk_abc_R(NA,NB,NC,NI,NJ,NK)-2.0_qs*W_ijk_abc_R(NA,NC,NB,NI,NJ,NK)
            ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDDO


         ELSE

         DO NI=1,VBMAX
         DO NJ=1,VBMAX
         DO NK=1,VBMAX
            DO NC=1,NUNOCC
            DO NB=1,NUNOCC
            DO NA=1,NUNOCC
             Wtmp_ijk_abc(NA,NB,NC,NI,NJ,NK)=Wtmp_ijk_abc(NA,NB,NC,NI,NJ,NK)-2.0_qs*W_ijk_abc(NA,NC,NB,NI,NJ,NK)
            ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDDO

         ENDIF

      END SUBROUTINE CALC_Ztmp_ijk_abc


      SUBROUTINE CALC_Wint_ijk_abc(KI,KJ,KK,KA,KB,KC,WDES,KPOINTS)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE mkpoints
         USE lattice
         USE constant
         USE base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE (kpoints_struct) KPOINTS
         INTEGER ITERATION,NBLOCK,MKIMAX,MKIMIN,NK,NC
         INTEGER KI,KJ,KK,KA,KB,KC,NI,NJ,NA,NB
         COMPLEX(q) :: VVOO(NUNOCC,NUNOCC,VBMAX,VBMAX)

         IF (LORBREAL) THEN
          Wint_ijk_abc_R(:,:,:,:,:,:)=(0.0_qs)
         ELSE
          Wint_ijk_abc(:,:,:,:,:,:)=(0.0_qs,0.0_qs)
         ENDIF
! ijk_abc=ijk_abc+ jik_bac
!         WRITE(*,*)'before calc_w_ijk'
         CALL CALC_W_ijk_abc(KJ,KI,KK,KB,KA,KC,WDES,KPOINTS)
!         WRITE(*,*)'after calc_w_ijk'
         IF (LORBREAL) THEN
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            DO NA=1,NUNOCC
            DO NB=1,NUNOCC
             Wint_ijk_abc_R(NA,NB,:,NI,NJ,:)=Wint_ijk_abc_R(NA,NB,:,NI,NJ,:)+W_ijk_abc_R(NB,NA,:,NJ,NI,:)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ELSE
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            DO NA=1,NUNOCC
            DO NB=1,NUNOCC
             Wint_ijk_abc(NA,NB,:,NI,NJ,:)=Wint_ijk_abc(NA,NB,:,NI,NJ,:)+W_ijk_abc(NB,NA,:,NJ,NI,:)
            ENDDO
            ENDDO
         ENDDO
         ENDDO

         ENDIF

! ijk_abc=ijk_abc+ kji_cba
         CALL CALC_W_ijk_abc(KK,KJ,KI,KC,KB,KA,WDES,KPOINTS)
         IF (LORBREAL) THEN
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            DO NA=1,NUNOCC
            DO NB=1,NUNOCC
             Wint_ijk_abc_R(NA,NB,:,NI,NJ,:)=Wint_ijk_abc_R(NA,NB,:,NI,NJ,:)+W_ijk_abc_R(:,NB,NA,:,NJ,NI)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ELSE
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            DO NA=1,NUNOCC
            DO NB=1,NUNOCC
             Wint_ijk_abc(NA,NB,:,NI,NJ,:)=Wint_ijk_abc(NA,NB,:,NI,NJ,:)+W_ijk_abc(:,NB,NA,:,NJ,NI)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF

! ijk_abc=ijk_abc+ ikj_acb
         CALL CALC_W_ijk_abc(KI,KK,KJ,KA,KC,KB,WDES,KPOINTS)
         IF (LORBREAL) THEN
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            DO NA=1,NUNOCC
            DO NB=1,NUNOCC
             Wint_ijk_abc_R(NA,NB,:,NI,NJ,:)=Wint_ijk_abc_R(NA,NB,:,NI,NJ,:)+W_ijk_abc_R(NA,:,NB,NI,:,NJ)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ELSE
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            DO NA=1,NUNOCC
            DO NB=1,NUNOCC
             Wint_ijk_abc(NA,NB,:,NI,NJ,:)=Wint_ijk_abc(NA,NB,:,NI,NJ,:)+W_ijk_abc(NA,:,NB,NI,:,NJ)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF

! ijk_abc=ijk_abc+ kij_cab
         CALL CALC_W_ijk_abc(KK,KI,KJ,KC,KA,KB,WDES,KPOINTS)
         IF (LORBREAL) THEN
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            DO NA=1,NUNOCC
            DO NB=1,NUNOCC
             Wint_ijk_abc_R(NA,NB,:,NI,NJ,:)=Wint_ijk_abc_R(NA,NB,:,NI,NJ,:)+W_ijk_abc_R(:,NA,NB,:,NI,NJ)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ELSE
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            DO NA=1,NUNOCC
            DO NB=1,NUNOCC
             Wint_ijk_abc(NA,NB,:,NI,NJ,:)=Wint_ijk_abc(NA,NB,:,NI,NJ,:)+W_ijk_abc(:,NA,NB,:,NI,NJ)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF

! ijk_abc=ijk_abc+ jki_bca
         CALL CALC_W_ijk_abc(KJ,KK,KI,KB,KC,KA,WDES,KPOINTS)
         IF (LORBREAL) THEN
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            DO NA=1,NUNOCC
            DO NB=1,NUNOCC
             Wint_ijk_abc_R(NA,NB,:,NI,NJ,:)=Wint_ijk_abc_R(NA,NB,:,NI,NJ,:)+W_ijk_abc_R(NB,:,NA,NJ,:,NI)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ELSE
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            DO NA=1,NUNOCC
            DO NB=1,NUNOCC
             Wint_ijk_abc(NA,NB,:,NI,NJ,:)=Wint_ijk_abc(NA,NB,:,NI,NJ,:)+W_ijk_abc(NB,:,NA,NJ,:,NI)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF

! ijk_abc=ijk_abc+ijk_abc
         CALL CALC_W_ijk_abc(KI,KJ,KK,KA,KB,KC,WDES,KPOINTS)
         IF (LORBREAL) THEN
         Wint_ijk_abc_R(:,:,:,:,:,:)=Wint_ijk_abc_R(:,:,:,:,:,:)+W_ijk_abc_R(:,:,:,:,:,:)
         ELSE
         Wint_ijk_abc(:,:,:,:,:,:)=Wint_ijk_abc(:,:,:,:,:,:)+W_ijk_abc(:,:,:,:,:,:)
         ENDIF


      END SUBROUTINE CALC_Wint_ijk_abc


      SUBROUTINE CALC_Wint_ijk_abc_kq(KI,KJ,KK,KA,KB,KC,NG,WDES,KPOINTS)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE mkpoints
         USE lattice
         USE constant
         USE base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE (kpoints_struct) KPOINTS
         INTEGER ITERATION,NBLOCK,MKIMAX,MKIMIN,NK,NC
         INTEGER KI,KJ,KK,KA,KB,KC,NI,NJ,NA,NB,NG
         COMPLEX(q) :: VVOO(NUNOCC,NUNOCC,VBMAX,VBMAX)

         IF (LORBREAL) THEN
         ELSE
          Wint_ijk_abc_kq(:,:,:,:,:,:,:)=(0.0_qs,0.0_qs)
         ENDIF
! ijk_abc=ijk_abc+ jik_bac
         CALL CALC_W_ijk_abc_kq(KJ,KI,KK,KB,KA,KC,NG,WDES,KPOINTS)
         IF (LORBREAL) THEN
         ELSE
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            DO NA=1,NUNOCC
            DO NB=1,NUNOCC
             Wint_ijk_abc_kq(NA,NB,:,NI,NJ,:,:)=Wint_ijk_abc_kq(NA,NB,:,NI,NJ,:,:)+W_ijk_abc_kq(NB,NA,:,NJ,NI,:,:)
            ENDDO
            ENDDO
         ENDDO
         ENDDO

         ENDIF

! ijk_abc=ijk_abc+ kji_cba
         CALL CALC_W_ijk_abc_kq(KK,KJ,KI,KC,KB,KA,NG,WDES,KPOINTS)
         IF (LORBREAL) THEN
         ELSE
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            DO NA=1,NUNOCC
            DO NB=1,NUNOCC
             Wint_ijk_abc_kq(NA,NB,:,NI,NJ,:,:)=Wint_ijk_abc_kq(NA,NB,:,NI,NJ,:,:)+W_ijk_abc_kq(:,NB,NA,:,NJ,NI,:)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF

! ijk_abc=ijk_abc+ ikj_acb
         CALL CALC_W_ijk_abc_kq(KI,KK,KJ,KA,KC,KB,NG,WDES,KPOINTS)
         IF (LORBREAL) THEN
         ELSE
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            DO NA=1,NUNOCC
            DO NB=1,NUNOCC
             Wint_ijk_abc_kq(NA,NB,:,NI,NJ,:,:)=Wint_ijk_abc_kq(NA,NB,:,NI,NJ,:,:)+W_ijk_abc_kq(NA,:,NB,NI,:,NJ,:)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF

! ijk_abc=ijk_abc+ kij_cab
         CALL CALC_W_ijk_abc_kq(KK,KI,KJ,KC,KA,KB,NG,WDES,KPOINTS)
         IF (LORBREAL) THEN
         ELSE
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            DO NA=1,NUNOCC
            DO NB=1,NUNOCC
             Wint_ijk_abc_kq(NA,NB,:,NI,NJ,:,:)=Wint_ijk_abc_kq(NA,NB,:,NI,NJ,:,:)+W_ijk_abc_kq(:,NA,NB,:,NI,NJ,:)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF

! ijk_abc=ijk_abc+ jki_bca
         CALL CALC_W_ijk_abc_kq(KJ,KK,KI,KB,KC,KA,NG,WDES,KPOINTS)
         IF (LORBREAL) THEN
         ELSE
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            DO NA=1,NUNOCC
            DO NB=1,NUNOCC
             Wint_ijk_abc_kq(NA,NB,:,NI,NJ,:,:)=Wint_ijk_abc_kq(NA,NB,:,NI,NJ,:,:)+W_ijk_abc_kq(NB,:,NA,NJ,:,NI,:)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF

! ijk_abc=ijk_abc+ijk_abc
         CALL CALC_W_ijk_abc_kq(KI,KJ,KK,KA,KB,KC,NG,WDES,KPOINTS)
         IF (LORBREAL) THEN
         ELSE
         Wint_ijk_abc_kq(:,:,:,:,:,:,:)=Wint_ijk_abc_kq(:,:,:,:,:,:,:)+W_ijk_abc_kq(:,:,:,:,:,:,:)
         ENDIF


      END SUBROUTINE CALC_Wint_ijk_abc_kq



      SUBROUTINE CALC_Z_ijk_abc(KI,KJ,KK,KA,KB,KC,WDES,KPOINTS)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE mkpoints
         USE lattice
         USE constant
         USE base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE (kpoints_struct) KPOINTS
         INTEGER ITERATION,NBLOCK,MKIMAX,MKIMIN
         INTEGER KI,KJ,KK,KA,KB,KC,NI,NJ,NK,NM,NA,NB,NC,KQ,KQ_,KF,KM,KMAX,KE,KE2,KM2
         INTEGER ierr
         COMPLEX(q) :: OVOV(VBMAX,NUNOCC,VBMAX,NUNOCC)
         REAL(q) :: OVOV_R(VBMAX,NUNOCC,VBMAX,NUNOCC)
         COMPLEX(q) KW
         INTEGER :: time_array1(8), time_array2(8)

         call date_and_time (values=time_array1)
         KW=(1.0_qs,0.0_qs)*WTKPT
         IF (LORBREAL) THEN
         W_ijk_abc_R=(0.0_qs)
         ELSE
         W_ijk_abc=(0.0_qs,0.0_qs)
         ENDIF

         KQ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,(KB))- &
               WDES%VKPT(:,KJ),KPOINTS_FULL)
         IF (KI==KA) THEN

         IF (.not. LUSESHARED) THEN
            CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
              PW_AI_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KB),1),OC_AI_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KB),1),&
              (VBMAX)*(NUNOCC),&
              PW_AI_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KC),2),OC_AI_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KC),2), &
              VBMAX*NUNOCC,&
              OVOV(1,1,1,1),OVOV_R(1,1,1,1),zero)
            !jbkc
         ELSE
            CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
              FTOD_PW_AI(1,1,1,INDEX_KQ_KA(KQ,KB),1),FTOD_OC_AI(1,1,1,INDEX_KQ_KA(KQ,KB),1),&
              (VBMAX)*(NUNOCC),&
              FTOD_PW_AI(1,1,1,INDEX_KQ_KA(KQ,KC),2),FTOD_OC_AI(1,1,1,INDEX_KQ_KA(KQ,KC),2), &
              VBMAX*NUNOCC,&
              OVOV(1,1,1,1),OVOV_R(1,1,1,1),zero)
            !jbkc
         ENDIF

            IF (LORBREAL) THEN
               DO NI=1,VBMAX
               DO NJ=1,VBMAX
               DO NK=1,VBMAX
                  DO NC=1,NUNOCC
                  DO NB=1,NUNOCC
                     W_ijk_abc_R(:,NB,NC,NI,NJ,NK)=W_ijk_abc_R(:,NB,NC,NI,NJ,NK)+(OVOV_R(NJ,NB,NK,NC))*(T1_R(:,NI,RKIofKI(KI)))
                  ENDDO
                  ENDDO
               ENDDO
               ENDDO
               ENDDO
            ELSE
               DO NI=1,VBMAX
               DO NJ=1,VBMAX
               DO NK=1,VBMAX
                  DO NC=1,NUNOCC
                  DO NB=1,NUNOCC
                     W_ijk_abc(:,NB,NC,NI,NJ,NK)=W_ijk_abc(:,NB,NC,NI,NJ,NK)+CONJG(OVOV(NJ,NB,NK,NC))*CONJG(T1(:,NI,RKIofKI(KI)))
                  ENDDO
                  ENDDO
               ENDDO
               ENDDO
               ENDDO
            ENDIF

         ENDIF

         KQ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,(KA))- &
               WDES%VKPT(:,KI),KPOINTS_FULL)

         IF (KJ==KB) THEN

         IF (.not. LUSESHARED) THEN
            CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
              PW_AI_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KA),1),OC_AI_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KA),1),&
              (VBMAX)*(NUNOCC),&
              PW_AI_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KC),2),OC_AI_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KC),2), &
              VBMAX*NUNOCC,&
              OVOV(1,1,1,1),OVOV_R(1,1,1,1),zero)
            !iakc
         ELSE
            CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
              FTOD_PW_AI(1,1,1,INDEX_KQ_KA(KQ,KA),1),FTOD_OC_AI(1,1,1,INDEX_KQ_KA(KQ,KA),1),&
              (VBMAX)*(NUNOCC),&
              FTOD_PW_AI(1,1,1,INDEX_KQ_KA(KQ,KC),2),FTOD_OC_AI(1,1,1,INDEX_KQ_KA(KQ,KC),2), &
              VBMAX*NUNOCC,&
              OVOV(1,1,1,1),OVOV_R(1,1,1,1),zero)
            !iakc
         ENDIF

            IF (LORBREAL) THEN
               DO NI=1,VBMAX
               DO NJ=1,VBMAX
               DO NK=1,VBMAX
                  DO NC=1,NUNOCC
                  DO NB=1,NUNOCC
                     W_ijk_abc_R(:,NB,NC,NI,NJ,NK)=W_ijk_abc_R(:,NB,NC,NI,NJ,NK)+(OVOV_R(NI,:,NK,NC))*(T1_R(NB,NJ,RKIofKI(KJ)))
                  ENDDO
                  ENDDO
               ENDDO
               ENDDO
               ENDDO
            ELSE
               DO NI=1,VBMAX
               DO NJ=1,VBMAX
               DO NK=1,VBMAX
                  DO NC=1,NUNOCC
                  DO NB=1,NUNOCC
                     W_ijk_abc(:,NB,NC,NI,NJ,NK)=W_ijk_abc(:,NB,NC,NI,NJ,NK)+CONJG(OVOV(NI,:,NK,NC))*CONJG(T1(NB,NJ,RKIofKI(KJ)))
                  ENDDO
                  ENDDO
               ENDDO
               ENDDO
               ENDDO
            ENDIF

         ENDIF


         KQ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,(KA))- &
               WDES%VKPT(:,KI),KPOINTS_FULL)

         IF (KK==KC) THEN

         IF (.not. LUSESHARED) THEN
            CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
              PW_AI_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KA),1),OC_AI_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KA),1),&
              (VBMAX)*(NUNOCC),&
              PW_AI_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KB),2),OC_AI_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KB),2), &
              VBMAX*NUNOCC,&
              OVOV(1,1,1,1),OVOV_R(1,1,1,1),zero)
            !iajb
         ELSE
            CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
              FTOD_PW_AI(1,1,1,INDEX_KQ_KA(KQ,KA),1),FTOD_OC_AI(1,1,1,INDEX_KQ_KA(KQ,KA),1),&
              (VBMAX)*(NUNOCC),&
              FTOD_PW_AI(1,1,1,INDEX_KQ_KA(KQ,KB),2),FTOD_OC_AI(1,1,1,INDEX_KQ_KA(KQ,KB),2), &
              VBMAX*NUNOCC,&
              OVOV(1,1,1,1),OVOV_R(1,1,1,1),zero)
            !iajb
         ENDIF

            IF (LORBREAL) THEN
               DO NI=1,VBMAX
               DO NJ=1,VBMAX
               DO NK=1,VBMAX
                  DO NC=1,NUNOCC
                  DO NB=1,NUNOCC
                     W_ijk_abc_R(:,NB,NC,NI,NJ,NK)=W_ijk_abc_R(:,NB,NC,NI,NJ,NK)+(OVOV_R(NI,:,NJ,NB))*(T1_R(NC,NK,RKIofKI(KC)))
                  ENDDO
                  ENDDO
               ENDDO
               ENDDO
               ENDDO
            ELSE
               DO NI=1,VBMAX
               DO NJ=1,VBMAX
               DO NK=1,VBMAX
                  DO NC=1,NUNOCC
                  DO NB=1,NUNOCC
                     W_ijk_abc(:,NB,NC,NI,NJ,NK)=W_ijk_abc(:,NB,NC,NI,NJ,NK)+CONJG(OVOV(NI,:,NJ,NB))*CONJG(T1(NC,NK,RKIofKI(KC)))
                  ENDDO
                  ENDDO
               ENDDO
               ENDDO
               ENDDO
            ENDIF

         ENDIF

         call date_and_time (values=time_array2)
         EMSCALCWINT=EMSCALCWINT+((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
                time_array2(8)-time_array1(8))

      END SUBROUTINE CALC_Z_ijk_abc


      SUBROUTINE CALC_W_ijk_abc(KI,KJ,KK,KA,KB,KC,WDES,KPOINTS)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE mkpoints
         USE lattice
         USE constant
         USE base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE (kpoints_struct) KPOINTS
         INTEGER ITERATION,NBLOCK,MKIMAX,MKIMIN
         INTEGER KI,KJ,KK,KA,KB,KC,NI,NJ,NK,NM,NA,NB,NC,KQ,KQ_,KF,KM,KMAX,KE,KE2,KM2
         INTEGER ierr
         COMPLEX(q) :: VVOV(NUNOCC,NUNOCC,VBMAX,NUNOCC)
         REAL(q) :: VVOV_R(NUNOCC,NUNOCC,VBMAX,NUNOCC)
         COMPLEX(qs) :: VVVO_S(NUNOCC,NUNOCC,NUNOCC,VBMAX)
         COMPLEX(qs) :: VVOV_S(NUNOCC,NUNOCC,VBMAX,NUNOCC)
         REAL(qs) :: VVVO_S_R(NUNOCC,NUNOCC,NUNOCC,VBMAX)

         COMPLEX(q) :: OOOV(VBMAX,VBMAX,VBMAX,NUNOCC)
         REAL(q) :: OOOV_R(VBMAX,VBMAX,VBMAX,NUNOCC)
         COMPLEX(qs) :: OVOO_S(VBMAX,NUNOCC,VBMAX,VBMAX)
         REAL(qs) :: OVOO_S_R(VBMAX,NUNOCC,VBMAX,VBMAX)
         COMPLEX(qs) :: VVOO_S(NUNOCC,NUNOCC,VBMAX,VBMAX)
         REAL(qs) :: VVOO_S_R(NUNOCC,NUNOCC,VBMAX,VBMAX)
         COMPLEX(q) KW
         INTEGER :: time_array1(8), time_array2(8)

         call date_and_time (values=time_array1)
         KW=(1.0_qs,0.0_qs)*WTKPT

         IF (LORBREAL) THEN
         W_ijk_abc_R=(0.0_qs)
         ELSE
         W_ijk_abc=(0.0_qs,0.0_qs)
         ENDIF

         KQ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,(KC))- &
               WDES%VKPT(:,KK),KPOINTS_FULL)
         KE=KPOINT_IN_FULL_GRID(WDES%VKPT(:,(KB))+ &
               WDES%VKPT(:,KQ),KPOINTS_FULL)
         KQ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,(KI))- &
               WDES%VKPT(:,KA),KPOINTS_FULL)
         KM=KPOINT_IN_FULL_GRID(-WDES%VKPT(:,(KQ))+ &
               WDES%VKPT(:,KB),KPOINTS_FULL)

         ! \sum_e
         KQ=KPOINT_IN_FULL_GRID(-WDES%VKPT(:,(KC))+ &
               WDES%VKPT(:,KK),KPOINTS_FULL)
         IF (.not. LUSESHARED) THEN
            CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
              PW_AB_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KB)),OC_AB_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KB)),&
              (NUNOCC)*(NUNOCC),&
              PW_AI_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KC),2),OC_AI_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KC),2), &
              VBMAX*NUNOCC,&
              VVOV(1,1,1,1),VVOV_R(1,1,1,1),zero)
! ebkc
         ELSE
            CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
              FTOD_PW_AB(1,1,1,INDEX_KQ_KA(KQ,KB)),FTOD_OC_AB(1,1,1,INDEX_KQ_KA(KQ,KB)),&
              (NUNOCC)*(NUNOCC),&
              FTOD_PW_AI(1,1,1,INDEX_KQ_KA(KQ,KC),2),FTOD_OC_AI(1,1,1,INDEX_KQ_KA(KQ,KC),2), &
              VBMAX*NUNOCC,&
              VVOV(1,1,1,1),VVOV_R(1,1,1,1),zero)
! ebkc
         ENDIF
         
         IF (LORBREAL) THEN
            DO NC=1,NUNOCC
            DO NB=1,NUNOCC
               VVVO_S_R(:,NB,NC,:)=(VVOV_R(:,NB,:,NC))
            ENDDO
            ENDDO
            IF (.not. LUSESHARED) THEN
               DO NA=1,NUNOCC
                  VVOO_S_R(:,NA,:,:)=T2_TMP_R(NA,:,:,:,RKIofKI(KI),RKIofKI(KJ),NTMP_KP(KA))
               ENDDO
            ELSE
               DO NA=1,NUNOCC
                  VVOO_S_R(:,NA,:,:)=T2_R(NA,:,:,:,RKIofKI(KI),INDEX_KA_KJ(KA,KJ))
               ENDDO
            ENDIF
            DO NK=1,VBMAX
            DO NJ=1,VBMAX
            DO NI=1,VBMAX
               CALL SGEMM('t','n',(NUNOCC),(NUNOCC*NUNOCC),(NUNOCC), &
                 (1.0_qs), VVOO_S_R(1,1,NI,NJ),(NUNOCC), &
                 VVVO_S_R(1,1,1,NK),(NUNOCC), &
                  (1.0_qs), W_ijk_abc_R(1,1,1,NI,NJ,NK),(NUNOCC))
            ENDDO
            ENDDO
            ENDDO
         ELSE
            DO NC=1,NUNOCC
            DO NB=1,NUNOCC
               VVVO_S(:,NB,NC,:)=(VVOV(:,NB,:,NC))
            ENDDO
            ENDDO
            IF (.not. LUSESHARED) THEN
               DO NA=1,NUNOCC
                  VVOO_S(:,NA,:,:)=T2_TMP(NA,:,:,:,RKIofKI(KI),RKIofKI(KJ),NTMP_KP(KA))
               ENDDO
            ELSE
               DO NA=1,NUNOCC
                  VVOO_S(:,NA,:,:)=T2(NA,:,:,:,RKIofKI(KI),INDEX_KA_KJ(KA,KJ))
               ENDDO
            ENDIF
            DO NK=1,VBMAX
            DO NJ=1,VBMAX                
            DO NI=1,VBMAX
               CALL CGEMM('t','n',(NUNOCC),(NUNOCC*NUNOCC),(NUNOCC), &
                 (1.0_qs,0.0_qs), VVOO_S(1,1,NI,NJ),(NUNOCC), &
                 VVVO_S(1,1,1,NK),(NUNOCC), &
                  (1.0_qs,0.0_qs), W_ijk_abc(1,1,1,NI,NJ,NK),(NUNOCC))
            ENDDO
            ENDDO
            ENDDO
         ENDIF

         KQ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KI)- &
               WDES%VKPT(:,KA),KPOINTS_FULL)

!        -\sum_m

         KQ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,(KM))- &
               WDES%VKPT(:,(KJ)),KPOINTS_FULL)

         IF (.not. LUSESHARED) THEN
            CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
                 FTOD_PW_IJ(1,1,1,RKQofKQ(KQ),RKIofKI(KM)),FTOD_OC_IJ(1,1,1,RKQofKQ(KQ),RKIofKI(KM)),&
                 VBMAX*VBMAX,&
                 PW_AI_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KC),2),OC_AI_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KC),2), &
                 VBMAX*NUNOCC,&
                 OOOV(1,1,1,1),OOOV_R(1,1,1,1),zero)
         !jmkc
         ELSE
            CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
                 FTOD_PW_IJ(1,1,1,RKQofKQ(KQ),RKIofKI(KM)),FTOD_OC_IJ(1,1,1,RKQofKQ(KQ),RKIofKI(KM)),&
                 VBMAX*VBMAX,&
                 FTOD_PW_AI(1,1,1,INDEX_KQ_KA(KQ,KC),2),FTOD_OC_AI(1,1,1,INDEX_KQ_KA(KQ,KC),2), &
                 VBMAX*NUNOCC,&
                 OOOV(1,1,1,1),OOOV_R(1,1,1,1),zero)
         !jmkc
         ENDIF

         IF (.not. LUSESHARED) THEN

         IF (LORBREAL) THEN
            DO NJ=1,VBMAX
            DO NK=1,VBMAX
               OVOO_S_R(:,:,NJ,NK)=(OOOV_R(NJ,:,NK,:))
            ENDDO
            ENDDO
            DO NI=1,VBMAX
            DO NJ=1,VBMAX
            DO NK=1,VBMAX
            DO NM=1,VBMAX
               DO NC=1,NUNOCC
                  W_ijk_abc_R(:,:,NC,NI,NJ,NK)=W_ijk_abc_R(:,:,NC,NI,NJ,NK)-OVOO_S_R(NM,NC,NJ,NK)*T2_TMP_R(:,:,NI,NM,RKIofKI(KI),RKIofKI(KM),NTMP_KP(KA))
               ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ELSE
            DO NJ=1,VBMAX
            DO NK=1,VBMAX
               OVOO_S(:,:,NJ,NK)=(OOOV(NJ,:,NK,:))
            ENDDO
            ENDDO
            DO NI=1,VBMAX
            DO NJ=1,VBMAX
            DO NK=1,VBMAX
            DO NM=1,VBMAX
               DO NC=1,NUNOCC
                  W_ijk_abc(:,:,NC,NI,NJ,NK)=W_ijk_abc(:,:,NC,NI,NJ,NK)-OVOO_S(NM,NC,NJ,NK)*T2_TMP(:,:,NI,NM,RKIofKI(KI),RKIofKI(KM),NTMP_KP(KA))
               ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDIF

         ELSE

         IF (LORBREAL) THEN
            DO NJ=1,VBMAX
            DO NK=1,VBMAX
               OVOO_S_R(:,:,NJ,NK)=(OOOV_R(NJ,:,NK,:))
            ENDDO
            ENDDO
            DO NI=1,VBMAX
            DO NJ=1,VBMAX
            DO NK=1,VBMAX
            DO NM=1,VBMAX
               DO NC=1,NUNOCC
                  W_ijk_abc_R(:,:,NC,NI,NJ,NK)=W_ijk_abc_R(:,:,NC,NI,NJ,NK)-OVOO_S_R(NM,NC,NJ,NK)*T2_R(:,:,NI,NM,RKIofKI(KI),INDEX_KA_KJ(KA,KM))
               ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ELSE
            DO NJ=1,VBMAX
            DO NK=1,VBMAX
               OVOO_S(:,:,NJ,NK)=(OOOV(NJ,:,NK,:))
            ENDDO
            ENDDO
            DO NI=1,VBMAX
            DO NJ=1,VBMAX
            DO NK=1,VBMAX
            DO NM=1,VBMAX
               DO NC=1,NUNOCC
                  W_ijk_abc(:,:,NC,NI,NJ,NK)=W_ijk_abc(:,:,NC,NI,NJ,NK)-OVOO_S(NM,NC,NJ,NK)*T2(:,:,NI,NM,RKIofKI(KI),INDEX_KA_KJ(KA,KM))
               ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDIF

         ENDIF
         call date_and_time (values=time_array2)
         EMSCALCWINT=EMSCALCWINT+((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
                time_array2(8)-time_array1(8))

      END SUBROUTINE CALC_W_ijk_abc


      SUBROUTINE CALC_W_ijk_abc_kq(KI,KJ,KK,KA,KB,KC,NG,WDES,KPOINTS)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE mkpoints
         USE lattice
         USE constant
         USE base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE (kpoints_struct) KPOINTS
         INTEGER ITERATION,NBLOCK,MKIMAX,MKIMIN
         INTEGER KI,KJ,KK,KA,KB,KC,NI,NJ,NK,NM,NA,NB,NC,KQ,KQ_,KF,KM,KMAX,KE,KE2,KM2,NG
         INTEGER ierr
         COMPLEX(q) :: VVOV(NUNOCC,NUNOCC,VBMAX,NUNOCC)
         REAL(q) :: VVOV_R(NUNOCC,NUNOCC,VBMAX,NUNOCC)
         COMPLEX(qs) :: VVVO_S(NUNOCC,NUNOCC,NUNOCC,VBMAX)
         COMPLEX(qs) :: VVOV_S(NUNOCC,NUNOCC,VBMAX,NUNOCC)
         REAL(qs) :: VVVO_S_R(NUNOCC,NUNOCC,NUNOCC,VBMAX)

         COMPLEX(q) :: OOOV(VBMAX,VBMAX,VBMAX,NUNOCC)
         REAL(q) :: OOOV_R(VBMAX,VBMAX,VBMAX,NUNOCC)
         COMPLEX(qs) :: OVOO_S(VBMAX,NUNOCC,VBMAX,VBMAX)
         REAL(qs) :: OVOO_S_R(VBMAX,NUNOCC,VBMAX,VBMAX)
         COMPLEX(qs) :: VVOO_S(NUNOCC,NUNOCC,VBMAX,VBMAX)
         REAL(qs) :: VVOO_S_R(NUNOCC,NUNOCC,VBMAX,VBMAX)
         COMPLEX(q) KW
         INTEGER :: time_array1(8), time_array2(8)
         COMPLEX(q) :: MPW_AB_TMP(NUNOCC,NUNOCC), MPW_AI_TMP(VBMAX,NUNOCC),MFTOD_PW_IJ(VBMAX,VBMAX)

         call date_and_time (values=time_array1)
         KW=(1.0_qs,0.0_qs)*WTKPT

         IF (LORBREAL) THEN
         ELSE
         W_ijk_abc_kq=(0.0_qs,0.0_qs)
         ENDIF

         KQ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,(KC))- &
               WDES%VKPT(:,KK),KPOINTS_FULL)
         KE=KPOINT_IN_FULL_GRID(WDES%VKPT(:,(KB))+ &
               WDES%VKPT(:,KQ),KPOINTS_FULL)
         KQ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,(KI))- &
               WDES%VKPT(:,KA),KPOINTS_FULL)
         KM=KPOINT_IN_FULL_GRID(-WDES%VKPT(:,(KQ))+ &
               WDES%VKPT(:,KB),KPOINTS_FULL)

         ! \sum_e
         KQ=KPOINT_IN_FULL_GRID(-WDES%VKPT(:,(KC))+ &
               WDES%VKPT(:,KK),KPOINTS_FULL)
         IF (.not. LUSESHARED) THEN
            MPW_AB_TMP(:,:)=PW_AB_TMP(NG,:,:,RKQofKQ(KQ),NTMP_KP(KB))
            MPW_AI_TMP(:,:)=PW_AI_TMP(NG,:,:,RKQofKQ(KQ),NTMP_KP(KC),2)
            CALL CONTR_FTOD(1,NHVECTOR, &
              MPW_AB_TMP(1,1),OC_AB_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KB)),&
              (NUNOCC)*(NUNOCC),&
              MPW_AI_TMP(1,1),OC_AI_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KC),2), &
              VBMAX*NUNOCC,&
              VVOV(1,1,1,1),VVOV_R(1,1,1,1),zero)
! ebkc
         ELSE
            MPW_AB_TMP(:,:)=FTOD_PW_AB(NG,:,:,INDEX_KQ_KA(KQ,KB))
            MPW_AI_TMP(:,:)=FTOD_PW_AI(NG,:,:,INDEX_KQ_KA(KQ,KC),2)
            CALL CONTR_FTOD(1,NHVECTOR, &
              MPW_AB_TMP(1,1),FTOD_OC_AB(1,1,1,INDEX_KQ_KA(KQ,KB)),&
              (NUNOCC)*(NUNOCC),&
              MPW_AI_TMP(1,1),FTOD_OC_AI(1,1,1,INDEX_KQ_KA(KQ,KC),2), &
              VBMAX*NUNOCC,&
              VVOV(1,1,1,1),VVOV_R(1,1,1,1),zero)
! ebkc
         ENDIF
         
         IF (LORBREAL) THEN
         ELSE
            DO NC=1,NUNOCC
            DO NB=1,NUNOCC
               VVVO_S(:,NB,NC,:)=(VVOV(:,NB,:,NC))
            ENDDO
            ENDDO
            IF (.not. LUSESHARED) THEN
               DO NA=1,NUNOCC
                  VVOO_S(:,NA,:,:)=T2_TMP(NA,:,:,:,RKIofKI(KI),RKIofKI(KJ),NTMP_KP(KA))
               ENDDO
            ELSE
               DO NA=1,NUNOCC
                  VVOO_S(:,NA,:,:)=T2(NA,:,:,:,RKIofKI(KI),INDEX_KA_KJ(KA,KJ))
               ENDDO
            ENDIF
            DO NK=1,VBMAX
            DO NJ=1,VBMAX                
            DO NI=1,VBMAX
               CALL CGEMM('t','n',(NUNOCC),(NUNOCC*NUNOCC),(NUNOCC), &
                 (1.0_qs,0.0_qs), VVOO_S(1,1,NI,NJ),(NUNOCC), &
                 VVVO_S(1,1,1,NK),(NUNOCC), &
                  (1.0_qs,0.0_qs), W_ijk_abc_kq(1,1,1,NI,NJ,NK,RKQofKQ(KQ)),(NUNOCC))
            ENDDO
            ENDDO
            ENDDO
         ENDIF

         KQ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KI)- &
               WDES%VKPT(:,KA),KPOINTS_FULL)

!        -\sum_m

         KQ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,(KM))- &
               WDES%VKPT(:,(KJ)),KPOINTS_FULL)

         IF (.not. LUSESHARED) THEN
            MFTOD_PW_IJ(:,:)=FTOD_PW_IJ(NG,:,:,RKQofKQ(KQ),RKIofKI(KM))
            MPW_AI_TMP(:,:)=PW_AI_TMP(NG,:,:,RKQofKQ(KQ),NTMP_KP(KC),2)
            CALL CONTR_FTOD(1,NHVECTOR, &
                 MFTOD_PW_IJ(1,1),FTOD_OC_IJ(1,1,1,RKQofKQ(KQ),RKIofKI(KM)),&
                 VBMAX*VBMAX,&
                 MPW_AI_TMP(1,1),OC_AI_TMP(1,1,1,RKQofKQ(KQ),NTMP_KP(KC),2), &
                 VBMAX*NUNOCC,&
                 OOOV(1,1,1,1),OOOV_R(1,1,1,1),zero)
         !jmkc
         ELSE
            MFTOD_PW_IJ(:,:)=FTOD_PW_IJ(NG,:,:,RKQofKQ(KQ),RKIofKI(KM))
            MPW_AI_TMP(:,:)=FTOD_PW_AI(NG,:,:,INDEX_KQ_KA(KQ,KC),2)
            CALL CONTR_FTOD(1,NHVECTOR, &
                 MFTOD_PW_IJ(1,1),FTOD_OC_IJ(1,1,1,RKQofKQ(KQ),RKIofKI(KM)),&
                 VBMAX*VBMAX,&
                 MPW_AI_TMP(1,1),FTOD_OC_AI(1,1,1,INDEX_KQ_KA(KQ,KC),2), &
                 VBMAX*NUNOCC,&
                 OOOV(1,1,1,1),OOOV_R(1,1,1,1),zero)
         !jmkc
         ENDIF

         IF (.not. LUSESHARED) THEN

         IF (LORBREAL) THEN
         ELSE
            DO NJ=1,VBMAX
            DO NK=1,VBMAX
               OVOO_S(:,:,NJ,NK)=(OOOV(NJ,:,NK,:))
            ENDDO
            ENDDO
            DO NI=1,VBMAX
            DO NJ=1,VBMAX
            DO NK=1,VBMAX
            DO NM=1,VBMAX
               DO NC=1,NUNOCC
                  W_ijk_abc_kq(:,:,NC,NI,NJ,NK,RKQofKQ(KQ))=W_ijk_abc_kq(:,:,NC,NI,NJ,NK,RKQofKQ(KQ))-OVOO_S(NM,NC,NJ,NK)*T2_TMP(:,:,NI,NM,RKIofKI(KI),RKIofKI(KM),NTMP_KP(KA))
               ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDIF

         ELSE

         IF (LORBREAL) THEN
         ELSE
            DO NJ=1,VBMAX
            DO NK=1,VBMAX
               OVOO_S(:,:,NJ,NK)=(OOOV(NJ,:,NK,:))
            ENDDO
            ENDDO
            DO NI=1,VBMAX
            DO NJ=1,VBMAX
            DO NK=1,VBMAX
            DO NM=1,VBMAX
               DO NC=1,NUNOCC
                  W_ijk_abc_kq(:,:,NC,NI,NJ,NK,RKQofKQ(KQ))=W_ijk_abc_kq(:,:,NC,NI,NJ,NK,RKQofKQ(KQ))-OVOO_S(NM,NC,NJ,NK)*T2(:,:,NI,NM,RKIofKI(KI),INDEX_KA_KJ(KA,KM))
               ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDIF

         ENDIF
         call date_and_time (values=time_array2)
         EMSCALCWINT=EMSCALCWINT+((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
                time_array2(8)-time_array1(8))

      END SUBROUTINE CALC_W_ijk_abc_kq



      Subroutine BCAST2ALL_T2_KA(KA,WDES)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: KA,KJ
         integer :: time_array1(8),time_array2(8),ems

            call BLACS_PINFO(ME,PROCS)
            call date_and_time (values=time_array1)
         DO KJ=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ)==0)) CYCLE
            IF (LORBREAL) THEN
               IF (PROCS_KA_KJ(KA,KJ)==ME) THEN
                  CALL SGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC),(VBMAX*VBMAX*REALNKPTS), &
                     T2_R(1,1,1,1,1,INDEX_KA_KJ(KA,KJ)), (NUNOCC*NUNOCC))
                  T2_TMP_R(:,:,:,:,:,RKIofKI(KJ),1)=T2_R(:,:,:,:,:,INDEX_KA_KJ(KA,KJ))
               ELSE
                  CALL SGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC),(VBMAX*VBMAX*REALNKPTS), &
                     T2_TMP_R(1,1,1,1,1,RKIofKI(KJ),1), (NUNOCC*NUNOCC),0,PROCS_KA_KJ(KA,KJ))
               ENDIF
            ELSE
               IF (PROCS_KA_KJ(KA,KJ)==ME) THEN
                  CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC),(VBMAX*VBMAX*REALNKPTS), &
                     T2(1,1,1,1,1,INDEX_KA_KJ(KA,KJ)), (NUNOCC*NUNOCC))
                  T2_TMP(:,:,:,:,:,RKIofKI(KJ),1)=T2(:,:,:,:,:,INDEX_KA_KJ(KA,KJ))
               ELSE
                  CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC),(VBMAX*VBMAX*REALNKPTS), &
                     T2_TMP(1,1,1,1,1,RKIofKI(KJ),1), (NUNOCC*NUNOCC),0,PROCS_KA_KJ(KA,KJ))
               ENDIF
            ENDIF

            call date_and_time (values=time_array2)
            ems=((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
                time_array2(8)-time_array1(8))
            TET2S=TET2S+ems
         ENDDO

      END Subroutine BCAST2ALL_T2_KA

      Subroutine BCAST2ALL_FTOD_AB_KA(KA,WDES)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: KA,KQ
         integer :: time_array1(8),time_array2(8),ems

            call BLACS_PINFO(ME,PROCS)
            call date_and_time (values=time_array1)
         DO KQ=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0.0_q)) CYCLE
            IF (PROCS_KQ_KA(KQ,KA)==ME) THEN
               CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR,(NUNOCC)*(NUNOCC), &
                  FTOD_PW_AB(1,1,1,INDEX_KQ_KA(KQ,KA)), NGVECTOR)
               PW_AB_TMP(:,:,:,RKQofKQ(KQ),1)=FTOD_PW_AB(:,:,:,INDEX_KQ_KA(KQ,KA))
            ELSE
               CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR,(NUNOCC)*(NUNOCC), &
                  PW_AB_TMP(1,1,1,RKQofKQ(KQ),1), NGVECTOR, 0, PROCS_KQ_KA(KQ,KA))
            ENDIF
            IF (ASSOCIATED(H)) THEN

            ENDIF
            call date_and_time (values=time_array2)
            ems=((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
                time_array2(8)-time_array1(8))
            TET2S=TET2S+ems
!            write(*,*)'done'
         ENDDO

      END Subroutine BCAST2ALL_FTOD_AB_KA

      Subroutine BCAST2ALL_FTOD_AI_KA(KA,WDES)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: KA,KQ,mncc
         integer :: time_array1(8),time_array2(8),ems

            call BLACS_PINFO(ME,PROCS)
            call date_and_time (values=time_array1)
         DO mncc=1,2
         DO KQ=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0.0_q)) CYCLE
            IF (PROCS_KQ_KA(KQ,KA)==ME) THEN
               CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR,(NUNOCC)*(VBMAX), &
                  FTOD_PW_AI(1,1,1,INDEX_KQ_KA(KQ,KA),mncc), NGVECTOR)
               PW_AI_TMP(:,:,:,RKQofKQ(KQ),1,mncc)=FTOD_PW_AI(:,:,:,INDEX_KQ_KA(KQ,KA),mncc)
            ELSE
               CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR,(NUNOCC)*(VBMAX), &
                  PW_AI_TMP(1,1,1,RKQofKQ(KQ),1,mncc), NGVECTOR, 0, PROCS_KQ_KA(KQ,KA))
            ENDIF
            IF (ASSOCIATED(H)) THEN

            ENDIF
            call date_and_time (values=time_array2)
            ems=((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
                time_array2(8)-time_array1(8))
            TET2S=TET2S+ems
         ENDDO
         ENDDO

      END Subroutine BCAST2ALL_FTOD_AI_KA


      Subroutine BCAST2ALL_T2(WDES)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: KA,KJ,MY_TMP
         COMPLEX(qs) :: MTMP(NUNOCC,NUNOCC,VBMAX,VBMAX,WDES%NKPTS)
         REAL(qs) :: MTMP_R(NUNOCC,NUNOCC,VBMAX,VBMAX,WDES%NKPTS)
         integer :: time_array1(8),time_array2(8),ems

            call BLACS_PINFO(ME,PROCS)
            call date_and_time (values=time_array1)
         IF (LORBREAL) THEN
            DO KA=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0.0_q)) CYCLE
            DO KJ=1,WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ)==0.0_q)) CYCLE
               IF (PROCS_KA_KJ(KA,KJ)==ME) THEN
                  CALL SGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC),(VBMAX*VBMAX*REALNKPTS), &
                     T2_R(1,1,1,1,1,INDEX_KA_KJ(KA,KJ)), (NUNOCC*NUNOCC))
                  DO MY_TMP=1,MY_NKPTS*2
                     IF (KP_TMP(MY_TMP+1)==KA) THEN
                        T2_TMP_R(:,:,:,:,:,RKIofKI(KJ),MY_TMP+1)=T2_R(:,:,:,:,:,INDEX_KA_KJ(KA,KJ))
                     ENDIF
                  ENDDO
               ELSE
                  CALL SGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC),(VBMAX*VBMAX*REALNKPTS), &
                     MTMP_R(1,1,1,1,1), (NUNOCC*NUNOCC),0,PROCS_KA_KJ(KA,KJ))
                  DO MY_TMP=1,MY_NKPTS*2
                     IF (KP_TMP(MY_TMP+1)==KA) THEN
                        T2_TMP_R(:,:,:,:,:,RKIofKI(KJ),MY_TMP+1)=MTMP_R(:,:,:,:,:)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            ENDDO
         ELSE
            DO KA=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0.0_q)) CYCLE
            DO KJ=1,WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ)==0.0_q)) CYCLE
               IF (PROCS_KA_KJ(KA,KJ)==ME) THEN
                  CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC),(VBMAX*VBMAX*REALNKPTS), &
                     T2(1,1,1,1,1,INDEX_KA_KJ(KA,KJ)), (NUNOCC*NUNOCC))
                  DO MY_TMP=1,MY_NKPTS*2
                     IF (KP_TMP(MY_TMP+1)==KA) THEN
                        T2_TMP(:,:,:,:,:,RKIofKI(KJ),MY_TMP+1)=T2(:,:,:,:,:,INDEX_KA_KJ(KA,KJ))
                     ENDIF
                  ENDDO
               ELSE
                  CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC),(VBMAX*VBMAX*REALNKPTS), &
                     MTMP(1,1,1,1,1), (NUNOCC*NUNOCC),0,PROCS_KA_KJ(KA,KJ))
                  DO MY_TMP=1,MY_NKPTS*2
                     IF (KP_TMP(MY_TMP+1)==KA) THEN
                        T2_TMP(:,:,:,:,:,RKIofKI(KJ),MY_TMP+1)=MTMP(:,:,:,:,:)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            ENDDO
         ENDIF

         call date_and_time (values=time_array2)
         ems=((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
             time_array2(8)-time_array1(8))
         TET2S=TET2S+ems

      END Subroutine BCAST2ALL_T2

      Subroutine BCAST2ALL_FTOD_AB(WDES)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: KA,KJ,MY_TMP,KQ
         COMPLEX(q) :: MTMP(NGVECTOR,NUNOCC,NUNOCC)
         REAL(q) :: MTMP_R(NGVECTOR,NUNOCC,NUNOCC)
         integer :: time_array1(8),time_array2(8),ems

            call BLACS_PINFO(ME,PROCS)
            call date_and_time (values=time_array1)
            DO KA=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0.0_q)) CYCLE
            DO KQ=1,WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0.0_q)) CYCLE
               IF (PROCS_KQ_KA(KQ,KA)==ME) THEN
                  CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR),(NUNOCC*NUNOCC), &
                     FTOD_PW_AB(1,1,1,INDEX_KQ_KA(KQ,KA)), (NGVECTOR))
                  DO MY_TMP=1,MY_NKPTS*2
                     IF (KP_TMP(MY_TMP+1)==KA) THEN
                        PW_AB_TMP(:,:,:,RKQofKQ(KQ),MY_TMP+1)=FTOD_PW_AB(:,:,:,INDEX_KQ_KA(KQ,KA))
                     ENDIF
                  ENDDO
               ELSE
                  CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR),(NUNOCC)*(NUNOCC), &
                     MTMP(1,1,1), (NGVECTOR),0,PROCS_KQ_KA(KQ,KA))
                  DO MY_TMP=1,MY_NKPTS*2
                     IF (KP_TMP(MY_TMP+1)==KA) THEN
                        PW_AB_TMP(:,:,:,RKQofKQ(KQ),MY_TMP+1)=MTMP(:,:,:)
                     ENDIF
                  ENDDO
               ENDIF
               IF (ASSOCIATED(H)) THEN
               ENDIF
            ENDDO
            ENDDO

         call date_and_time (values=time_array2)
         ems=((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
             time_array2(8)-time_array1(8))
         TET2S=TET2S+ems

      END Subroutine BCAST2ALL_FTOD_AB

      Subroutine BCAST2ALL_FTOD_AI(WDES)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: KA,KJ,MY_TMP,KQ
         COMPLEX(q) :: MTMP(NGVECTOR,VBMAX,NUNOCC)
         REAL(q) :: MTMP_R(NGVECTOR,VBMAX,NUNOCC)
         integer :: time_array1(8),time_array2(8),ems,mncc

            call BLACS_PINFO(ME,PROCS)
            call date_and_time (values=time_array1)
            DO KA=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0.0_q)) CYCLE
            DO KQ=1,WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0.0_q)) CYCLE
               IF (PROCS_KQ_KA(KQ,KA)==ME) THEN
                  DO mncc=1,2
                     CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR),(NUNOCC*VBMAX), &
                        FTOD_PW_AI(1,1,1,INDEX_KQ_KA(KQ,KA),mncc), (NGVECTOR))
                     DO MY_TMP=1,MY_NKPTS*2
                        IF (KP_TMP(MY_TMP+1)==KA) THEN
                           PW_AI_TMP(:,:,:,RKQofKQ(KQ),MY_TMP+1,mncc)=FTOD_PW_AI(:,:,:,INDEX_KQ_KA(KQ,KA),mncc)
                        ENDIF
                     ENDDO
                  ENDDO
               ELSE
                  DO mncc=1,2
                     CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR),(NUNOCC)*(VBMAX), &
                        MTMP(1,1,1), (NGVECTOR),0,PROCS_KQ_KA(KQ,KA))
                     DO MY_TMP=1,MY_NKPTS*2
                        IF (KP_TMP(MY_TMP+1)==KA) THEN
                           PW_AI_TMP(:,:,:,RKQofKQ(KQ),MY_TMP+1,mncc)=MTMP(:,:,:)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
               IF (ASSOCIATED(H)) THEN
               ENDIF
            ENDDO
            ENDDO

         call date_and_time (values=time_array2)
         ems=((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
             time_array2(8)-time_array1(8))
         TET2S=TET2S+ems

      END Subroutine BCAST2ALL_FTOD_AI

    SUBROUTINE CHECK_SHIFTED_KPOINTS(WDES,W,KPOINTS)
         USE constant
         USE full_kpoints
         USE mkpoints
         USE wave
         implicit NONE
         TYPE(wavedes) WDES
         TYPE(wavespin) W
         TYPE (kpoints_struct) KPOINTS
         integer :: MKI,KI
         LOGICAL :: GAMMA_FOUND
         LOGICAL :: SECOND_GAMMA
         REAL(q) :: SECOND_GAMMA_WEIGHT
         REAL(q) :: GAMMA_WEIGHT
         INTEGER :: NRKI, NRKQ
         
         GAMMA_FOUND=.FALSE.
         SHIFTED_KPOINTS=.TRUE.
         SECOND_GAMMA=.FALSE.
         do MKI=1,WDES%NKPTS
            if ( (ABS(W%WDES%VKPT(1,MKI))==0.0_q) .and. (ABS(W%WDES%VKPT(2,MKI))==0.0_q) .and. (ABS(W%WDES%VKPT(3,MKI))==0.0_q)) then
               IF (GAMMA_FOUND) SECOND_GAMMA=.TRUE.
               IF (SECOND_GAMMA) SECOND_GAMMA_WEIGHT=W%WDES%WTKPT(MKI)
               IF (.NOT. GAMMA_FOUND) THEN
                  GAMMA_WEIGHT=W%WDES%WTKPT(MKI)
                  GAMMA_FOUND=.true.
               ENDIF
            endif
         enddo
         if ((GAMMA_FOUND) .and. (GAMMA_WEIGHT/=0.0_q)) THEN
            SHIFTED_KPOINTS=.FALSE.
         endif
         IF (SECOND_GAMMA) SHIFTED_KPOINTS=.TRUE.

         REALNKPTS=WDES%NKPTS
         IF (SHIFTED_KPOINTS) REALNKPTS=WDES%NKPTS/2

         ALLOCATE(RKQofKQ(WDES%NKPTS))
         ALLOCATE(RKIofKI(WDES%NKPTS))
         ALLOCATE(KQofRKQ(REALNKPTS))
         ALLOCATE(KIofRKI(REALNKPTS))
         NRKI=0
         NRKQ=0
         RKQofKQ=0
         RKIofKI=0
         KQofRKQ=0
         KIofRKI=0
         IF (SHIFTED_KPOINTS) THEN
            DO KI=1,WDES%NKPTS
               IF ((WDES%WTKPT(KI)/=0.00_q))  THEN
                  NRKI=NRKI+1
                  RKIofKI(KI)=NRKI
                  KIofRKI(NRKI)=KI
               ENDIF
               IF ((WDES%WTKPT(KI)==0.00_q))  THEN
                  NRKQ=NRKQ+1
                  RKQofKQ(KI)=NRKQ
                  KQofRKQ(NRKQ)=KI
               ENDIF

            ENDDO
         ELSE
            DO KI=1,WDES%NKPTS
               RKIofKI(KI)=KI
               KIofRKI(KI)=KI
               RKQofKQ(KI)=KI
               KQofRKQ(KI)=KI
            ENDDO
         ENDIF

        WTKPT=1.0_q/REALNKPTS 
         
    END SUBROUTINE CHECK_SHIFTED_KPOINTS

    SUBROUTINE SET_GVEC(WDES1, LATT_CUR, POTFAK, GVX, GVY, GVZ ,ENCUT, ENCUTSOFT)
      USE constant
      USE wave
      USE lattice
      USE full_kpoints
      IMPLICIT NONE

      TYPE (wavedes1) WDES1
      TYPE (latt) :: LATT_CUR
      REAL(q) :: POTFAK(WDES1%NGVECTOR)
      REAL(q) :: GVX(WDES1%NGVECTOR)
      REAL(q) :: GVY(WDES1%NGVECTOR)
      REAL(q) :: GVZ(WDES1%NGVECTOR)
      REAL(q), OPTIONAL :: ENCUT, ENCUTSOFT
   ! local
      INTEGER    NI,NP
      REAL(q) :: DKX,DKY,DKZ,GX,GY,GZ,GSQU,GSQUP,GQUAD,SCALE,  RHOEFF, FOMEGAP2
      REAL(q) :: E, FSG_TMP
      REAL(q) :: Q1,Q2

      ! input needed for attenuated Coulomb kernel
      IF (PRESENT(ENCUT).AND. PRESENT(ENCUTSOFT)) THEN
         Q1=SQRT(ENCUTSOFT/HSQDTM)
         Q2=SQRT(ENCUT/HSQDTM)
      ENDIF

      ! effective electron density in a.u.^-3
      RHOEFF=(HFSCREEN*HFSCREEN*AUTOA*AUTOA/(4._q*EXP(LOG(3._q/PI)/3._q)))**3
      ! plasma frequency
      FOMEGAP2=16._q*PI*RHOEFF/(AUTOA*AUTOA*AUTOA*AUTOA)

      SCALE=EDEPS/LATT_CUR%OMEGA/TPI**2

! 
!  Calculate the Madelungterm for different potentials (Yukawa, F12, Coulomb)
!

      IF (KPOINTS_FULL%WTKPT(1)==0) THEN
         WRITE(*,*) 'internal error in  SET_GFAC_WAVEFUN: division by zero in SCALEFSG, and FSG anyway most likely wrong'
         WRITE(*,*) ' needs to be fixed in  CALCULATE_LOCAL_FIELD_PREPARE as well'
         STOP
      ENDIF

      DKX=WDES1%VKPT(1)*LATT_CUR%B(1,1)+ &
          WDES1%VKPT(2)*LATT_CUR%B(1,2)+ &
          WDES1%VKPT(3)*LATT_CUR%B(1,3)
      DKY=WDES1%VKPT(1)*LATT_CUR%B(2,1)+ &
          WDES1%VKPT(2)*LATT_CUR%B(2,2)+ &
          WDES1%VKPT(3)*LATT_CUR%B(2,3)
      DKZ=WDES1%VKPT(1)*LATT_CUR%B(3,1)+ &
          WDES1%VKPT(2)*LATT_CUR%B(3,2)+ &
          WDES1%VKPT(3)*LATT_CUR%B(3,3)

      NP=WDES1%NGVECTOR

      DO NI=1,NP

         GX=(WDES1%IGX(NI)*LATT_CUR%B(1,1)+WDES1%IGY(NI)* &
              LATT_CUR%B(1,2)+WDES1%IGZ(NI)*LATT_CUR%B(1,3))
         GY=(WDES1%IGX(NI)*LATT_CUR%B(2,1)+WDES1%IGY(NI)* &
              LATT_CUR%B(2,2)+WDES1%IGZ(NI)*LATT_CUR%B(2,3))
         GZ=(WDES1%IGX(NI)*LATT_CUR%B(3,1)+WDES1%IGY(NI)* &
              LATT_CUR%B(3,2)+WDES1%IGZ(NI)*LATT_CUR%B(3,3))

         POTFAK(NI)=((DKX+GX)**2+(DKY+GY)**2+(DKZ+GZ)**2)**0.5
         GVX(NI)=DKX+GX
         GVY(NI)=DKY+GY
         GVZ(NI)=DKZ+GZ

      ENDDO

    END SUBROUTINE SET_GVEC

#endif // gammareal
#endif // scaLAPACK

END MODULE bracketsT
EOF
    echo "Adding file  bracketst.F " 
cp .tmp_vaspcc4s_patch  bracketst.F
 
#######################
# new file   cc4s.F
#######################
cat > .tmp_vaspcc4s_patch<<"EOF"
!*********************************************************************************
!
! MODULE :: CC4S
! AUTHOR :: Andreas Grüneis ( andreas.grueneis@tuwien.ac.at )
!
! DESCRIPTION: 
!
! THIS MODULE PREPARES ALL REQUIRED OUTPUT FOR THE CC4S CODE TO PERFORM CCSD(T)-LEVEL CALCULATIONS.
!
!*********************************************************************************
#include "symbol.inc"
MODULE CC4S
   USE prec
   USE fock
   USE chi_base
   USE lattice
   USE wpot
   USE hamil
   USE pawm
   USE us
   USE pot
   IMPLICIT NONE

#ifdef scaLAPACK

!**********************************************************************
   PRIVATE :: CALC_CVERTEX,REDISTRIBUTE_CVERTEX_1D,INIT_BLACS_1D,INIT_BLACS_2D
!*********************************************************************************
! Allocation of Coulomb-Vertex Arrays
!*********************************************************************************
   !CoulombVertex and SingularVectors
   COMPLEX(q) , ALLOCATABLE :: CVERTEX_2D(:,:,:,:,:,:), SVEC_2D(:,:), SVM_2D(:,:), SVM_TMP_2D(:,:)
   COMPLEX(q) , ALLOCATABLE :: TSVEC_2D(:,:)
   COMPLEX(q) , ALLOCATABLE :: TSVEC_1D(:,:)
  ! CoulombVertex in optimized auxilliary field basis
   COMPLEX(q) , ALLOCATABLE :: OPTCVERTEX_1D(:,:,:,:,:,:), TMPCVERTEX_2D(:,:), TMPCVERTEX_1D(:,:)
   DOUBLE PRECISION, ALLOCATABLE :: SVALUES(:)
   INTEGER, PRIVATE, SAVE :: NGVECTOR, NOPTAUX, NOPTAUX_EDIFF
   !Scalapack array descriptors
   integer, dimension(9)   :: DESC_CVERTEX_1D
   integer, dimension(9)   :: DESC_CVERTEX_2D
   integer, dimension(9)   :: DESC_TMPCVERTEX_2D
   integer, dimension(9)   :: DESC_OPTCVERTEX_1D
   integer, dimension(9)   :: DESC_TMPCVERTEX_1D
   integer, dimension(9)   :: DESC_SVM_2D
   integer, dimension(9)   :: DESC_TSVEC_2D
   integer, dimension(9)   :: DESC_TSVEC_1D
   !BLACS related variables and function
   !PROCS.. number of processors, ME... local processor number, NPROW... number of rows in process grid
   !NPCOL... number of columns in process grid, myrow,mycol... my coordinates in process grid
   !mb... blocking size of rows for block cyclic distribution
   !nb... blocking size of columns for block cyclic distribution
   INTEGER, ALLOCATABLE :: PRO_NI(:)
   INTEGER, PRIVATE, SAVE :: PROCS,ME,NPROW_1D,NPCOL_1D,MYROW_1D,MYCOL_1D,CONTXT_1D
   INTEGER, PRIVATE, SAVE :: NPROW_2D, NPCOL_2D, MYROW_2D, MYCOL_2D, CONTXT_2D
   INTEGER, PRIVATE, SAVE :: VBMAX, VBMAXUP, VBMAXDN
   INTEGER, EXTERNAL :: NUMROC, BLACS_PNUM
   LOGICAL :: DUMPFNO
   TYPE(one_center_handle), POINTER, PRIVATE, SAVE :: H
!*********************************************************************************
!
! For shifted (non-gamma-centered k-meshes) we need a set of the following
! work arrays to handle all k-points in the employed k-meshes
! KI refers to k-point of one-electron state
! KQ refers to momentum transfer of two one-electron states
! N.B.: {KQ} will always be gamma-centered
!       {KI} is defined by the KPOINTS file in VASP
!
!*********************************************************************************
   LOGICAL :: SHIFTED_KPOINTS
   REAL(q) :: WTKPT
   INTEGER :: REALNKPTS
   INTEGER, ALLOCATABLE :: RKQofKQ(:), RKIofKI(:)
   INTEGER, ALLOCATABLE :: KQofRKQ(:), KIofRKI(:)
!*********************************************************************************
! NFREEZE specifies the number of frozen occupied orbitals, ignored for embedding calculations
! and controlled by the NBANDSLOW flag
!*********************************************************************************
   INTEGER :: NFREEZE
   INTEGER :: NBANDSDUMP
   INTEGER :: NBANDSDUMPLOC, NFREEZELOC
   INTEGER :: SRCPROC
! intermediate lists not used anymore but needed for all FTODDUMP outputs 
   INTEGER, ALLOCATABLE :: N_ORD(:)
   INTEGER, ALLOCATABLE :: SP_ORD(:)
! At some point we will support non-canonical references in CC4S and then we need to compute
! the Fock matrix
   GDEF,ALLOCATABLE ::  FOCKM(:,:,:,:) ! Hamiltonian
!*********************************************************************************
!
! For the calculation of the electronic transition structure factor we need
! the following work arrays.
!
!*********************************************************************************
   COMPLEX(q) , ALLOCATABLE, PRIVATE, SAVE :: SINGULAR_VECTORS(:,:)
   REAL(q) , ALLOCATABLE, PRIVATE, SAVE :: GVEC_FULL(:,:,:)
   COMPLEX(q) , ALLOCATABLE, PRIVATE, SAVE :: POTFAK_FULL(:,:)
   INTEGER, DIMENSION(9)   :: DESC_SV
!*********************************************************************************
! Some control flags 
!*********************************************************************************
   LOGICAL :: LLFTODDUMP, LSKIPFOCK
!*********************************************************************************
! The D2PAW_VVOO array stores the Delta Kernel two-electron integrals
!*********************************************************************************
   GDEF, ALLOCATABLE :: D2PAW_VVOO(:, :, :,:,:,:,:) ! delta kernel integrals


   CONTAINS

!***********************************************************************
!
! Main CC4S interface routine
!
!***********************************************************************
      SUBROUTINE CC4S_INTERFACE(P,WDES,W,LATT_CUR,LATT_INI,T_INFO,IO,KPOINTS,WGW, ENCUTGW, ENCUTGWSOFT, LMAXMP2, INFO, &
                           HAMILTONIAN,SYMM,GRID,NONLR_S,NONL_S,LMDIM,CDIJ,CQIJ,SV,E, &
                           GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C,&
                           CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
                           CRHODE, N_MIX_PAW, RHOLM, CHDEN)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE base
         USE mkpoints
         USE mpimy
         USE ini
         USE hamil_struct_def
         USE pawm
         USE pead, ONLY : LUSEPEAD,LPEAD_NO_SCF,W_STORE,PEAD_EIGENVALUES,PEAD_ACC_ADD_CPROJ,PEAD_ACC_ADD_PW
         USE subrot, ONLY : SHIFT_BANDS_BETWEEN_LOW_HIGH, SHIFT_BACK_BANDS_BETWEEN_LOW_HIGH
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavedes) WGW
         TYPE(wavespin) W
         TYPE(type_info) T_INFO
         TYPE(potcar) P(T_INFO%NTYP)
         TYPE(latt) LATT_CUR
         TYPE(latt) LATT_INI
         TYPE (info_struct) INFO
         INTEGER LMAXMP2
         REAL(q) :: ENCUTGW, ENCUTGWSOFT,change,denom
         TYPE (in_struct) IO
         TYPE (kpoints_struct) KPOINTS
         TYPE (ham_handle)  HAMILTONIAN
         TYPE (symmetry)    SYMM
         TYPE (grid_3d)     GRID       ! grid for wavefunctions
         TYPE (nonlr_struct) NONLR_S
         TYPE (nonl_struct) NONL_S
         INTEGER  LMDIM
         OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
         OVERLAP  CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
         RGRID       SV(DIMREAL(WDES%GRID%MPLWV),WDES%NCDIJ)
         TYPE (energy)      E
         TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
         TYPE (grid_3d)     GRID_SOFT  ! grid for soft chargedensity
         TYPE (grid_3d)     GRIDUS     ! temporary grid in us.F
         TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
         TYPE (transit)     SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
         COMPLEX(q) CHTOT(GRIDC%MPLWV,WDES%NCDIJ) ! charge density
         COMPLEX(q) CHTOTL(GRIDC%MPLWV,WDES%NCDIJ)! old charge-density
         RGRID      DENCOR(GRIDC%RL%NP)
         COMPLEX(q)  CVTOT(GRIDC%MPLWV,WDES%NCDIJ) ! local potential
         COMPLEX(q)  CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
         INTEGER     IRDMAX
         OVERLAP  CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
         INTEGER N_MIX_PAW
         REAL(q)  RHOLM(N_MIX_PAW,WDES%NCDIJ),RHOLM_LAST(N_MIX_PAW,WDES%NCDIJ)
         COMPLEX(q)  CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)

         !local variables
         INTEGER :: NI,NJ,NR,NS,KI,KJ,KR,KS,KQ,KQ_,direction
         integer :: iteration,TWOE4ORBITAL_COLS,TWOE4ORBITAL_ROWS,niteration,readerr
         character(255) :: neci_inputfile
         LOGICAL :: LFOCKMCAR
         LOGICAL :: LFOCK
         REAL(q) :: EXX

! Initialisation of some important variables
         NBANDSDUMP=WDES%NB_TOT
         CALL CC4S_INCAR_READER(IO%IU5, IO%IU6, IO%IU0)
         NBANDSDUMP=NBANDSDUMP-NFREEZE

         CALL CHECK_FULL_KPOINTS ! all set up properly ?
         CALL CHECK_SHIFTED_KPOINTS(WDES,W,KPOINTS)
         IF ((IO%IU0>0) .and. (SHIFTED_KPOINTS)) write(IO%IU0,*)'You are using a shifted k-mesh.'

! Compute Coulomb Vertex

         VBMAX=INFO%NELECT/2
         IF (WDES%ISPIN==2) VBMAX=INFO%NELECT
         !Initialize the 1D process grid
         CALL INIT_BLACS_1D()

         !Initialize the 2D process grid
         CALL INIT_BLACS_2D(WDES)

         CALL SORT_ORBITALINDICES(W,WDES)

         IF (IO%IU0>=0) write(*,*)''




!***********************************************************************
! Compute Fock Matrix. At the moment we assumae canonical orbitals.
! Non-canonical orbitals are not yet supported and LSKIPFOCK=FALSE.
!***********************************************************************
         IF (LSKIPFOCK) THEN
            IF (IO%IU0>=0) write(*,*)'* Assuming canonical HF orbitals. Fock matrix will not be recalculated.'
            ALLOCATE(FOCKM(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN))
            FOCKM=zero
            DO NI=1,WDES%NB_TOT
               FOCKM(NI,NI,:,:)=W%CELTOT(NI,:,:)
            ENDDO
         ELSE
            IF (IO%IU0>=0) write(*,*)'* Computing Fock matrix.'
            LFOCK=.TRUE.
            ! just make sure that data distribution is over bands
            CALL REDIS_PW_OVER_BANDS(WDES, W)
            CALL ON_SYMMETRY
            CALL POTENTIAL_AND_CHARGE(1)
            CALL CALC_FOCKMATRIX(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
                               LMDIM,CDIJ,CQIJ, 1,SV,T_INFO,P,IO%IU0,E%EXHF)
         ENDIF

!         IF (IO%IU0>=0) write(*,*)'- Writing Fock Matrix.'
!         CALL WRITEFOCKMATRIX_YAML(W,WDES)


!***********************************************************************
! Calculate Coulomb Vertex
!***********************************************************************
         CALL CALC_CVERTEX(WDES,WGW,W,P,T_INFO,LATT_CUR,LMAXMP2,ENCUTGW,ENCUTGWSOFT,IO,FSG_STORE(1))
!***********************************************************************
! Write G-mesh to disk
!***********************************************************************
         CALL WRITE_GVEC_POTFAK_YAML(WDES,LATT_CUR)
 

         IF (IO%IU0>=0) write(*,*)'* Computing Delta-kernel two-electron integrals.'
!***********************************************************************
! Delta-kernel two-electron integrals are computed in the PAW framework
!***********************************************************************
         CALL CALC_D2PAW(LATT_CUR, W, WDES, T_INFO, P, IO)
         IF (IO%IU0>=0) write(*,*)'- Writing Delta-kernel two-electron integrals.'
         CALL WRITE_D2PAW_YAML(W, WDES) 

         IF (IO%IU0>=0) write(*,*)'* Computing optimized auxilliary field.'
         IF (IO%IU0>=0) write(IO%IU6,*)'* Computing optimized auxilliary field.'
         CALL CALC_OPTCVERTEX(W,WDES,INFO,IO)
         IF (IO%IU0>=0) write(*,*)'* Estimating truncation of optimized auxilliary field.'
         IF (IO%IU0>=0) write(IO%IU6,*)'Estimating truncation of optimized auxilliary field.'
         CALL ESTIMATE_ACC_NOPTAUX(IO,W,WDES,INFO) 
         IF (IO%IU0>=0) write(*,*)'- Writing CoulombVertexSingularVectors.'

!***********************************************************************
! Write singular vectors for the optimized-auxilliary-field to plane-wave transformation to disk
!***********************************************************************
         CALL WRITE_TSVEC_1D_GVEC_YAML(WDES)

! A sanity test of the computed Coulomb Vertex in the optimized auxilliary field
!         CALL TEST_OPTCVERTEX(W,WDES,INFO,IO)

         IF (LLFTODDUMP) THEN
            CALL WRITEBINCVERTEX2DISK_FTODDUMP(W,WDES)
         ENDIF


!***********************************************************************
! Write Coulomb Vertex to disk
!***********************************************************************
         CALL WRITEBINCVERTEX2DISK_YAML(W,WDES)

!***********************************************************************
! Exit BLACS cleanly
!***********************************************************************
         call BLACS_GRIDEXIT(CONTXT_1D)
         call BLACS_GRIDEXIT(CONTXT_2D)

      CONTAINS

!***********************************************************************
! Several routines needed for calculation of Fock matrix in case non-canonical references shall be used.
! These routines have been copied from different VASP modules. This part needs to be cleaned at some point.
!***********************************************************************
         SUBROUTINE ON_SYMMETRY
            USE pead, ONLY : LPEAD_SYM_RED
            LOGICAL :: LGAMMA

            IF (W%WDES%NKPTS == 1 .AND. SUM(W%WDES%VKPT(:,1)*W%WDES%VKPT(:,1))<G2ZERO)  THEN
               LGAMMA = .TRUE.
            ELSE
               LGAMMA = .FALSE.
            ENDIF

            IF (SYMM%ISYM>=0.AND. .NOT. LGAMMA) THEN

            ENDIF
         END SUBROUTINE ON_SYMMETRY

         SUBROUTINE CALC_FOCKMATRIX(HAMILTONIAN, &
              GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
              LMDIM,CDIJ,CQIJ,IFLAG,SV,T_INFO,P,IU0,EXHF, & 
              CHAMHF,LFIRST,LLAST,NKSTART,NKSTOP,NBANDS_MAX,EXHF_ACFDT)
#ifdef _OPENACC
           USE mopenacc
#endif
           USE prec
           USE wave_high
           USE lattice
           USE mpimy
           USE mgrid
           USE nonl_high
           USE hamil_struct_def
           USE constant
           USE jacobi
           USE scala
           USE main_mpi
           USE fock
           USE pseudo
           USE ini
           USE sym_grad
           USE fileio
           USE openmp, ONLY : omp_nthreads
#ifdef fock_dblbuf
           USE fock_dbl
#endif

           IMPLICIT NONE
           TYPE (ham_handle)  HAMILTONIAN
           TYPE (grid_3d)     GRID
           TYPE (latt)        LATT_CUR
           TYPE (nonlr_struct) NONLR_S
           TYPE (nonl_struct) NONL_S
           TYPE (wavespin)    W
           TYPE (wavedes)     WDES
           TYPE (symmetry) :: SYMM
           INTEGER LMDIM
           OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ),CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
           INTEGER            IFLAG            ! determines mode of diagonalisation
           RGRID   SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ) ! local potential
           TYPE (type_info)   T_INFO
           TYPE (potcar)      P(T_INFO%NTYP)
           INTEGER IU0
           REAL(q) EXHF
! if CHAMHF is present, CHAMHF is added to the subspace matrix at each iteration
! in addition calculation of Fock part is bypassed
           GDEF, OPTIONAL :: CHAMHF(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)
! if LFIRST is supplied CHAMHF is set to CONJG(CHAMHF)-CHAM
! if LLAST  is supplied CHAMHF is set to CONJG(CHAM)
           LOGICAL, OPTIONAL :: LFIRST, LLAST
           INTEGER, OPTIONAL :: NKSTART, NKSTOP     ! start k-point
           INTEGER, OPTIONAL :: NBANDS_MAX  ! maximum band index
           LOGICAL :: LSCAAWARE_LOCAL
           REAL(q), OPTIONAL :: EXHF_ACFDT
    
! local
    ! work arrays for ZHEEV (blocksize times number of bands)
           INTEGER, PARAMETER :: LWORK=32
           GDEF       CWRK(LWORK*WDES%NB_TOT)
           REAL(q)    R(WDES%NB_TOT)
#ifndef USE_ZHEEVX
           REAL(q)    RWORK(3*WDES%NB_TOT)
#else
           REAL(q)    RWORK(7*WDES%NB_TOT), ABSTOL, VL, VU
           INTEGER IWORK(5*WDES%NB_TOT), INFO(WDES%NB_TOT),IL, IU, NB_CALC
#endif
    ! work arrays (do max of 16 strips simultaneously)

           TYPE (wavedes1)    WDES1          ! descriptor for one k-point
           TYPE (wavefun1)    W1             ! current wavefunction
           TYPE (wavefuna)    WA             ! array to store wavefunction
           TYPE (wavespin)    WFOCK          ! array to store the Fock (exchange contribution)
           TYPE (wavefuna)    WNONL          ! array to hold non local part D * wave function character
           TYPE (wavefuna)    WOVL           ! array to hold non local part Q * wave function character
           TYPE (wavefuna)    WHAM           ! array to store accelerations for a selected block
           COMPLEX(q)         CDCHF          ! HF double counting energy
           INTEGER :: ICALL=0                ! number of calls
           GDEF,ALLOCATABLE,TARGET::  CHAM(:,:),COVL(:,:) ! Hamiltonian and overlap matrix
           TYPE (wavefun1)    WTMP(NSTRIP_STANDARD)

    ! array for asyncronous data redistribution
           TYPE (REDIS_PW_CTR),POINTER :: H_PW1, H_PW2

           LOGICAL, EXTERNAL :: USEFOCK_CONTRIBUTION

           INTEGER :: NB_TOT, NBANDS, NSTRIP, ISP, NK, N, NP, NPOS, NSTRIP_ACT, &
                NPOS_RED, NSTRIP_RED, IFAIL, MY_NKSTART, MY_NKSTOP, NSIM_LOCAL
!$  INTEGER NSIM_FOCK

           PROFILING_START('eddiag')

           LscaAWARE_LOCAL=__IF_ACC_OFF__(LSCAAWARE.AND.(.NOT.(IFLAG==4)))
           IF (PRESENT(CHAMHF) .OR. IFLAG==1 .OR. IFLAG==23) LSCAAWARE_LOCAL=.FALSE.

           IF (.NOT. LscaAWARE_LOCAL) THEN
              NB_TOT=WDES%NB_TOT
              ALLOCATE(FOCKM(NB_TOT,NB_TOT,WDES%NKPTS,WDES%ISPIN))
           ELSE
              CALL INIT_scala(WDES%COMM_KIN, NB_TOT)
              ALLOCATE(FOCKM(SCALA_NP(),SCALA_NQ(),WDES%NKPTS,WDES%ISPIN))
           ENDIF   


           CDCHF=0         ! double counting HF
           IF (PRESENT(EXHF_ACFDT)) EXHF_ACFDT=0
           ICALL=ICALL+1

           NB_TOT=WDES%NB_TOT
           NBANDS=WDES%NBANDS
           IF (PRESENT(NBANDS_MAX) .AND. IFLAG == 0) NBANDS=NBANDS_MAX

           NSTRIP=NSTRIP_STANDARD
!$  NSTRIP=NSTRIP*omp_nthreads
 
    ! allocate work space
           IF (.NOT. LscaAWARE_LOCAL) THEN
              ALLOCATE(CHAM(NB_TOT,NB_TOT))
!$ACC ENTER DATA CREATE(CHAM) IF(ACC_EXEC_ON)
!$ACC ENTER DATA COPYIN(CHAMHF) ASYNC(ACC_ASYNC_Q) IF(ACC_EXEC_ON.AND.PRESENT(CHAMHF))
           ELSE
              CALL INIT_scala(WDES%COMM_KIN, NB_TOT)
              ALLOCATE(CHAM(SCALA_NP(),SCALA_NQ()))
           ENDIF
!$ACC ENTER DATA CREATE(WDES1) IF(ACC_EXEC_ON)
           CALL SETWDES(WDES,WDES1,0)
!$ACC ENTER DATA CREATE(WHAM,WNONL) IF(ACC_EXEC_ON)
           CALL NEWWAVA(WHAM, WDES1, NSTRIP)
           CALL NEWWAVA_PROJ(WNONL, WDES1)
!$ACC ENTER DATA CREATE(W1) IF(ACC_EXEC_ON)
           CALL NEWWAV_R(W1, WDES1)

           IF (IFLAG==5) THEN
              ALLOCATE(COVL(NB_TOT,NB_TOT))
!$ACC ENTER DATA CREATE(COVL,WOVL) IF(ACC_EXEC_ON)
              CALL NEWWAVA_PROJ(WOVL, WDES1)
           ENDIF
#ifdef _OPENACC
           CALL ACC_COPYIN_TYPED_VAR(GRID)
           CALL ACC_COPYIN_TYPED_VAR(NONL_S)
           CALL ACC_COPYIN_TYPED_VAR(NONLR_S)
!$ACC ENTER DATA COPYIN(SV,CDIJ,CQIJ) IF(ACC_EXEC_ON)
#endif
           IF (PRESENT(NKSTART)) THEN
              MY_NKSTART=NKSTART
           ELSE
              MY_NKSTART=1
           ENDIF
           IF (PRESENT(NKSTOP)) THEN
              MY_NKSTOP=NKSTOP
           ELSE
              MY_NKSTOP=WDES%NKPTS
           ENDIF

!=======================================================================
! start with HF part and store results WFOCK
!=======================================================================
           IF (USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF))) THEN

              CALL ALLOCW(WDES,WFOCK)
#ifdef _OPENACC
              CALL ACC_COPYIN_TYPED_VAR(WFOCK)
#endif
#ifndef fock_dblbuf
              NSIM_LOCAL=(W%WDES%NSIM*2+W%WDES%NB_PAR-1)/W%WDES%NB_PAR
              IF (LPEAD_NO_SCF()) THEN
                 DO N=1,NSIM_LOCAL
                    CALL NEWWAV(WTMP(N) , WDES1, .FALSE.)
                 ENDDO
              ENDIF

              DO ISP=1,WDES%ISPIN
              DO NK=MY_NKSTART,MY_NKSTOP
#ifdef MPI
                 IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
                 CALL SETWDES(WDES,WDES1,NK)

                 IF (NONLR_S%LREAL) THEN
                    CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
                 ELSE
                    CALL PHASE(WDES,NONL_S,NK)
                 ENDIF

#ifndef _OPENMP
                 DO NPOS=1,NBANDS,NSIM_LOCAL
                    NSTRIP_ACT=MIN(NBANDS+1-NPOS,NSIM_LOCAL)
                    IF (LPEAD_NO_SCF()) THEN
             ! this is a little bit tricky/dirty
             ! if local field effects are not taken into account in the pead,
             ! (i.e. Hamiltonian constructed from initial wavefunctions)
             ! W_STORE (original wavefunctions) must be passed to FOCK_ACC and WTMP 
             ! (wavefunction on which the action is calculated)
             ! must be supplied as an additional argument
                       DO N=NPOS,NPOS+NSTRIP_ACT-1
                          NP=N-NPOS+1
                          CALL W1_COPY(ELEMENT(W,WDES1,N,ISP),WTMP(NP))
                       ENDDO
                       CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W_STORE,   &
                            NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
                            WFOCK%CW(:,NPOS:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF, WTMP, LSYMGRAD=LSYMGRAD)
                    ELSE
                       IF (PRESENT(EXHF_ACFDT)) THEN
                          CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W,  &
                               NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
                               WFOCK%CW(:,NPOS:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF,  LSYMGRAD=LSYMGRAD, EXHF_ACFDT=EXHF_ACFDT)
                       ELSE
                          CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W,  &
                               NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
                               WFOCK%CW(:,NPOS:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF,  LSYMGRAD=LSYMGRAD)
                       ENDIF
                    ENDIF
                 ENDDO
#else
                 IF (LPEAD_NO_SCF()) THEN
                    DO NPOS=1,NBANDS,NSIM_LOCAL
                       NSTRIP_ACT=MIN(NBANDS+1-NPOS,NSIM_LOCAL)
                ! this is a little bit tricky/dirty
                ! if local field effects are not taken into account in the pead,
                ! (i.e. Hamiltonian constructed from initial wavefunctions)
                ! W_STORE (original wavefunctions) must be passed to FOCK_ACC and WTMP 
                ! (wavefunction on which the action is calculated)
                ! must be supplied as an additional argument
                       DO N=NPOS,NPOS+NSTRIP_ACT-1
                          NP=N-NPOS+1
                          CALL W1_COPY(ELEMENT(W,WDES1,N,ISP),WTMP(NP))
                       ENDDO
                       CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W_STORE,   &
                            NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
                            WFOCK%CW(:,NPOS:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF, WTMP, &
                            LSYMGRAD=LSYMGRAD)
                    ENDDO
                 ELSE
                    NSIM_FOCK=NBLOCK_FOCK
                    DO NPOS=1,MIN(W%WDES%NB_TOTK(NK,ISP),NBANDS*W%WDES%NB_PAR),NSIM_FOCK
                       NSTRIP_ACT=MIN(NSIM_FOCK,W%WDES%NB_TOTK(NK,ISP)-NPOS+1)
                       IF (PRESENT(EXHF_ACFDT)) THEN
                          CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W, &
                               NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
                               WFOCK%CW(:,:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF, &
                               LSYMGRAD=LSYMGRAD, EXHF_ACFDT=EXHF_ACFDT, &
                               LCOMMENSURATE=.FALSE.)
                       ELSE
                          CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W, &
                               NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
                               WFOCK%CW(:,:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF, &
                               LSYMGRAD=LSYMGRAD, LCOMMENSURATE=.FALSE.)
                       ENDIF
                    ENDDO
                 ENDIF
#endif
              ENDDO
              ENDDO

              CALLMPI( M_sum_z(WDES%COMM_KIN,CDCHF,1))
              CALLMPI( M_sum_z(WDES%COMM_KINTER,CDCHF,1))
              EXHF=CDCHF
 
              IF (PRESENT(EXHF_ACFDT)) THEN
                 CALLMPI( M_sum_d(WDES%COMM_KIN,EXHF_ACFDT,1))
                 CALLMPI( M_sum_d(WDES%COMM_KINTER,EXHF_ACFDT,1))
              ENDIF
#else
              IF (LPEAD_NO_SCF()) THEN
                 IF (PRESENT(EXHF_ACFDT)) THEN
                    CALL FOCK_ALL_DBLBUF(WDES,W_STORE,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ, &
                   &   EX=EXHF,EX_ACFDT=EXHF_ACFDT,NBMAX=NBANDS*W%WDES%NB_PAR,NKSTART=MY_NKSTART,NKSTOP=MY_NKSTOP, &
                   &   LSYMGRAD=LSYMGRAD,XI=WFOCK,WP=W)
                 ELSE
                    CALL FOCK_ALL_DBLBUF(WDES,W_STORE,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ, &
                   &   EX=EXHF,NBMAX=NBANDS*W%WDES%NB_PAR,NKSTART=MY_NKSTART,NKSTOP=MY_NKSTOP, &
                   &   LSYMGRAD=LSYMGRAD,XI=WFOCK,WP=W)
                 ENDIF
              ELSE
                 IF (PRESENT(EXHF_ACFDT)) THEN
                    CALL FOCK_ALL_DBLBUF(WDES,W,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ, &
                   &   EX=EXHF,EX_ACFDT=EXHF_ACFDT,NBMAX=NBANDS*W%WDES%NB_PAR,NKSTART=MY_NKSTART,NKSTOP=MY_NKSTOP, &
                   &   LSYMGRAD=LSYMGRAD,XI=WFOCK)
                 ELSE
                    CALL FOCK_ALL_DBLBUF(WDES,W,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ, &
                   &   EX=EXHF,NBMAX=NBANDS*W%WDES%NB_PAR,NKSTART=MY_NKSTART,NKSTOP=MY_NKSTOP, &
                   &   LSYMGRAD=LSYMGRAD,XI=WFOCK)
                 ENDIF
              ENDIF
#endif
              IF (LSYMGRAD) &
                   CALL APPLY_SMALL_SPACE_GROUP_OP( W, WFOCK, NONLR_S, NONL_S,P, T_INFO%NIONS, LATT_CUR, SYMM, CQIJ, .FALSE. , -1, MY_NKSTART)
           ENDIF
!=======================================================================
           spin:  DO ISP=1,WDES%ISPIN
           kpoint: DO NK=MY_NKSTART,MY_NKSTOP
#ifdef MPI
              IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
!=======================================================================
              IF (LscaAWARE_LOCAL) CALL INIT_scala(WDES%COMM_KIN, WDES%NB_TOTK(NK,ISP))

              CALL SETWDES(WDES,WDES1,NK)
!=======================================================================
!  IFLAG=0 calculate eigenvalues  only
!=======================================================================
              IF (IFLAG==0) THEN
                 W%CELEN(:,NK,ISP)=0
                 DO N=1,NBANDS
             ! transform wavefunction to real space
             ! and calculate eigenvalues calling ECCP, no redistribution !
                    CALL SETWAV(W, W1, WDES1, N, ISP) ! allocation for W1%CR done above
                    CALL FFTWAV_W1(W1)
!            IF (ASSOCIATED(HAMILTONIAN%AVEC)) THEN
!               CALL ECCP_VEC(WDES1,W1,W1,LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP),HAMILTONIAN%AVEC, W1%CELEN)
                    IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
                       CALL ECCP_TAU(WDES1,W1,W1,LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP),LATT_CUR,HAMILTONIAN%MU(:,ISP),W1%CELEN)    
                    ELSE
                       CALL ECCP(WDES1,W1,W1,LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP), W1%CELEN)
                    ENDIF
                    IF (USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF))) THEN
                       W1%CELEN=W1%CELEN+W1_DOT( ELEMENT( W, WDES1, N, ISP), ELEMENT (WFOCK, WDES1, N, ISP))
                    ENDIF

                    W%CELEN(N,NK,ISP)=W1%CELEN
                 ENDDO
                 CALL PEAD_EIGENVALUES(W,NK,ISP)
                 CYCLE kpoint
              ENDIF

              WA=ELEMENTS(W, WDES1, ISP)
!=======================================================================
!  IFLAG /= 0 calculate Hamiltonian CHAM
!=======================================================================
       !  caclulate D |cfin_n> (D = non local strength of PP)
              IF (WDES%DO_REDIS .AND. LASYNC) THEN
                 CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW1)
                 CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW2)
                 DO NPOS=1,NSTRIP
                    CALL REDIS_PW_START(WDES, WA%CW(1,NPOS), NPOS, H_PW1)
                 ENDDO
              ENDIF

              CALL OVERL(WDES1,.TRUE.,LMDIM,CDIJ(1,1,1,ISP),WA%CPROJ(1,1),WNONL%CPROJ(1,1))
              DO N=1,NBANDS
                 CALL PEAD_ACC_ADD_CPROJ(WNONL%CPROJ(:,N),N,NK,ISP)
              ENDDO

              IF (IFLAG==5) THEN
                 CALL OVERL(WDES1,WDES1%LOVERL,LMDIM,CQIJ(1,1,1,ISP),WA%CPROJ(1,1),WOVL%CPROJ(1,1))
              ENDIF

       ! redistribute the wavefunction characters
              CALL REDISTRIBUTE_PROJ(WA)
              CALL REDISTRIBUTE_PROJ(WNONL)
              IF (IFLAG==5) CALL REDISTRIBUTE_PROJ(WOVL)
!$ACC KERNELS PRESENT(CHAM) __IF_ASYNC__
              CHAM=0
!$ACC END KERNELS
              strip: DO NPOS=1,NBANDS,NSTRIP
                 NSTRIP_ACT=MIN(NBANDS+1-NPOS,NSTRIP)

          !  calculate V_{local} |phi> + T | phi >
          !  for a block containing NSTRIP wavefunctions

          ! set Fock contribution
                 IF (USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF))) THEN
#if PGI_BEFORE_19_1
                    CALL __ZCOPY__(SIZE(WFOCK%CW,1)*NSTRIP_ACT,WFOCK%CW(1,NPOS,NK,ISP),1,WHAM%CW(1,1),1)
#else
!$ACC KERNELS PRESENT(WHAM,WFOCK) __IF_ASYNC__
                    WHAM%CW(:,1:NSTRIP_ACT)=WFOCK%CW(:,NPOS:NPOS+NSTRIP_ACT-1,NK,ISP)
!$ACC END KERNELS
#endif
                 ENDIF

                 DO N=NPOS,NPOS+NSTRIP_ACT-1
                    NP=N-NPOS+1
                    CALL SETWAV(W, W1, WDES1, N, ISP)
                    CALL FFTWAV_W1(W1)
                    IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
                       CALL HAMILT_LOCAL_TAU(W1, SV, LATT_CUR, HAMILTONIAN%MU, ISP, WHAM%CW(:,NP), &
                      &   USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF)), IFLAG/=12) 
                    ELSE
                       CALL HAMILT_LOCAL(W1, SV, ISP, WHAM%CW(:,NP), USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF)), IFLAG/=12)
                    ENDIF
                    CALL PEAD_ACC_ADD_PW(WHAM%CW(:,NP),N,NK,ISP)
                    IF (WDES%DO_REDIS.AND. LASYNC) CALL REDIS_PW_START(WDES, WHAM%CW(1,NP), N, H_PW2)
                 ENDDO
          ! redistribute wavefunctions
          ! after this redistributed up to and including 1...NPOS+NSTRIP_ACT
                 IF (WDES%DO_REDIS) THEN
                    IF (LASYNC) THEN
                       DO N=NPOS,NPOS+NSTRIP_ACT-1
                          NP=N-NPOS+1
                          CALL REDIS_PW_STOP (WDES, WA%CW(1,N), N, H_PW1)
                          IF (N+NSTRIP<=NBANDS) &
                               CALL REDIS_PW_START(WDES, WA%CW(1,N+NSTRIP), N+NSTRIP, H_PW1)
                          CALL REDIS_PW_STOP (WDES, WHAM%CW(1,NP), N, H_PW2)
                       ENDDO
                    ELSE
                       CALL REDISTRIBUTE_PW( ELEMENTS( WA, NPOS, NPOS-1+NSTRIP_ACT))
                       CALL REDISTRIBUTE_PW( ELEMENTS( WHAM, 1, NSTRIP_ACT))
                    ENDIF
                 ENDIF

                 NPOS_RED  =(NPOS-1)*WDES%NB_PAR+1
                 NSTRIP_RED=NSTRIP_ACT*WDES%NB_PAR

                 IF (.NOT. LscaAWARE_LOCAL) THEN
                    CALL ORTH1('U', &
                      WA%CW_RED(1,1),WHAM%CW(1,1),WA%CPROJ_RED(1,1), &
                      WNONL%CPROJ_RED(1,NPOS_RED),NB_TOT, &
                      NPOS_RED, NSTRIP_RED, WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CHAM(1,1))
                 ELSE
                    CALL ORTH1_DISTRI('U', &
                      WA%CW_RED(1,1),WHAM%CW(1,1),WA%CPROJ_RED(1,1), &
                      WNONL%CPROJ_RED(1,NPOS_RED),NB_TOT, &
                      NPOS_RED, NSTRIP_RED, WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CHAM(1,1), & 
                      WDES%COMM_KIN, WDES%NB_TOTK(NK,ISP))
                 ENDIF
              ENDDO strip

              IF (WDES%DO_REDIS .AND. LASYNC) THEN
                 CALL REDIS_PW_DEALLOC(H_PW1)
                 CALL REDIS_PW_DEALLOC(H_PW2)
              ENDIF

              IF (.NOT. LscaAWARE_LOCAL) THEN
                 CALLMPI( M_sum_g(WDES%COMM_KIN,CHAM(1,1),NB_TOT*NB_TOT))
                 ! add lower triangle
!$ACC PARALLEL LOOP GANG PRESENT(CHAM) __IF_ASYNC__
                 DO N=1,NB_TOT
!$ACC LOOP VECTOR
                    DO NP=N+1,NB_TOT
                       CHAM(NP,N)=GCONJG(CHAM(N,NP))
                    ENDDO
                 ENDDO
                 IF (IFLAG==23) CALL RESTRICT_TO_OCCUPIED_ONLY(WDES%NB_TOTK(NK,ISP), CHAM, W%FERTOT(:,NK, ISP))
              ENDIF

              IF (PRESENT(CHAMHF).AND.PRESENT(LFIRST)) THEN
!$ACC KERNELS PRESENT(CHAMHF,CHAM) __IF_ASYNC__
                 IF (LFIRST) THEN
                    CHAMHF(:,:,NK,ISP)=GCONJG(CHAMHF(:,:,NK,ISP))-CHAM(:,:)
                 ENDIF
                 CHAM(:,:)=CHAM(:,:)+CHAMHF(:,:,NK,ISP)
!$ACC END KERNELS
              ENDIF
#ifdef debug
!$ACC UPDATE SELF(CHAM) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
              IF (IU0>=0) CALL DUMP_HAM( "Hamilton matrix subrot",WDES, CHAM)
#endif
!-----------------------------------------------------------------------
! calculate the overlap matrix
!-----------------------------------------------------------------------
              IF (IFLAG==5) THEN
!$ACC KERNELS PRESENT(COVL) __IF_ASYNC__
                 COVL=(0._q,0._q)
!$ACC END KERNELS
                 DO NPOS=1,NB_TOT-NSTRIP_STANDARD_GLOBAL,NSTRIP_STANDARD_GLOBAL
                    CALL ORTH1('U',WA%CW_RED(1,1),WA%CW_RED(1,NPOS),WA%CPROJ_RED(1,1), &
                         WOVL%CPROJ_RED(1,NPOS),NB_TOT, &
                         NPOS,NSTRIP_STANDARD_GLOBAL,WDES1%NPL_RED,WDES1%NPRO_O_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,COVL(1,1))
                 ENDDO

                 CALL ORTH1('U',WA%CW_RED(1,1),WA%CW_RED(1,NPOS),WA%CPROJ_RED(1,1), &
                      WOVL%CPROJ_RED(1,NPOS),NB_TOT, &
                      NPOS,NB_TOT-NPOS+1,WDES1%NPL_RED,WDES1%NPRO_O_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,COVL(1,1))
                 CALLMPI( M_sum_g(WDES%COMM_KIN,COVL(1,1),NB_TOT*NB_TOT))
#ifdef debug
!$ACC UPDATE SELF(COVL) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
                 IF (IU0>=0) CALL DUMP_HAM( "Overlap matrix",WDES, COVL)
#endif
              ENDIF
!=======================================================================
! IFLAG =2
!  simply copy eigenvalues
!=======================================================================
              IF (IFLAG==12) THEN
!$ACC UPDATE SELF(CHAM) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
                 IF (IU0>=0) CALL DUMP_HAM( "Hamilton matrix",WDES, CHAM)
                 CALL vtutor%stopCode()
              ENDIF

              IF (.NOT. LscaAWARE_LOCAL) THEN
!$ACC PARALLEL LOOP PRESENT(WDES,W,CHAM) __IF_ASYNC__
                 DO N=1,WDES%NB_TOTK(NK,ISP)
                    W%CELTOT(N,NK,ISP)=CHAM(N,N)
                 ENDDO
!$ACC UPDATE SELF(W%CELTOT(1:WDES%NB_TOTK(NK,ISP),NK,ISP)) __IF_ASYNC__
              ENDIF
!=======================================================================
! IFLAG =4 use Loewdin perturbation to get rotation matrix
! this preserves the ordering of the eigenvalues
! MIND: does not work for real matrices
!=======================================================================
              IF (IFLAG==4) THEN
                 CALL LOEWDIN_DIAG(WDES%NB_TOTK(NK,ISP), NB_TOT, CHAM)

                 CALL ORSP(WDES%NB_TOTK(NK,ISP), NB_TOT, NB_TOT, CHAM)
!test writegamma
!          IF (IU0>=0) CALL WRITEGAMMA(NK, WDES%NB_TOTK(NK,ISP), NB_TOT, CHAM, .TRUE.)
#ifdef debug
                  IF (IU0>=0) CALL DUMP_HAM( "Loewdin",WDES, CHAM)
#endif
              ELSE
!=======================================================================
! IFLAG > 1 and IFLAG <4
! diagonalization of CHAM
! we have lots of choices for the parallel version
! this  makes things rather complicated
! to allow for reasonable simple programming, once the diagonalisation
! has been done I jump to line 100
!=======================================================================
                 IF (IFLAG==1) GOTO 1000

#ifndef gammareal
                 IF (.NOT. LscaAWARE_LOCAL) THEN
NOACC !$OMP PARALLEL DO DEFAULT(SHARED) &
NOACC !$OMP PRIVATE(N)
!$ACC KERNELS PRESENT(WDES,CHAM) __IF_ASYNC__
                    DO N=1,WDES%NB_TOTK(NK,ISP)
                       IF (ABS(AIMAG(CHAM(N,N)))>1E-2_q .AND. IU0>=0) THEN
NOACC !$OMP CRITICAL (omp_wrt_stdout)
NOACC                     WRITE(IU0,'(A,I5,E14.3)')'WARNING in EDDIAG: sub space matrix is not hermitian',N,AIMAG(CHAM(N,N))
NOACC !$OMP END CRITICAL (omp_wrt_stdout)
                       ENDIF
                       CHAM(N,N)= REAL( CHAM(N,N) ,KIND=q)
                    ENDDO
!$ACC END KERNELS
NOACC !$OMP END PARALLEL DO
                 ELSE
                    CALL BG_CHANGE_DIAGONALE(WDES%NB_TOTK(NK,ISP),CHAM(1,1),IU0)
                 ENDIF
#endif

          !
          ! parallel versions
          ! if fast Jacobi method exists use it (T3D, T3E only)
          ! use the first line in that case
                 CALL SHIFT_BANDS_BETWEEN_LOW_HIGH( WDES, WDES%NB_TOTK(NK,ISP), LscaAWARE_LOCAL, CHAM )

                 IFAIL=0

#if defined(MPI)
                 IF ( LJACOBI .AND. IFLAG /=13 ) THEN
                    IF (IU0>=0) WRITE(IU0,*)'jacobi called'
                    CALL jacDSSYEV(WDES%COMM_KIN, CHAM(1,1), R, NB_TOT)
                    CALLMPI( M_sum_g(WDES%COMM_KIN, CHAM(1,1),NB_TOT*NB_TOT))
                    CALLMPI( M_sum_g(WDES%COMM_KIN, R , NB_TOT))

                    GOTO 100
                 ENDIF

          ! use scaLAPACK if available in parallel version
                 IF (__IF_ACC_OFF__(LscaLAPACK.AND.IFLAG/=5)) THEN
                    IF (.NOT. LscaAWARE_LOCAL) THEN
                       CALL pDSSYEX_ZHEEVX(WDES%COMM_KIN, CHAM(1,1), R,  NB_TOT, WDES%NB_TOTK(NK,ISP))
                       CALLMPI( M_sum_g(WDES%COMM_KIN, CHAM(1,1),NB_TOT*NB_TOT))
                    ELSE
                       CALL BG_pDSSYEX_ZHEEVX(WDES%COMM_KIN, CHAM(1,1), R,  WDES%NB_TOTK(NK,ISP))
                    ENDIF
                    stotmg
                    GOTO 100
                 ENDIF
#endif
          !
          !  seriell codes
          !
#ifdef  gammareal
                 IF (IFLAG == 5) THEN
                    CALL __DSYGV__ &
                         (1,'V','U',WDES%NB_TOTK(NK,ISP),CHAM(1,1),NB_TOT,COVL(1,1),NB_TOT, &
                         R,CWRK,LWORK*NB_TOT, IFAIL)
                 ELSE
#ifndef USE_ZHEEVX
                    CALL __DSYEV__ &
                         ('V','U',WDES%NB_TOTK(NK,ISP),CHAM(1,1),NB_TOT, &
                         R,CWRK,LWORK*NB_TOT, IFAIL)
#else
                    ABSTOL=1E-10_q
                    VL=0 ; VU=0 ; IL=0 ; IU=0
                    ALLOCATE(COVL(NB_TOT,NB_TOT))
!$ACC ENTER DATA CREATE(COVL) __IF_ASYNC__
                    CALL __DSYEVX__ &
                         ( 'V', 'A', 'U', WDES%NB_TOTK(NK,ISP), CHAM(1,1) , NB_TOT, VL, VU, IL, IU, &
                         ABSTOL , NB_CALC , R, COVL(1,1), NB_TOT, CWRK, &
                         LWORK*NB_TOT, RWORK, IWORK, INFO, IFAIL)
                    CALL __DCOPY__(NB_TOT*NB_TOT,COVL,1,CHAM,1)
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC EXIT DATA DELETE(COVL) IF(ACC_EXEC_ON)
                    DEALLOCATE(COVL)
#endif
                 ENDIF
#else
                 IF (IFLAG == 5) THEN
                    CALL __ZHEGV__ &
                         (1,'V','U',WDES%NB_TOTK(NK,ISP),CHAM(1,1),NB_TOT,COVL(1,1),NB_TOT, &
                         R,CWRK,LWORK*NB_TOT, RWORK, IFAIL)
                 ELSE
#ifndef USE_ZHEEVX
                    CALL __ZHEEV__ &
                         ('V','U',WDES%NB_TOTK(NK,ISP),CHAM(1,1),NB_TOT, &
                         R,CWRK,LWORK*NB_TOT, RWORK, IFAIL)
#else
                    ABSTOL=1E-10_q
                    VL=0 ; VU=0 ; IL=0 ; IU=0
                    ALLOCATE(COVL(NB_TOT,NB_TOT))
!$ACC ENTER DATA CREATE(COVL) __IF_ASYNC__
                    CALL __ZHEEVX__ &
                         ( 'V', 'A', 'U', WDES%NB_TOTK(NK,ISP), CHAM(1,1) , NB_TOT, VL, VU, IL, IU, &
                         ABSTOL , NB_CALC , R, COVL(1,1), NB_TOT, CWRK, &
                         LWORK*NB_TOT, RWORK, IWORK, INFO, IFAIL)
                    CALL __ZCOPY__(NB_TOT*NB_TOT,COVL,1,CHAM,1)
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC EXIT DATA DELETE(COVL) IF(ACC_EXEC_ON)
                    DEALLOCATE(COVL)
#endif
                 ENDIF
#endif
          ! T3D uses a global sum which does not guarantee to give the same results on all nodes
          ! the following line is required to make the code waterproof (we had problems)
          ! since we now use a propritary sum (see mpi.F) we should not require
          ! this broadcast anymore
          ! stotmg
          ! CALLMPI( M_bcast_g(WDES%COMM, CHAM(1,1), NB_TOT*NB_TOT))

100              CONTINUE
                 stotmg

                 IF (IFAIL/=0) THEN
                    CALL vtutor%error("ERROR in EDDIAG: call to ZHEEV/ZHEEVX/DSYEV/DSYEVX failed! error code &
                       &was " // str(IFAIL))
                 ENDIF

          ! shift eigenvalues back
                 IF (IFLAG==23) CALL SHIFT_UNOCCUPIED_BACK(WDES%NB_TOTK(NK,ISP), R, W%FERTOT(:,NK, ISP))
                 CALL SHIFT_BACK_BANDS_BETWEEN_LOW_HIGH(WDES, WDES%NB_TOTK(NK,ISP), R )

                 DO N=1,WDES%NB_TOTK(NK,ISP)
                    W%CELTOT(N,NK,ISP)=R(N)
                 ENDDO
              ENDIF
!=======================================================================
! IFLAG > 2
! rotate wavefunctions
!=======================================================================
              IF (IFLAG==2) GOTO 1000

              IF (.NOT. LscaAWARE_LOCAL) THEN
                 CALL LINCOM('F',WA%CW_RED(:,:),WA%CPROJ_RED(:,:),CHAM(1,1), &
                   WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), & 
                   WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
                   WA%CW_RED(:,:),WA%CPROJ_RED(:,:))
              ELSE
                 CALL LINCOM_DISTRI('F',WA%CW_RED(1,1),WA%CPROJ_RED(1,1),CHAM(1,1), &
                   WDES%NB_TOTK(NK,ISP), & 
                   WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
                   WDES%COMM_KIN, NBLK )
              ENDIF

1000          CONTINUE
       
              !  back redistribution over bands
              IF (WDES%DO_REDIS) THEN
                 CALL REDISTRIBUTE_PROJ( ELEMENTS( W, WDES1, ISP))
                 IF (LASYNC) THEN
                    W%OVER_BAND=.TRUE.
                 ELSE
                    CALL REDISTRIBUTE_PW( ELEMENTS( W, WDES1, ISP))
                 ENDIF
                 DWRITE "redis ok"
              ENDIF
#ifdef ACC_DEBUG
!$ACC UPDATE SELF(W%CW(:,:,NK,ISP),W%CPROJ(:,:,NK,ISP)) __IF_ASYNC__
#endif
       ! return Hamiltonian (IFLAG==1) or rotation matrix (IFLAG=3,4)
              IF (PRESENT(LLAST)) THEN       
                 IF (LLAST) THEN
!$ACC KERNELS PRESENT(CHAMHF,CHAM) __IF_ASYNC__
                    IF (IFLAG==1) CHAM=GCONJG(CHAM)
                    CHAMHF(:,:,NK,ISP)=CHAM(:,:)
!$ACC END KERNELS
                 ENDIF
              ENDIF

              FOCKM(:,:,NK,ISP)=CHAM(:,:)
           ENDDO kpoint
           ENDDO spin

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)

!    CALLMPI( M_sum_z(WDES%COMM_KIN,CDCHF,1))
!    CALLMPI( M_sum_z(WDES%COMM_KINTER,CDCHF,1))
!    EXHF=CDCHF
!
!    IF (PRESENT(EXHF_ACFDT)) THEN
!       CALLMPI( M_sum_d(WDES%COMM_KIN,EXHF_ACFDT,1))
!       CALLMPI( M_sum_d(WDES%COMM_KINTER,EXHF_ACFDT,1))
!    ENDIF

    ! need to correct CELTOT at this point so that it is correct on all nodes
           IF (IFLAG==0) THEN
               CALL MRG_CEL(WDES,W)
           ENDIF
#ifdef MPI
           IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
              CALL KPAR_SYNC_CELTOT(WDES,W)
           ENDIF

           IF (IFLAG>=3.AND.(LHFCALC.OR.LUSEPEAD()).AND.WDES%COMM_KINTER%NCPU.GT.1) THEN
              CALL KPAR_SYNC_FERTOT(WDES,W)
              CALL KPAR_SYNC_CELTOT(WDES,W)
              CALL KPAR_SYNC_AUXTOT(WDES,W)
              CALL KPAR_SYNC_WAVEFUNCTIONS(WDES,W)
           ENDIF
#endif
!$ACC UPDATE DEVICE(W%CELTOT) __IF_ASYNC__

    ! deallocation ...
!$ACC WAIT IF(ACC_EXEC_ON)
!$ACC EXIT DATA COPYOUT(CHAMHF) IF(ACC_EXEC_ON.AND.PRESENT(CHAMHF)) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(SV,CDIJ,CQIJ) __IF_ASYNC__

!$ACC EXIT DATA DELETE(CHAM) IF((.NOT.LscaAWARE_LOCAL).AND.ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
           DEALLOCATE(CHAM)

           CALL DELWAV_R(W1)
!$ACC EXIT DATA DELETE(W1) __IF_ASYNC__

           CALL DELWAVA(WHAM)
           CALL DELWAVA_PROJ(WNONL)
!$ACC EXIT DATA DELETE(WHAM,WNONL) __IF_ASYNC__

           IF (USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF))) THEN 
#ifdef _OPENACC
              CALL ACC_DELETE_TYPED_VAR(WFOCK)
#endif
              CALL DEALLOCW(WFOCK)
              IF (LPEAD_NO_SCF()) THEN
                 DO N=1,NSIM_LOCAL
                    CALL DELWAV(WTMP(N), .FALSE.)
                 ENDDO
              ENDIF
           ENDIF

           IF (IFLAG==5) THEN
!$ACC EXIT DATA DELETE(COVL) __IF_ASYNC__
              DEALLOCATE(COVL)
              CALL DELWAVA_PROJ(WOVL)
!$ACC EXIT DATA DELETE(WOVL) __IF_ASYNC__
           ENDIF

           IF (PRESENT(CHAMHF).AND.PRESENT(LFIRST)) THEN
              LFIRST=.FALSE.
           ENDIF
#ifdef _OPENACC
           CALL ACC_DELETE_TYPED_VAR(NONLR_S)
           CALL ACC_DELETE_TYPED_VAR(NONL_S)
           CALL ACC_DELETE_TYPED_VAR(GRID)
           CALL ACC_DELETE_TYPED_VAR(WDES1)
!$ACC WAIT IF(ACC_EXEC_ON)
#endif
           PROFILING_STOP('eddiag')

           RETURN

  
         END SUBROUTINE CALC_FOCKMATRIX


         SUBROUTINE RESTRICT_TO_OCCUPIED_ONLY(NB_TOT, CHAM, FERTOT)
#ifdef _OPENACC
           USE mopenacc_struct_def
#endif
           INTEGER :: NB_TOT
           GDEF    :: CHAM(:,:)
           REAL(q) :: FERTOT(:)
           INTEGER :: NB_OCC, I, J

           ! seek first occupancy that differs from 1.00
           DO NB_OCC=1, NB_TOT
              IF (ABS((FERTOT(NB_OCC))-1.0_q)>1E-5_q) EXIT
           ENDDO

           ! now set the lower (upper) triangle starting with row NB_OCC to zero
           ! loop over row index
!$ACC PARALLEL LOOP GANG PRESENT(CHAM) __IF_ASYNC__
           DO I=NB_OCC, NB_TOT
              ! loop over column index
!$ACC LOOP VECTOR
              DO J=1,I-1
                 CHAM(I,J)=0
                 CHAM(J,I)=0
              ENDDO
              ! shift diagonal elements by 10 eV
              ! so that the diagonalization does not mix occupied and unoccupied manyfold
              ! this shift needs to be removed by  SHIFT_UNOCCUPIED_BACK
              CHAM(I,I)=CHAM(I,I)+10.0_q
           ENDDO

         END SUBROUTINE RESTRICT_TO_OCCUPIED_ONLY

         SUBROUTINE SHIFT_UNOCCUPIED_BACK(NB_TOT, R, FERTOT)
           INTEGER :: NB_TOT
           REAL(q) :: R(:)
           REAL(q) :: FERTOT(:)
           INTEGER :: NB_OCC, I

           ! seek first occupancy that differs from 1.00
           DO NB_OCC=1, NB_TOT
              IF (ABS((FERTOT(NB_OCC)-1.0_q))>1E-5_q) EXIT
           ENDDO

           ! now set the lower (upper) triangle starting with row NB_OCC to zero
           ! loop over row index
           DO I=NB_OCC, NB_TOT
              R(I)=R(I)-10.0_q
           ENDDO
         END SUBROUTINE SHIFT_UNOCCUPIED_BACK


!***********************************************************************
!
! update potential and charge
!
!***********************************************************************

         SUBROUTINE POTENTIAL_AND_CHARGE(NELM_HF)
           IMPLICIT NONE
           INTEGER NELM_HF
           INTEGER :: IRDMAA             ! temporary
           REAL(q) :: XCSIF(3,3)         ! temporary (stress tensor)

    ! update charge
           CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
                GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
                LATT_CUR, P, SYMM, T_INFO, &
                CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

    !  calculate potential (Hartree + ionic only)
           CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES,  &
                INFO,P,T_INFO,E,LATT_CUR, &
                CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF)

    ! add the one center augmentation related terms  
           CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
                LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)
    ! finally add one center terms
           CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
                WDES%NCDIJ, LMDIM, CDIJ(1,1,1,1),  RHOLM, CRHODE(1,1,1,1), &
                E,  LMETA=.FALSE., LASPH=INFO%LASPH, LCOREL= .FALSE.  )
         
         END SUBROUTINE POTENTIAL_AND_CHARGE




      END SUBROUTINE CC4S_INTERFACE

!*********************************************************************************
!
! The <p|G|q> quantities are most important for the construction of the two electron
! four orbital integrals <pr|qs>.  <pr|qs>=4pi e^2 \sum_G <p|-G|q> <r|G|s> /(G+q)^2
!
!*********************************************************************************

      SUBROUTINE CALC_CVERTEX(WDES,WGW,W,P,T_INFO,LATT_CUR,LMAXMP2,ENCUTGW,ENCUTGWSOFT,IO,FSG)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE mkpoints
         USE lattice
         USE constant
         USE base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavedes) WGW
         TYPE(wavespin) W
         TYPE(type_info) T_INFO
         TYPE(potcar) P(T_INFO%NTYP)
         TYPE(latt) LATT_CUR
         INTEGER LMAXMP2
         TYPE (in_struct) IO
         ! local variables
         TYPE(wavespin) WHF
         TYPE(wavedes1) WGWQ
         TYPE(wavedes1) WDESKI,WDESKA,WDESKB
         TYPE(wavefun1), ALLOCATABLE :: WI(:), WA(:), WB(:)
         INTEGER KQ,KI,KA,KB,KI_IN_FULL_ORIG,KQ_
         INTEGER NSTRIP, NSTRIPA, ISP, NFFT
         INTEGER NBI,NBA,NBAA, i, nbi_start, nbi_end,rnba
         INTEGER NP ! number of plane waves for the overlap density i*(r)a(r) GCHGIA(r)
         REAL(q) :: FSG ! singularity correction
         REAL(q) :: ENCUTGW,ENCUTGWSOFT
         REAL(q) :: POTFAK(GRIDHF%MPLWV)
         COMPLEX(q) CPHASE(GRIDHF%MPLWV)
         COMPLEX(q) CPHASE2(GRIDHF%MPLWV)
         LOGICAL LPHASE
         LOGICAL LPHASE2
         COMPLEX(q) :: GWORK(MAX(GRIDHF%MPLWV,WGW%GRID%MPLWV)) !work array for calculating GCHGIA
         COMPLEX(q), ALLOCATABLE :: GCHGIA(:,:,:)  ! charge
         GDEF      , ALLOCATABLE :: CRHOIA(:,:)    ! one-center charge
         GDEF      , ALLOCATABLE :: CRHOIB(:,:)    ! one-center charge
         GDEF      , ALLOCATABLE :: CRHOLM(:)      ! augmentation occupancy matrix
         COMPLEX(q), ALLOCATABLE :: TMP_CVERTEX(:,:,:)
         Real(q) :: mem_req

         CALL CHECK_FULL_KPOINTS

         IF (LMAXMP2>=0) THEN
           CALL SET_UP_ONE_CENTER_H(WDES,P,T_INFO,LMAXMP2,H)
         ENDIF
!         FSG=zero
         WHF=W
         WHF%WDES => WDES_FOCK
         NSTRIP=30
         CALL SETWDES(WHF%WDES,WDESKI,0)
         CALL SETWDES(WHF%WDES,WDESKA,0)
         CALL SETWDES(WHF%WDES,WDESKB,0)

         MCALPHA=0
         IF (MCALPHA/=0) THEN
            FSG=0
         ENDIF

         VBMAXDN=LAST_FILLED_XI_NOMOD(W,1,1)
         VBMAXUP=LAST_FILLED_XI_NOMOD(W,1,WDES%ISPIN)
         IF (WDES%ISPIN==1)  VBMAXDN=0
         !WRITE(*,*)'vbmaxdn, vbmaxup ',VBMAXDN,VBMAXUP
         ALLOCATE(WI(NSTRIP),WA(NSTRIP),WB(NSTRIP))
         DO NBI=1,NSTRIP !VBMAX
            CALL NEWWAV(WI(NBI),WDESKI,.TRUE.)
         ENDDO
         DO NBA=1,NSTRIP
            CALL NEWWAV(WA(NBA),WDESKA,.TRUE.)
         ENDDO
         DO NBA=1,NSTRIP
            CALL NEWWAV(WB(NBA),WDESKB,.TRUE.)
         ENDDO

         DO KQ=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0.00_q)) CYCLE
            NGVECTOR=(WGW%NGVECTOR(KQ))
         ENDDO
         NGVECTOR=NGVECTOR

         ALLOCATE(POTFAK_FULL(NGVECTOR,REALNKPTS),GVEC_FULL(3,NGVECTOR,REALNKPTS))
         CALL CALC_GVEC_FULL(WGW,WDES,LATT_CUR)

         IF (ALLOCATED(CVERTEX_2D)) THEN
            DEALLOCATE(CVERTEX_2D)
         ENDIF

         IF (IO%IU0>=0) write(IO%IU6,*)'Allocating memory for  Coulomb-Vertex'
         CALL SETUP_CVERTEX_2D(WDES)

         call BLACS_GRIDINFO(CONTXT_1D, NPROW_1D, NPCOL_1D, MYROW_1D, MYCOL_1D)
         SRCPROC=MOD(NFREEZE,PROCS)
         NBANDSDUMPLOC=numroc(NBANDSDUMP,1,MYCOL_1D,SRCPROC,PROCS)
         NFREEZELOC=numroc(NFREEZE,1,MYCOL_1D,0,PROCS)
         ALLOCATE(TMP_CVERTEX(NGVECTOR,NBANDSDUMPLOC,NSTRIP))

         IF (IO%IU0>0) THEN
            WRITE(IO%IU0,*)'* Computing Coulomb Vertex.'
            WRITE(IO%IU6,*)'* Computing Coulomb Vertex.'
         ENDIF

         spin: DO ISP=1,WDES%ISPIN
         kqloop: DO KQ=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0.00_q)) CYCLE
            CALL SETWDES(WGW,WGWQ,KQ)

            NP=WGWQ%NGVECTOR
            IF (NP>NGVECTOR) THEN
               WRITE(*,*)'Internal error in "CALC_CVERTEX": NP larger than NGVECTOR'
               EXIT
            ENDIF
            ALLOCATE(GCHGIA(NP,NSTRIP,2),CRHOLM(AUG_DES%NPRO*WDES%NRSPINORS))
            kiloop: DO KI=1,WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI)==0)) CYCLE
               tmp_CVERTEX=(0._q,0._q)

               KI_IN_FULL_ORIG=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,KI),KPOINTS_FULL_ORIG)
               ! collect NSTRIP nbi bands at k_i
               CALL SETWDES(WHF%WDES,WDESKI,KI)

               DO NBI_START=1,NBANDSDUMP,NSTRIP
                  call BLACS_GRIDINFO(CONTXT_1D, NPROW_1D, NPCOL_1D, MYROW_1D, MYCOL_1D)
                  NBI_END=NBI_START+MIN(NBANDSDUMP-NBI_start+1,NSTRIP)-1

                  CALL W1_GATHER_GLB(WHF,NBI_START+NFREEZE,NBI_END+NFREEZE,ISP,WI)
                  ! k_b = k_i - k_q - G
                  KB=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KQ)+WDES%VKPT(:,KI),KPOINTS_FULL)
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB)==0)) Write(*,*)'error in calc_cvertex with shifted k-mesh'
                  ! k_a = k_i + k_q - G
                  KA=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KI)-WDES%VKPT(:,KQ),KPOINTS_FULL)
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0)) Write(*,*)'error in calc_cvertex with shifted k-mesh'

                  CALL SETWDES(WHF%WDES,WDESKA,KA)

                  ! CPHASE(r) = e^iGr, where G = k_i - k_q - k_b
                  CALL SETPHASE(WDES%VKPT(:,KI)-WDES%VKPT(:,KQ)-WDES%VKPT(:,KA),GRIDHF,CPHASE,LPHASE)

                  CALL SET_GFAC_WITHOUT_WEIGHT(GRIDHF,LATT_CUR,KI,KA,FSG,POTFAK)
                  ! 1/(G+q)**2

                  IF (ENCUTGW /= ENCUTGWSOFT .AND. ENCUTGW > 0 .AND. ENCUTGWSOFT > 0) THEN
                     CALL SET_GFAC_WAVEFUN(WGWQ,LATT_CUR,FSG,POTFAK,ENCUTGW,ENCUTGWSOFT)
                  ELSE
                     CALL SET_GFAC_WAVEFUN(WGWQ,LATT_CUR,FSG,POTFAK)
                  ENDIF

                  IF (MCALPHA/=0) THEN
!                 setup routine for multipole correction. experimental feature.
                     CALL FOCK_MULTIPOLE_CORR_SETUP(LATT_CUR, GRIDHF)
                     CALL FOCK_MULTIPOLE_SETUP_WGW(WGW, LATT_CUR, GRIDHF)
                  ENDIF

                  POTFAK=SQRT(POTFAK)
                  POTFAK_FULL(1:NGVECTOR,RKQofKQ(KQ))=POTFAK(1:NGVECTOR)*SQRT(1.0_q*WGWQ%GRID%NPLWV)

                  ! loop over all bands
                  DO NBA=1,NBANDSDUMPLOC,NSTRIP
                     NSTRIPA=MIN(NBANDSDUMPLOC+1-NBA,NSTRIP)
                  ! FFT{psi_a} to real space
                  DO NBAA=1,NSTRIPA !copy and fourier transform NSTRIP wave functions
                     CALL W1_COPY( ELEMENT(WHF,WDESKA,NFREEZELOC+NBA+NBAA-1,ISP),WA(NBAA))
                     CALL FFTWAV_W1(WA(NBAA))
                  ENDDO
                  ! loop over valence bands only
                  DO NBI=1,MIN(NBANDSDUMP-NBI_start+1,NSTRIP)
                     GCHGIA=0

                     !loop over all bands in NSTRIP
                     DO NBAA=1,NSTRIPA

                     ! calculate rho(r)=psi_i(r)* psi_a(r) for,
                     ! one center terms and, on the plane wave grid.
                        CALL FOCK_CHARGE_NOINT( WI(NBI),WA(NBAA), GWORK(1), &
                             CRHOLM, SIZE(CRHOLM))
                        ! Set phase e^iGr, where G = k_i - k_q - k_b
                        IF (LPHASE) THEN
                           CALL APPLY_PHASE( GRIDHF, CPHASE(1), GWORK(1), &
                             GWORK(1))
                        ENDIF

                     ! FFT{rho} to reciprocal space
                        CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), &
                          GWORK(1),GCHGIA(1,NBAA,1),WGWQ%GRID,.FALSE.)
                        NFFT=NFFT+1
                     ! multiply with potential factor
                        IF (MCALPHA/=0) THEN
                           CALL APPLY_GFAC_MULTIPOLE_WAVEFUN(WGWQ,GCHGIA(1,NBAA,1), &
                            POTFAK(1))
                        ELSE
                           CALL APPLY_GFAC_WAVEFUN(WGWQ,GCHGIA(1,NBAA,1), &
                            POTFAK(1))
                        ENDIF

                     ENDDO !NBAA (loop over all bands in NSTRIP)


                     DO NBAA=1,NSTRIPA
                        TMP_CVERTEX(1:NP,NBAA+NBA-1,NBI)=(GCHGIA(1:NP,NBAA,1))*SQRT(1.0_q/GRIDHF%NPLWV)

                     ENDDO

                     ENDDO !NBI (loop over nstrip nbi bands only)
                  ENDDO !NBA (loop over all bands)

                  CALL REDISTRIBUTE_CVERTEX_1D(IO,WDES,KI,KQ,ISP,TMP_CVERTEX,NBI_START,NBI_END)
               ENDDO !NBI_start (loop over all bands)
            ENDDO kiloop

            DEALLOCATE(GCHGIA,CRHOLM)
         ENDDO kqloop
         ENDDO spin
         iF (ALLOCATED(CRHOIA)) DEALLOCATE(CRHOIA)
         IF (ALLOCATED(CRHOIB)) DEALLOCATE(CRHOIB)
         IF (ALLOCATED(tmp_CVERTEX)) DEALLOCATE(tmp_CVERTEX)

         DO NBI=1,NSTRIP
            CALL DELWAV(WI(NBI),.TRUE.)
         ENDDO
         DO NBA=1,NSTRIP
            CALL DELWAV(WA(NBA),.TRUE.)
         ENDDO
         DO NBA=1,NSTRIP
            CALL DELWAV(WB(NBA),.TRUE.)
         ENDDO
         DEALLOCATE(WI,WA,WB)

         RETURN
      END SUBROUTINE CALC_CVERTEX

      SUBROUTINE INSULATOR_CHECK(W,WDES)
         IMPLICIT NONE
         type(wavespin) :: W
         TYPE(wavedes) WDES
         INTEGER :: NI,KI,ISP

         DO ISP=1,WDES%ISPIN
         DO NI=1,WDES%NBANDS
            DO KI=1,WDES%NKPTS
               IF ((W%FERTOT(NI,KI,ISP)>0.01_q) .and. (W%FERTOT(NI,KI,ISP)<0.99_q)) THEN
                  WRITE(*,*)"THERE SEEMS TO BE A PROBLEM WITH THE OCCUPATION NUMBERS!!!"
                  WRITE(*,*) "BAND NUMBER:",NI,"IS SOMEWHERE BETWEEN EMPTY AND&
                  OCCUPIED!!!!! THIS ROUTINE CAN'T DEAL WITH PARTIAL&
                  OCCUPANCIES."
                  WRITE(*,*) "PLEASE CHECK THE HF GROUNDSTATE CALCULATION! "
                  WRITE(*,*) "SETTING SIGMA TO A SMALLER VALUE IN THE HF CALCULATION COULD HELP (SIGMA=0.000001)."
                  EXIT
               ENDIF
            ENDDO
         ENDDO
         ENDDO

      END SUBROUTINE INSULATOR_CHECK


!***********************************************************************
!This subroutine redistributes the Coulomb Vertex
!from the column-only process grid (CONTXT) to the 2D
!process grid (CONTXT_GRID)
!***********************************************************************

      SUBROUTINE REDISTRIBUTE_CVERTEX_1D(IO,WDES,KI,KQ,ISP,TMP_CVERTEX,NBI_START,NBI_END)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE (in_struct) IO
         INTEGER :: KI,KQ,ISP, nbi_start, nbi_end
         COMPLEX(q), TARGET :: TMP_CVERTEX(:,:,:)
         INTEGER :: CVERTEX_1D_ROWS, CVERTEX_1D_COLS,cc
         INTEGER :: NBI,MNBJ,NBJ,CVERTEX_ROWS_BR, MB, NB

         ! Prepare array descriptors for ScaLAPACK
         call BLACS_GRIDINFO(CONTXT_1D, NPROW_1D, NPCOL_1D, MYROW_1D, MYCOL_1D)
         MB=NGVECTOR !MIN(((NGVECTOR)/NPROW),50)   !Row blocking size
         CVERTEX_1D_rows = NGVECTOR !numroc(NGVECTOR,mb,MYROW,0,NPROW)
         DESC_CVERTEX_1D(1) = 1              ! descriptor type
         DESC_CVERTEX_1D(2) = CONTXT_1D         ! blacs context
         DESC_CVERTEX_1D(3) = NGVECTOR    ! global number of rows
         DESC_CVERTEX_1D(4) = (NBANDSDUMP) ! global number of cols
         DESC_CVERTEX_1D(5) = NGVECTOR     ! row block size
         DESC_CVERTEX_1D(6) = 1             ! col block size
         DESC_CVERTEX_1D(7) = 0              ! initial process row
         DESC_CVERTEX_1D(8) = SRCPROC              ! initial process col
         DESC_CVERTEX_1D(9) = MAX(1,(NGVECTOR)) ! leading dimension of local array
         !Distribute the fourier-transformed overlap integrals

         IF (ME==0) WRITE(IO%IU6,*)'nbi_start,nbi_end',(NBI_START+NFREEZE),(NBI_END+NFREEZE)

         DO NBI=NBI_START,NBI_END,1
            CALL PZGEMR2D(NGVECTOR,(NBANDSDUMP),TMP_CVERTEX(1,1,(NBI-nbi_start+1)),1,1,&
              DESC_CVERTEX_1D,CVERTEX_2D(1,1,NBI,RKIofKI(KI),RKQofKQ(KQ),ISP),1,1,DESC_CVERTEX_2D,CONTXT_1D)
         ENDDO


      END SUBROUTINE REDISTRIBUTE_CVERTEX_1D

!***********************************************************************
!This routine initializes a process grid which is used for the
!pre-calculation of the Coulomb Vertex
!This process grid has NPROCS(number of processors used) columns and 1 row.
!The assigned context variable is called CONTXT_COLS.
!***********************************************************************

      SUBROUTINE INIT_BLACS_1D()
         implicit none
         INTEGER :: a_PRCS, i
         REAL :: NPCOL_TMP

         !first we create a column-only process grid in order to
         !calculate the <i|-G|a> and <j|G|b> quantities
         CALL BLACS_PINFO(ME,PROCS)
         NPROW_1D=1
         NPCOL_1D=PROCS
         CALL BLACS_PINFO(ME,PROCS)
         CALL BLACS_GET     (0, 0, CONTXT_1D)
         CALL BLACS_GRIDINIT(CONTXT_1D, 'R', NPROW_1D, NPCOL_1D)
         CALL BLACS_GRIDINFO(CONTXT_1D, NPROW_1D, NPCOL_1D, MYROW_1D, MYCOL_1D)

      END SUBROUTINE INIT_BLACS_1D

!***********************************************************************
!This routine initializes the process grid, which is used for the
!matrix-matrix multiplications and their
!block-cyclic data distribution. Note that the context variable
!CONTXT_GRID2 is used for this grid.
!The routine tries to create the most quadratic grid possible for the given
!number of processors.
!***********************************************************************

      SUBROUTINE INIT_BLACS_2D(WDES)
         implicit none
         TYPE(wavedes) WDES
         INTEGER :: A_PRCS, i
         REAL :: NPCOL_TMP

         call BLACS_PINFO(ME,PROCS)
         
         !since, a quadratic process grid is going to be needed for
         !the matrix-matrix multiplications, its optimal row and column
         !numbers are estimated for the given number of processors
         A_PRCS=CEILING(SQRT(Real(PROCS,kind=q)))
         IF (a_PRCS==SQRT(Real(PROCS,kind=q))) THEN
            NPROW_2D=A_PRCS
            NPCOL_TMP=A_PRCS
         ENDIF
         IF (A_PRCS/=SQRT(Real(PROCS,kind=q))) THEN
            DO i=1,CEILING(SQRT(Real(PROCS,kind=q)))
               NPCOL_TMP=PROCS/Real(i,kind=q)
               IF ((NPCOL_TMP-INT(NPCOL_TMP))==0) THEN
                  NPROW_2D=i
               ENDIF
            ENDDO
            NPCOL_TMP=PROCS/NPROW_2D
         ENDIF
         NPCOL_2D=NPCOL_TMP
         IF (ABS(NPROW_2D-A_PRCS)>(A_PRCS/2.0)) THEN
            IF ((MYROW_1D==0) .AND. (MYCOL_1D==0)) WRITE(*,*)'The allocated number of CPUs does not allow for the use of an efficient process grid.'
            IF ((MYROW_1D==0) .AND. (MYCOL_1D==0)) WRITE(*,*)'Suggested number of CPUs: 4, 16, 32, ....'
         ENDIF

   
         CALL BLACS_PINFO(ME,PROCS)
         CALL BLACS_GET     (0, 0, CONTXT_2D)
         CALL BLACS_GRIDINIT(CONTXT_2D, 'R', NPROW_2D, NPCOL_2D)
         CALL BLACS_GRIDINFO(CONTXT_2D, NPROW_2D, NPCOL_2D, MYROW_2D, MYCOL_2D) ! all set up properly?
         
         IF ((MYROW_2D==0) .AND. (MYCOL_2D==0)) THEN
           WRITE(*,'(A,I3,A,I3,A)')'The allocated processors form a',NPROW_2D,'x',NPCOL_2D,' grid.'
         ENDIF

      END SUBROUTINE INIT_BLACS_2D


!***********************************************************************
! TEST_OPTCVERTEX performs a sanity check of the Coulomb Vertex and is useful for debugging.
!***********************************************************************
      SUBROUTINE TEST_OPTCVERTEX(W,WDES,INFO,IO)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE mkpoints
         USE lattice
         USE ini
         USE constant
         USE base
         IMPLICIT NONE
         TYPE(wavespin) W
         TYPE(wavedes) WDES
         TYPE (info_struct) INFO
         TYPE (in_struct) IO
         ! local variables
         INTEGER :: TSVEC_2D_ROWS, TSVEC_2D_COLS, SVM_2D_ROWS, SVM_2D_COLS
         INTEGER :: NG, NI, NJ, NA, NB
         COMPLEX(q) :: vijab, vijba, EMP2

         EMP2=0.0_q
         DO NI=1,4
         DO NJ=1,4
         DO NA=5,NBANDSDUMP
         DO NB=5,NBANDSDUMP

            vijab=(0.0_q,0.0_q)
            vijba=(0.0_q,0.0_q)
         DO NG=1,NGVECTOR
#ifdef gammareal
            vijab=vijab+REAL(CVERTEX_2D(NG,NA,NI,1,1,1),kind=8)*REAL(CVERTEX_2D(NG,NJ,NB,1,1,1),kind=8)
            vijab=vijab+AIMAG(CVERTEX_2D(NG,NA,NI,1,1,1))*AIMAG(CVERTEX_2D(NG,NJ,NB,1,1,1))
            vijba=vijba+REAL(CVERTEX_2D(NG,NB,NI,1,1,1),kind=8)*REAL(CVERTEX_2D(NG,NJ,NA,1,1,1),kind=8)
            vijba=vijba+AIMAG(CVERTEX_2D(NG,NB,NI,1,1,1))*AIMAG(CVERTEX_2D(NG,NJ,NA,1,1,1))
#else
            vijab=vijab+CVERTEX_2D(NG,NA,NI,1,1,1)*CONJG(CVERTEX_2D(NG,NJ,NB,1,1,1))
            vijba=vijba+CVERTEX_2D(NG,NB,NI,1,1,1)*CONJG(CVERTEX_2D(NG,NJ,NA,1,1,1))
#endif

         ENDDO
            EMP2=EMP2+vijab*conjg(2*vijab-vijba)/(W%CELTOT(NI,1,1)+W%CELTOT(NJ,1,1)-W%CELTOT(NA,1,1)-W%CELTOT(NB,1,1))
         ENDDO
         ENDDO
         ENDDO
         ENDDO

         WRITE(*,*)'EMP2=',EMP2


         EMP2=0.0_q
         DO NI=1,4
         DO NJ=1,4
         DO NA=5,NBANDSDUMP
         DO NB=5,NBANDSDUMP

            vijab=(0.0_q,0.0_q)
            vijba=(0.0_q,0.0_q)
         DO NG=1,NOPTAUX
#ifdef gammareal
            vijab=vijab+REAL(OPTCVERTEX_1D(NG,NA,NI,1,1,1),kind=8)*REAL(OPTCVERTEX_1D(NG,NJ,NB,1,1,1),kind=8)
            vijab=vijab+AIMAG(OPTCVERTEX_1D(NG,NA,NI,1,1,1))*AIMAG(OPTCVERTEX_1D(NG,NJ,NB,1,1,1))
            vijba=vijba+REAL(OPTCVERTEX_1D(NG,NB,NI,1,1,1),kind=8)*REAL(OPTCVERTEX_1D(NG,NJ,NA,1,1,1),kind=8)
            vijba=vijba+AIMAG(OPTCVERTEX_1D(NG,NB,NI,1,1,1))*AIMAG(OPTCVERTEX_1D(NG,NJ,NA,1,1,1))
#else
            vijab=vijab+OPTCVERTEX_1D(NG,NA,NI,1,1,1)*CONJG(OPTCVERTEX_1D(NG,NJ,NB,1,1,1))
            vijba=vijba+OPTCVERTEX_1D(NG,NB,NI,1,1,1)*CONJG(OPTCVERTEX_1D(NG,NJ,NA,1,1,1))
#endif
         ENDDO

            EMP2=EMP2+vijab*conjg(2*vijab-vijba)/(W%CELTOT(NI,1,1)+W%CELTOT(NJ,1,1)-W%CELTOT(NA,1,1)-W%CELTOT(NB,1,1))
         ENDDO
         ENDDO
         ENDDO
         ENDDO

         WRITE(*,*)'EMP2=',EMP2
      END SUBROUTINE TEST_OPTCVERTEX


!***********************************************************************
!
! Compute optimized auxilliary field for the Coulomb Vertex
! using a singular-value decomposition algorithm.
!
!***********************************************************************

      SUBROUTINE CALC_OPTCVERTEX(W,WDES,INFO,IO)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE mkpoints
         USE lattice
         USE ini
         USE constant
         USE base
         IMPLICIT NONE
         TYPE(wavespin) W
         TYPE(wavedes) WDES
         TYPE (info_struct) INFO
         TYPE (in_struct) IO
         ! local variables
! data needed for sv decomposition
         INTEGER :: LWORK, LRWORK, LIWORK, PINFO
         COMPLEX(q), ALLOCATABLE ::  WORK(:)
         COMPLEX(q), ALLOCATABLE :: RWORK(:)
         INTEGER, ALLOCATABLE :: IWORK(:)
!!
         INTEGER :: TSVEC_2D_ROWS, TSVEC_2D_COLS, SVM_2D_ROWS, SVM_2D_COLS
         INTEGER :: NP, KP, KQ,ISP, I, NG, NGP
         INTEGER :: RNG, RNGP
#ifdef gammareal
         REAL(q), ALLOCATABLE :: RSVM_2D(:,:), RSVALUES(:), RSVEC_2D(:,:)
         REAL(q), ALLOCATABLE ::  REALWORK(:)
#endif

! first build SVM by contracting coulomb vertices

         call BLACS_GRIDINFO(CONTXT_2D, NPROW_2D, NPCOL_2D, MYROW_2D, MYCOL_2D)
         CALL SETUP_SV_DATA_2D(WDES)
         SVM_2D=(0.0_q,0.0_q)

         DO ISP=1,WDES%ISPIN
         DO KP=1,REALNKPTS
         DO KQ=1,REALNKPTS
         DO NP=1,NBANDSDUMP

#ifdef gammareal
         !contract ftod over ab and add to G,G'-Matrix
            CALL PZGEMM('n','c',(NGVECTOR),(NGVECTOR),&
                NBANDSDUMP,(1.0_q,0.0_q), CVERTEX_2D(1,1,NP,KP,KQ,ISP),1,1,&
                DESC_CVERTEX_2D, CVERTEX_2D(1,1,NP,KP,KQ,ISP),1,1,&
                DESC_CVERTEX_2D,(1.0_q,0.0_q), SVM_2D(1,1),1,1,DESC_SVM_2D)

#else
         !contract ftod over ab and add to G,G'-Matrix
            CALL PZGEMM('n','c',(NGVECTOR),(NGVECTOR),&
                NBANDSDUMP,(1.0_q,0.0_q), CVERTEX_2D(1,1,NP,KP,KQ,ISP),1,1,&
                DESC_CVERTEX_2D, CVERTEX_2D(1,1,NP,KP,KQ,ISP),1,1,&
                DESC_CVERTEX_2D,(1.0_q,0.0_q), SVM_2D(1,1),1,1,DESC_SVM_2D)

#endif

         ENDDO
         ENDDO
         ENDDO
         ENDDO

#ifdef gammareal
         SVM_2D=CMPLX(REAL(SVM_2D), 0.0_q)
#else
         SVM_2D=CONJG(SVM_2D)
#endif

         IF (ALLOCATED(WORK)) DEALLOCATE(WORK)
         IF (ALLOCATED(RWORK)) DEALLOCATE(RWORK)
         IF (ALLOCATED(IWORK)) DEALLOCATE(IWORK)
         IF (ALLOCATED(SVALUES)) DEALLOCATE(SVALUES)
         ALLOCATE(SVALUES(DESC_SVM_2D(3)))
         ALLOCATE(WORK(1))
         ALLOCATE(RWORK(1))
         ALLOCATE(IWORK(1))
         WORK=(0.0_q,0.0_q)
         RWORK=(0.0_q,0.0_q)
         IWORK=0
         SVALUES=0.0_q
         SVEC_2D=(0.0_q,0.0_q)

! #ifdef gammareal
!        ! WRITE(*,*)'before alloc0'
!!         WRITE(*,*)'SVM_2D',SVM_2D
!         ALLOCATE(REALWORK(1))
!         ALLOCATE(RSVALUES(SIZE(SVALUES,1)))
!         ALLOCATE(RSVM_2D(SIZE(SVM_2D,1),SIZE(SVM_2D,2)))
!         ALLOCATE(RSVEC_2D(SIZE(SVEC_2D,1),SIZE(SVEC_2D,2)))
!         RSVM_2D(:,:)=REAL(REAL(SVM_2D(:,:)))
!         RSVEC_2D=0.0_q
!         REALWORK=zero

        
!         CALL PDSYEV ( 'V', 'U', NGVECTOR, RSVM_2D(1,1), 1, 1, DESC_SVM_2D,RSVALUES(1),&
!                RSVEC_2D(1,1), 1, 1, DESC_SVM_2D, REALWORK(1), -1,PINFO )
!        !write svm and svec and svalues to complex arrays
! #else
! prepare singular values decomposition of SVM
         CALL PZHEEVD( 'V', 'L', NGVECTOR, SVM_2D(1,1), 1, 1, DESC_SVM_2D, SVALUES(1),&
                SVEC_2D(1,1), 1, 1, DESC_SVM_2D, WORK(1),-1,RWORK(1),-1,IWORK(1), LIWORK, &
                PINFO)
! #endif


         CALL START_TIMING("PD")

! #ifdef gammareal
!        LWORK=INT(REALWORK(1))
!        DEALLOCATE(REALWORK)
!        ALLOCATE(REALWORK(LWORK))
!        REALWORK=zero
!        CALL PDSYEV( 'V', 'U', NGVECTOR, RSVM_2D(1,1), 1, 1,DESC_SVM_2D(1), RSVALUES(1),&
!                RSVEC_2D(1,1), 1, 1, DESC_SVM_2D(1), REALWORK,LWORK,PINFO)
!        !write svm and svec and svalues to complex arrays
!        SVEC_2D(:,:)=CMPLX(RSVEC_2D(:,:),0.0_q)
!!        SVM_2D(:,:)=CMPLX(RSVM_2D(:,:),0.0_q)
!        SVALUES(:)=RSVALUES(:)
! #else

         LWORK=INT(WORK(1))
         LRWORK=INT(RWORK(1)) !2*NGVECTOR+2*NGVECTOR-2
         LIWORK=INT(IWORK(1))
         DEALLOCATE(WORK)
         DEALLOCATE(RWORK)
         DEALLOCATE(IWORK)
         ALLOCATE(WORK(LWORK))
         ALLOCATE(RWORK(LRWORK))
         ALLOCATE(IWORK(LIWORK))
         WORK=(0.0_q,0.0_q)
         RWORK=(0.0_q,0.0_q)
         IWORK=0



! do singular values decomposition of SVM
         CALL PZHEEVD( 'V', 'L', NGVECTOR, SVM_2D(1,1), 1, 1, DESC_SVM_2D(1), SVALUES(1),&
                SVEC_2D(1,1), 1, 1, DESC_SVM_2D(1), WORK, LWORK,RWORK,LRWORK,IWORK,LIWORK,PINFO)
! #endif

         CALL STOP_TIMING("PD",IO%IU6,"PZHEEVD")
 
! write singular values
!         IF ((IO%IU6>0)) THEN
!            OPEN(unit = 9,file = "SVALUES_OPTIMIZED_AUXFIELD")
!            DO NG=1,NGVECTOR
!               write(9,*)NG,SVALUES(NGVECTOR-NG+1)
!            ENDDO
!            CLOSE(9)
!         ENDIF

! determine NOPTAUX according to EDIFF cutoff variable and singular values
         NOPTAUX=NGVECTOR
         DO I=1,NGVECTOR
            IF (SVALUES(I)<INFO%EDIFF) NOPTAUX=NGVECTOR-I
         ENDDO

         NOPTAUX=NGVECTOR

!         NOPTAUX=NGVECTOR
!         IF (IO%IU0>=0) write(IO%IU0,*)'Truncating all singular values below ',INFO%EDIFF,'results in ',NOPTAUX,' optimal auxiliary field vectors.'

! Allocate coulomb vertex for optimized aux field now
         CALL SETUP_TMPCVERTEX_2D(WDES)
         CALL SETUP_TSVEC_2D(WDES)

! unfortunaetly the singular vectors are given in ascending order and we need them in descending order
! there are probably better ways of changing the order than the following one:
         SVM_2D(:,:)=SVEC_2D(:,:)
         SVM_TMP_2D=(0.0_q,0.0_q)

         SVM_2D_ROWS = numroc(NGVECTOR,DESC_SVM_2D(5),MYROW_2D,0,NPROW_2D)
         SVM_2D_COLS = numroc(NGVECTOR,DESC_SVM_2D(6),MYCOL_2D,0,NPCOL_2D)
         DO NG=1,SVM_2D_ROWS
            DO NGP=1,SVM_2D_COLS
               CALL LOC2GLOB(NG,MYROW_2D,DESC_SVM_2D(3),NPROW_2D,DESC_SVM_2D(5),RNG)
               CALL LOC2GLOB(NGP,MYCOL_2D,DESC_SVM_2D(4),NPCOL_2D,DESC_SVM_2D(6),RNGP)
               IF ((-RNG+(NGVECTOR)+1)==RNGP) SVM_TMP_2D(NG,NGP)=(1.0_q,0.0_q)
            ENDDO
         ENDDO

! bring singular vectors in descending order w.r.t. the singular values
         CALL PZGEMM('n','n',NGVECTOR,NGVECTOR,&
                NGVECTOR,(1.0_q,0.0_q),SVM_2D(1,1),1,1,&
                DESC_SVM_2D,SVM_TMP_2D(1,1),1,1,&
                DESC_SVM_2D,(0.0_q,0.0_q), SVEC_2D(1,1),1,1,DESC_SVM_2D)


         TSVEC_2D_ROWS = numroc(DESC_TSVEC_2D(3),DESC_TSVEC_2D(5),MYROW_2D,0,NPROW_2D)
         TSVEC_2D_COLS = numroc(DESC_TSVEC_2D(4),DESC_TSVEC_2D(6),MYCOL_2D,0,NPCOL_2D)

         ! write slice of SVEC_2D to TSVEC_2D
         TSVEC_2D(1:TSVEC_2D_ROWS,1:TSVEC_2D_COLS)=SVEC_2D(1:TSVEC_2D_ROWS,1:TSVEC_2D_COLS)

! ag test : use this code block to bypass transformation into aux-field
!           representation of coulomb vertex
!         TSVEC_2D=(0.0_q,0.0_q)
!         DO NG=1,TSVEC_2D_ROWS
!            DO NGP=1,TSVEC_2D_COLS
!               CALL LOC2GLOB(NG,MYROW_2D,DESC_TSVEC_2D(3),NPROW_2D,DESC_TSVEC_2D(5),RNG)
!               CALL LOC2GLOB(NGP,MYCOL_2D,DESC_TSVEC_2D(4),NPCOL_2D,DESC_TSVEC_2D(6),RNGP)
!               IF (RNG==RNGP) TSVEC_2D(NG,NGP)=(1.0_q,0.0_q)
!            ENDDO
!         ENDDO
! ag test end



! transpose TSVEC_2D and redistribute using Scalapack.
         CALL SETUP_TSVEC_1D_USING_TSVEC_2D(WDES)

! transform CVERTEX_2D using TSVEC_2D to OPTCVERTEX2D
         DO ISP=1,WDES%ISPIN
         DO KQ=1,REALNKPTS
         DO KP=1,REALNKPTS
         DO NP=1,NBANDSDUMP

            CALL PZGEMM('t','n',(NOPTAUX),(NBANDSDUMP),&
                NGVECTOR,(1.0_q,0.0_q), TSVEC_2D(1,1),1,1,&
                DESC_TSVEC_2D, CVERTEX_2D(1,1,NP,KP,KQ,ISP),1,1,&
                DESC_CVERTEX_2D,(0.0_q,0.0_q), TMPCVERTEX_2D(1,1),1,1,DESC_TMPCVERTEX_2D)

            CVERTEX_2D(:,:,NP,KP,KQ,ISP)=TMPCVERTEX_2D(:,:)

         ENDDO
         ENDDO
         ENDDO
         ENDDO

      END SUBROUTINE CALC_OPTCVERTEX
 

      SUBROUTINE SETUP_TSVEC_2D(WDES)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: TSVEC_2D_ROWS, TSVEC_2D_COLS, NI, TMPMYCOL, RNI, MB, NB

         call BLACS_GRIDINFO(CONTXT_2D, NPROW_2D, NPCOL_2D, MYROW_2D, MYCOL_2D)
         MB=MIN((NGVECTOR/NPROW_2D),(NGVECTOR/NPCOL_2D))   !Row blocking size
         NB=MIN((NGVECTOR/NPCOL_2D),MB)   !Column blocking size , make squares


         TSVEC_2D_rows = numroc(NGVECTOR,mb,MYROW_2D,0,NPROW_2D)
         TSVEC_2D_cols = numroc(NOPTAUX,nb,MYCOL_2D,0,NPCOL_2D)
         DESC_TSVEC_2D(1) = 1              ! descriptor type
         DESC_TSVEC_2D(2) = contxt_2D         ! blacs context
         DESC_TSVEC_2D(3) = NGVECTOR           ! global number of rows
         DESC_TSVEC_2D(4) = NOPTAUX   ! global number of cols
         DESC_TSVEC_2D(5) = Mb       ! row block size
         DESC_TSVEC_2D(6) = Nb       ! column block size
         DESC_TSVEC_2D(7) = 0              ! initial process row
         DESC_TSVEC_2D(8) = 0              ! initial process column
         DESC_TSVEC_2D(9) = MAX(1,TSVEC_2D_ROWS) ! leading dimension of local array

         ALLOCATE(TSVEC_2D(TSVEC_2D_rows,TSVEC_2D_cols))

      END SUBROUTINE SETUP_TSVEC_2D

      SUBROUTINE SETUP_TSVEC_1D_USING_TSVEC_2D(WDES)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: TSVEC_1D_ROWS, TSVEC_1D_COLS, NI, TMPMYCOL, RNI, MB, NB

         call BLACS_GRIDINFO(CONTXT_1D, NPROW_1D, NPCOL_1D, MYROW_1D, MYCOL_1D)
         MB=NGVECTOR   !Row blocking size
         NB=1   !Column blocking size , make squares

         TSVEC_1D_rows = numroc(NGVECTOR,mb,MYROW_1D,0,NPROW_1D)
         TSVEC_1D_cols = numroc(NOPTAUX,nb,MYCOL_1D,0,NPCOL_1D)
         DESC_TSVEC_1D(1) = 1              ! descriptor type
         DESC_TSVEC_1D(2) = contxt_1D         ! blacs context
         DESC_TSVEC_1D(3) = NGVECTOR           ! global number of rows
         DESC_TSVEC_1D(4) = NOPTAUX   ! global number of cols
         DESC_TSVEC_1D(5) = Mb       ! row block size
         DESC_TSVEC_1D(6) = Nb       ! column block size
         DESC_TSVEC_1D(7) = 0              ! initial process row
         DESC_TSVEC_1D(8) = 0              ! initial process column
         DESC_TSVEC_1D(9) = MAX(1,TSVEC_1D_ROWS) ! leading dimension of local array

         ALLOCATE(TSVEC_1D(TSVEC_1D_rows,TSVEC_1D_cols))

         CALL PZGEMR2D(NGVECTOR,NOPTAUX,TSVEC_2D(1,1),1,1,&
           DESC_TSVEC_2D,TSVEC_1D(1,1),1,1,DESC_TSVEC_1D,CONTXT_1D)

      END SUBROUTINE SETUP_TSVEC_1D_USING_TSVEC_2D


!***********************************************************************
!Allocate the arrays for the Coulomb Vertex singular vectors (SV_2d )
!and setup their descriptors which are required by the scalaLAPACK routines
!***********************************************************************

      SUBROUTINE SETUP_SV_DATA_2D(WDES)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: SVM_2D_ROWS, SVM_2D_COLS, NI, TMPMYCOL, RNI, MB, NB

         call BLACS_GRIDINFO(CONTXT_2D, NPROW_2D, NPCOL_2D, MYROW_2D, MYCOL_2D)
         MB=MIN((NGVECTOR/NPROW_2D),(NGVECTOR/NPCOL_2D))   !Row blocking size
         NB=MIN((NGVECTOR/NPCOL_2D),MB)   !Column blocking size , make squares

         SVM_2D_rows = numroc(NGVECTOR,mb,MYROW_2D,0,NPROW_2D)
         SVM_2D_cols = numroc(NGVECTOR,nb,MYCOL_2D,0,NPCOL_2D)
         If ((MYROW_2D==0) .AND. (MYCOL_2D==0)) THEN
         ENDIF
         DESC_SVM_2D(1) = 1              ! DESCriptor type
         DESC_SVM_2D(2) = contxt_2D         ! blacs context
         DESC_SVM_2D(3) = NGVECTOR       ! global number of rows
         DESC_SVM_2D(4) = NGVECTOR   ! global number of cols
         DESC_SVM_2D(5) = mb       ! row block size
         DESC_SVM_2D(6) = nb       ! column block size
         DESC_SVM_2D(7) = 0              ! initial process row
         DESC_SVM_2D(8) = 0              ! initial process column
         DESC_SVM_2D(9) = MAX(1,SVM_2D_ROWS) ! leading dimension of local array

         ALLOCATE(    SVM_2D(SVM_2D_rows,SVM_2D_cols))
         ALLOCATE(SVM_TMP_2D(SVM_2D_rows,SVM_2D_cols))
         ALLOCATE(   SVEC_2D(SVM_2D_rows,SVM_2D_cols))

      END SUBROUTINE SETUP_SV_DATA_2D


!***********************************************************************
!Allocate the arrays to store the optimized Coulomb Vertex in aux-field basis
!in the block-cyclic data distribution
!and setup their descriptors which are required by the scalaLAPACK routines
!***********************************************************************

      SUBROUTINE SETUP_OPTCVERTEX_1D(WDES)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: OPTCVERTEX_1D_rows, OPTCVERTEX_1D_cols,cc,NI, TMPMYCOL, RNI, MB, NB

         call BLACS_GRIDINFO(CONTXT_1D, NPROW_1D, NPCOL_1D, MYROW_1D, MYCOL_1D)
         MB=NOPTAUX !MIN(((NGVECTOR)/NPROW),50)   !Row blocking size
         NB=1 !MIN((NPROW*WDES%NBANDS),MB)   !column blocking size

         OPTCVERTEX_1D_rows = numroc(NOPTAUX,mb,MYROW_1D,0,NPROW_1D)
         OPTCVERTEX_1D_cols = numroc(NBANDSDUMP,nb,MYCOL_1D,0,NPCOL_1D)
         DESC_OPTCVERTEX_1D(1) = 1              ! DESCriptor type
         DESC_OPTCVERTEX_1D(2) = CONTXT_1D         ! blacs context
         DESC_OPTCVERTEX_1D(3) = NOPTAUX       ! global number of rows
         DESC_OPTCVERTEX_1D(4) = (NBANDSDUMP)   ! global number of cols
         DESC_OPTCVERTEX_1D(5) = MB       ! row block size
         DESC_OPTCVERTEX_1D(6) = NB       ! column block size
         DESC_OPTCVERTEX_1D(7) = 0              ! initial process row
         DESC_OPTCVERTEX_1D(8) = 0              ! initial process column
         DESC_OPTCVERTEX_1D(9) = MAX(1,OPTCVERTEX_1D_ROWS) ! leading dimension of local array

         ALLOCATE(OPTCVERTEX_1D(OPTCVERTEX_1D_ROWS,OPTCVERTEX_1D_COLS,NBANDSDUMP,REALNKPTS,REALNKPTS,WDES%ISPIN))

         IF (.not. ALLOCATED(PRO_NI)) ALLOCATE(PRO_NI((NBANDSDUMP)))

         PRO_NI=-99
         DO NI=1,(NBANDSDUMP)
            CALL GLOB2LOC(NI,DESC_OPTCVERTEX_1D(4),NPCOL_1D,DESC_OPTCVERTEX_1D(6),TMPMYCOL,RNI)
            IF (TMPMYCOL==MYCOL_1D) THEN
               PRO_NI(NI)=ME
            ELSE
               PRO_NI(NI)=0
            ENDIF
         ENDDO
         CALL M_sum_i(WDES%COMM_INTER, PRO_NI,(NBANDSDUMP) )

      END SUBROUTINE SETUP_OPTCVERTEX_1D

      SUBROUTINE SETUP_TMPCVERTEX_1D(WDES)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: TMPCVERTEX_1D_rows, TMPCVERTEX_1D_cols,cc,NI, TMPMYCOL, RNI, MB, NB

         call BLACS_GRIDINFO(CONTXT_1D, NPROW_1D, NPCOL_1D, MYROW_1D, MYCOL_1D)
         MB=NGVECTOR !MIN(((NGVECTOR)/NPROW),50)   !Row blocking size
         NB=1 !MIN((NPROW*WDES%NBANDS),MB)   !column blocking size

         TMPCVERTEX_1D_rows = numroc(NGVECTOR,mb,MYROW_1D,0,NPROW_1D)
         TMPCVERTEX_1D_cols = numroc(NBANDSDUMP,nb,MYCOL_1D,0,NPCOL_1D)
         DESC_TMPCVERTEX_1D(1) = 1              ! DESCriptor type
         DESC_TMPCVERTEX_1D(2) = CONTXT_1D         ! blacs context
         DESC_TMPCVERTEX_1D(3) = NGVECTOR       ! global number of rows
         DESC_TMPCVERTEX_1D(4) = (NBANDSDUMP)   ! global number of cols
         DESC_TMPCVERTEX_1D(5) = MB       ! row block size
         DESC_TMPCVERTEX_1D(6) = NB       ! column block size
         DESC_TMPCVERTEX_1D(7) = 0              ! initial process row
         DESC_TMPCVERTEX_1D(8) = 0              ! initial process column
         DESC_TMPCVERTEX_1D(9) = MAX(1,TMPCVERTEX_1D_ROWS) ! leading dimension of local array

         ALLOCATE(TMPCVERTEX_1D(TMPCVERTEX_1D_ROWS,TMPCVERTEX_1D_COLS))

         IF (.not. ALLOCATED(PRO_NI)) ALLOCATE(PRO_NI((NBANDSDUMP)))

         PRO_NI=-99
         DO NI=1,(NBANDSDUMP)
            CALL GLOB2LOC(NI,DESC_TMPCVERTEX_1D(4),NPCOL_1D,DESC_TMPCVERTEX_1D(6),TMPMYCOL,RNI)
            IF (TMPMYCOL==MYCOL_1D) THEN
               PRO_NI(NI)=ME
            ELSE
               PRO_NI(NI)=0
            ENDIF
         ENDDO
         CALL M_sum_i(WDES%COMM_INTER, PRO_NI,(NBANDSDUMP) )

      END SUBROUTINE SETUP_TMPCVERTEX_1D


      SUBROUTINE SETUP_TMPCVERTEX_2D(WDES)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: TMPCVERTEX_2D_rows, TMPCVERTEX_2D_cols, MB, NB

         call BLACS_GRIDINFO(CONTXT_2D, NPROW_2D, NPCOL_2D, MYROW_2D, MYCOL_2D)
         MB=MIN(((NGVECTOR)/NPROW_2D),50)   !Row blocking size
         NB=MIN((NBANDSDUMP/NPCOL_2D),MB)   !column blocking size

         TMPCVERTEX_2D_rows = numroc(NGVECTOR,mb,MYROW_2D,0,NPROW_2D)
         TMPCVERTEX_2D_cols = numroc(NBANDSDUMP,nb,MYCOL_2D,0,NPCOL_2D)
         DESC_TMPCVERTEX_2D(1) = 1              ! descriptor type
         DESC_TMPCVERTEX_2D(2) = CONTXT_2D         ! blacs context
         DESC_TMPCVERTEX_2D(3) = NGVECTOR       ! global number of rows
         DESC_TMPCVERTEX_2D(4) = (NBANDSDUMP)   ! global number of cols
         DESC_TMPCVERTEX_2D(5) = mb       ! row block size
         DESC_TMPCVERTEX_2D(6) = nb       ! column block size
         DESC_TMPCVERTEX_2D(7) = 0              ! initial process row
         DESC_TMPCVERTEX_2D(8) = 0              ! initial process column
         DESC_TMPCVERTEX_2D(9) = MAX(1,TMPCVERTEX_2D_rows) ! leading dimension of local array

         ALLOCATE(TMPCVERTEX_2D(TMPCVERTEX_2D_rows,TMPCVERTEX_2D_cols))

      END SUBROUTINE SETUP_TMPCVERTEX_2D

      SUBROUTINE SETUP_CVERTEX_2D(WDES)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: CVERTEX_2D_rows, CVERTEX_2D_cols, MB, NB

         call BLACS_GRIDINFO(CONTXT_2D, NPROW_2D, NPCOL_2D, MYROW_2D, MYCOL_2D)
         MB=MIN(((NGVECTOR)/NPROW_2D),50)   !Row blocking size
         NB=MIN((NBANDSDUMP/NPCOL_2D),MB)   !column blocking size

         CVERTEX_2D_rows = numroc(NGVECTOR,MB,MYROW_2D,0,NPROW_2D)
         CVERTEX_2D_cols = numroc(NBANDSDUMP,NB,MYCOL_2D,0,NPCOL_2D)
         DESC_CVERTEX_2D(1) = 1              ! DESCriptor type
         DESC_CVERTEX_2D(2) = CONTXT_2D         ! blacs context
         DESC_CVERTEX_2D(3) = NGVECTOR       ! global number of rows
         DESC_CVERTEX_2D(4) = (NBANDSDUMP)   ! global number of cols
         DESC_CVERTEX_2D(5) = mb       ! row block size
         DESC_CVERTEX_2D(6) = nb       ! column block size
         DESC_CVERTEX_2D(7) = 0              ! initial process row
         DESC_CVERTEX_2D(8) = 0              ! initial process column
         DESC_CVERTEX_2D(9) = MAX(1,CVERTEX_2D_rows) ! leading dimension of local array

         ALLOCATE(CVERTEX_2D(CVERTEX_2D_rows,CVERTEX_2D_cols,NBANDSDUMP,REALNKPTS,REALNKPTS,WDES%ISPIN))

      END SUBROUTINE SETUP_CVERTEX_2D


!***********************************************************************
! TRANSFORM LOCAL TO GLOBAL ARRAY INDEX FOR BLOCK-CYCLIC ARRAY
! DISTRIBUTION
!***********************************************************************
      SUBROUTINE LOC2GLOB(li,p,n,np,nb,gi)
         IMPLICIT NONE
         INTEGER :: li   ! local index
         INTEGER :: p    ! index in processor grid (either MYROW or MYCOL)
         INTEGER :: n    ! global array dimension
         INTEGER :: np   ! processor array dimension
         INTEGER :: nb   ! blocking size
         INTEGER :: gi   ! global index
         INTEGER :: litmp

         litmp = li-1
         gi=(((litmp/nb)*np)+p)*nb+mod(litmp,nb)+1
         RETURN
      END SUBROUTINE LOC2GLOB

!***********************************************************************
! TRANSFORM GLOBAL TO LOCAL ARRAY INDEX FOR BLOCK-CYCLIC ARRAY
! DISTRIBUTION
!***********************************************************************
      SUBROUTINE GLOB2LOC(i,n,np,nb,p,il)
         implicit none
         integer :: i    ! global array index, input
         integer :: n    ! global array dimension, input
         integer :: np   ! processor array dimension, input
         integer :: nb   ! block size, input
         integer :: p    ! processor array index, output
         integer :: il   ! local array index, output
         integer :: im1

         im1 = i-1
         p   = mod((im1/nb),np)
         il  = (im1/(np*nb))*nb + mod(im1,nb) + 1

         RETURN
      END SUBROUTINE


!***********************************************************************
!
! Write Cc4s v0.5 interface files for disk 
!
!***********************************************************************
      SUBROUTINE WRITEBINCVERTEX2DISK_FTODDUMP(W,WDES)
         use mkpoints
         IMPLICIT NONE
         type(wavespin) :: W
         TYPE(wavedes) WDES
         INTEGER :: NKI,NKJ,NKA,NI,NJ,NA,NKQ,NB,NKB,NG,NKK,NK
         INTEGER :: I,J,A,B,ISP1,ISP2,ISP3,ISP4,SP
         REAL(q) :: FA,FB
         COMPLEX(q), ALLOCATABLE :: CVERTEX_TMP(:,:,:), CVERTEX_TMP_SINGLE(:)
         COMPLEX(q) :: vijab,vijba, emp2
         LOGICAL :: ex, CONSISTENT
         INTEGER :: nbstart, nastart
         INTEGER*8 :: JUNKSIZE,IREC
         CHARACTER*8 :: MAGIC
         integer :: TMPMYCOL,MYJ
         INTEGER :: SPA,SPB,SI,SJ,NOCC
         REAL(q) :: EFERMI

         NOCC=0
         DO SI=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)
            SPA=SP_ORD(SI)
            I=N_ORD(SI)

            IF ((W%FERTOT(I,1,SPA)>0.001_q)) THEN
               NOCC=NOCC+1
            ENDIF
         ENDDO

         ALLOCATE(CVERTEX_TMP_SINGLE(NOPTAUX))
         ALLOCATE(CVERTEX_TMP(NOPTAUX,(NBANDSDUMP),WDES%ISPIN))

         IF (ME==0)  open(unit = 7,file = "FTODDUMP",FORM='UNFORMATTED',access='stream',STATUS='REPLACE')

         MAGIC='cc4sFTOD'
         IREC=1
         IF (ME==0) WRITE(7) MAGIC
         IREC=IREC+1
         IF (ME==0) WRITE(7) (NOCC-NFREEZE),INT(NBANDSDUMP)*WDES%ISPIN-(NOCC-NFREEZE*WDES%ISPIN)
         IREC=IREC+1
         IF (ME==0) WRITE(7) NOPTAUX,INT(WDES%ISPIN)
         IREC=IREC+1
         IF (ME==0) WRITE(7) INT(REALNKPTS),INT(REALNKPTS)
         IREC=IREC+1

         IF (ME==0) WRITE(7) 'FTODreal'
         IREC=IREC+1
         JUNKSIZE=(2+(((NBANDSDUMP)*WDES%ISPIN)**2*NOPTAUX))*8
         IF (ME==0) WRITE(7) JUNKSIZE
         IREC=IREC+1

         DO SI=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)

            SPA=SP_ORD(SI)
            I=N_ORD(SI)

            IF (I <= NFREEZE) CYCLE
            CVERTEX_TMP=(0.0_q,0.0_q)

            DO J=1,(NBANDSDUMP+NFREEZE)
               IF (J <= NFREEZE) CYCLE
               IF (PRO_NI(J-NFREEZE)==ME) THEN
                  CALL GLOB2LOC(J-NFREEZE,DESC_OPTCVERTEX_1D(4),NPCOL_1D,DESC_OPTCVERTEX_1D(6),TMPMYCOL,MYJ)
                  CVERTEX_TMP(:,J-NFREEZE,SPA)=OPTCVERTEX_1D(:,MYJ,I-NFREEZE,1,1,SPA)
               ENDIF
            ENDDO
            CALL M_sum_z(WDES%COMM_INTER, CVERTEX_TMP(1,1,SPA),NOPTAUX*(NBANDSDUMP))

            IF (ME==0) THEN
               DO SJ=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)
                  SPB=SP_ORD(SJ)
                  J=N_ORD(SJ)
                  IF (J <= NFREEZE) CYCLE
                  CVERTEX_TMP_SINGLE(:)=(0.0_q,0.0_q)
                  IF (SPA==SPB) CVERTEX_TMP_SINGLE(:)=CVERTEX_TMP(:,J-NFREEZE,SPA)
                     WRITE(7) REAL(CVERTEX_TMP_SINGLE(:),kind=8)
               ENDDO
            ENDIF

         ENDDO


         IF (ME==0) WRITE(7) 'FTODimag'
         IREC=IREC+1
         JUNKSIZE=(2+(((NBANDSDUMP)*WDES%ISPIN)**2*NOPTAUX))*8
         IF (ME==0) WRITE(7) JUNKSIZE
         IREC=IREC+1

         DO SI=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)

            SPA=SP_ORD(SI)
            I=N_ORD(SI)

            IF (I <= NFREEZE) CYCLE
            CVERTEX_TMP=(0.0_q,0.0_q)

            DO J=1,(NBANDSDUMP+NFREEZE)
               IF (J <= NFREEZE) CYCLE
               IF (PRO_NI(J-NFREEZE)==ME) THEN
                  CALL GLOB2LOC(J-NFREEZE,DESC_OPTCVERTEX_1D(4),NPCOL_1D,DESC_OPTCVERTEX_1D(6),TMPMYCOL,MYJ)
                  CVERTEX_TMP(:,J-NFREEZE,SPA)=OPTCVERTEX_1D(:,MYJ,I-NFREEZE,1,1,SPA)
               ENDIF
            ENDDO
            CALL M_sum_z(WDES%COMM_INTER, CVERTEX_TMP(1,1,SPA),NOPTAUX*(NBANDSDUMP))

            IF (ME==0) THEN
               DO SJ=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)
                  SPB=SP_ORD(SJ)
                  J=N_ORD(SJ)
                  IF (J <= NFREEZE) CYCLE
                  CVERTEX_TMP_SINGLE(:)=(0.0_q,0.0_q)
                  IF (SPA==SPB) CVERTEX_TMP_SINGLE(:)=CVERTEX_TMP(:,J-NFREEZE,SPA)
                     WRITE(7) REAL(DIMAG(CVERTEX_TMP_SINGLE(:)),kind=8)
               ENDDO
            ENDIF

         ENDDO


         IF (ME==0) WRITE(7) 'FTODepsi'
         IREC=IREC+1
         JUNKSIZE=(2+(NBANDSDUMP)*WDES%ISPIN)*8
         IF (ME==0) WRITE(7) JUNKSIZE
         IREC=IREC+1

         DO SI=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)
            SPA=SP_ORD(SI)
            I=N_ORD(SI)
            IF (I <= NFREEZE) CYCLE
            IF (ME==0) WRITE(7) REAL(W%CELTOT(I,1,SPA),kind=8)
            IREC=IREC+1
         ENDDO



         IF (ME==0) CLOSE(7)

      END SUBROUTINE WRITEBINCVERTEX2DISK_FTODDUMP

!***********************************************************************
!
! Write several input files need by Cc4s to disk including:
! CoulombVertex
! EigenEnergies
! States
! AuxiliaryField
!
!***********************************************************************

      SUBROUTINE WRITEBINCVERTEX2DISK_YAML(W,WDES)
         use mkpoints
         IMPLICIT NONE
         type(wavespin) :: W
         TYPE(wavedes) WDES
         INTEGER :: NKI,NKJ,NKA,NI,NJ,NA,NKQ,NB,NKB,NG,NKK,NK
         INTEGER :: I,J,A,B,ISP1,ISP2,ISP3,ISP4,SP
         REAL(q) :: FA,FB
         COMPLEX(q), ALLOCATABLE :: CVERTEX_TMP(:,:,:), CVERTEX_TMP_SINGLE(:)
         COMPLEX(q) :: vijab,vijba, emp2
         LOGICAL :: ex, CONSISTENT
         INTEGER :: nbstart, nastart
         INTEGER*8 :: JUNKSIZE,IREC
         CHARACTER*8 :: MAGIC
         integer :: TMPMYCOL,MYJ
         INTEGER :: SPA,SPB,SI,SJ,NOCC
         REAL(q) :: EFERMI

         NOCC=0
         DO SI=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)
            SPA=SP_ORD(SI)
            I=N_ORD(SI)

            IF ((W%FERTOT(I,1,SPA)>0.001_q)) THEN
               NOCC=NOCC+1
            ENDIF
         ENDDO
         EFERMI=(W%CELTOT(NOCC,1,SPA)+W%CELTOT(NOCC+1,1,SPA))/2.0_q/0.0367493221756387
         EFERMI=(W%CELTOT(NOCC,1,SPA)+W%CELTOT(NOCC+1,1,SPA))/2.0_q

         ALLOCATE(CVERTEX_TMP_SINGLE(NOPTAUX))
         ALLOCATE(CVERTEX_TMP(NOPTAUX,(NBANDSDUMP),WDES%ISPIN))


         IF (ME==0) write(*,*)'- Writing CoulombVertex .'
         IF (ME==0)  THEN

            OPEN(unit = 7,file = "CoulombVertex.yaml",FORM='FORMATTED',access='stream',STATUS='REPLACE')
            WRITE(7,'(A)')    'version: 100'
            WRITE(7,'(A)')    'type: Tensor'
            WRITE(7,'(A)')    'scalarType: Complex64'
            WRITE(7,'(A)')    'dimensions:'
            WRITE(7,'(A,I6)') '- length: ',NOPTAUX
            WRITE(7,'(A)')    '  type: AuxiliaryField'
            WRITE(7,'(A,I6)') '- length: ',(NBANDSDUMP)*WDES%ISPIN
            WRITE(7,'(A)')    '  type: State'
            WRITE(7,'(A,I6)') '- length: ',(NBANDSDUMP)*WDES%ISPIN
            WRITE(7,'(A)')    '  type: State'
            WRITE(7,'(A)')    'elements:'
            WRITE(7,'(A)')    '  type: IeeeBinaryFile'
            WRITE(7,'(A)')    'unit: 0.1917011272153577       # = sqrt(Eh/eV)'
            WRITE(7,'(A)')    'metaData:'
#ifdef gammareal
            WRITE(7,'(A)')    '  halfGrid: 1'
#else
            WRITE(7,'(A)')    '  halfGrid: 0'
#endif
            CLOSE(7)

         ENDIF


         IF (ME==0) OPEN(unit = 7,file = "CoulombVertex.elements",FORM='UNFORMATTED',access='stream',STATUS='REPLACE')

         DO SI=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)

            SPA=SP_ORD(SI)
            I=N_ORD(SI)

            IF (I <= NFREEZE) CYCLE
            CVERTEX_TMP=(0.0_q,0.0_q)

            DO J=1,(NBANDSDUMP+NFREEZE)
               IF (J <= NFREEZE) CYCLE
               IF (PRO_NI(J-NFREEZE)==ME) THEN
                  CALL GLOB2LOC(J-NFREEZE,DESC_OPTCVERTEX_1D(4),NPCOL_1D,DESC_OPTCVERTEX_1D(6),TMPMYCOL,MYJ)
                  CVERTEX_TMP(:,J-NFREEZE,SPA)=OPTCVERTEX_1D(:,MYJ,I-NFREEZE,1,1,SPA)
               ENDIF
            ENDDO
            CALL M_sum_z(WDES%COMM_INTER, CVERTEX_TMP(1,1,SPA),NOPTAUX*(NBANDSDUMP))

            IF (ME==0) THEN
               DO SJ=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)
                  SPB=SP_ORD(SJ)
                  J=N_ORD(SJ)
                  IF (J <= NFREEZE) CYCLE
                  CVERTEX_TMP_SINGLE(:)=(0.0_q,0.0_q)
                  IF (SPA==SPB) CVERTEX_TMP_SINGLE(:)=CVERTEX_TMP(:,J-NFREEZE,SPA)
                     WRITE(7) CVERTEX_TMP_SINGLE(:)
               ENDDO
            ENDIF

         ENDDO
         IF (ME==0) CLOSE(7)


         IF (ME==0) write(*,*)'- Writing EigenEnergies .'

         IF (ME==0)  THEN

            OPEN(unit = 7,file = "EigenEnergies.yaml",FORM='FORMATTED',access='stream',STATUS='REPLACE')
            WRITE(7,'(A)')    'version: 100'
            WRITE(7,'(A)')    'type: Tensor'
            WRITE(7,'(A)')    'scalarType: Real64'
            WRITE(7,'(A)')    'dimensions:'
            WRITE(7,'(A,I6)') '- length: ',((NBANDSDUMP)*WDES%ISPIN)
            WRITE(7,'(A)')    '  type: State'
            WRITE(7,'(A)')    'elements:'
            WRITE(7,'(A)')    '  type: TextFile'
            WRITE(7,'(A)')    'unit: 0.03674932217563878       # = (Eh/eV)'
            WRITE(7,'(A)')    'metaData:'
            WRITE(7,'(A,E22.15)')    '  fermiEnergy: ',EFERMI
            WRITE(7,'(A)')    '  energies:'

            DO SI=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)
               SPA=SP_ORD(SI)
               I=N_ORD(SI)
               IF (I <= NFREEZE) CYCLE
                IF (ME==0) WRITE(7,*) '  - ',REAL(W%CELTOT(I,1,SPA),kind=8)
            ENDDO

            CLOSE(7)

         ENDIF

         IF (ME==0)  THEN
            OPEN(unit = 7,file = "EigenEnergies.elements",FORM='FORMATTED',access='stream',STATUS='REPLACE')

            DO SI=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)
               SPA=SP_ORD(SI)
               I=N_ORD(SI)
               IF (I <= NFREEZE) CYCLE
                IF (ME==0) WRITE(7,*) REAL(W%CELTOT(I,1,SPA),kind=8)
            ENDDO


            CLOSE(7)
         ENDIF


         IF (ME==0) write(*,*)'- Writing States .'
!this needs to be adjusted for spin and k-points in future
         IF (ME==0)  THEN

            OPEN(unit = 7,file = "State.yaml",FORM='FORMATTED',access='stream',STATUS='REPLACE')
            WRITE(7,'(A)')    'version: 100'
            WRITE(7,'(A)')    'dimensionType: State'
            WRITE(7,'(A)')    'properties:'
            WRITE(7,'(A)')    '  CrystalMomentum:'
            WRITE(7,'(A)')    '  - - 0'
            WRITE(7,'(A)')    '    - 0'
            WRITE(7,'(A)')    '    - 0'
            WRITE(7,'(A)')    '  Spin:'
            WRITE(7,'(A)')    '  - 0'
            WRITE(7,'(A)')    'propertyIndices:'
            WRITE(7,'(A)')    '  CrystalMomentum:'
            DO SI=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)
               SPA=SP_ORD(SI)
               I=N_ORD(SI)
               IF (I <= NFREEZE) CYCLE
                           WRITE(7,'(A)')    '  - 0'
            ENDDO
            WRITE(7,'(A)')    '  Spin:'
            DO SI=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)
               SPA=SP_ORD(SI)
               I=N_ORD(SI)
               IF (I <= NFREEZE) CYCLE
                           WRITE(7,'(A)')    '  - 0'
            ENDDO

         ENDIF

         IF (ME==0) write(*,*)'- Writing AuxiliaryField .'
!this needs to be adjusted for k-points in future
         IF (ME==0)  THEN

            OPEN(unit = 7,file = "AuxiliaryField.yaml",FORM='FORMATTED',access='stream',STATUS='REPLACE')
            WRITE(7,'(A)')    'version: 100'
            WRITE(7,'(A)')    'dimensionType: AuxiliaryField'
            WRITE(7,'(A)')    'properties:'
            WRITE(7,'(A)')    '  CrystalMomentum:'
            WRITE(7,'(A)')    '  - - 0'
            WRITE(7,'(A)')    '    - 0'
            WRITE(7,'(A)')    '    - 0'
            WRITE(7,'(A)')    'propertyIndices:'
            WRITE(7,'(A)')    '  CrystalMomentum:'
            DO SI=1,NOPTAUX
                           WRITE(7,'(A)')    '  - 0'
            ENDDO

         ENDIF


      END SUBROUTINE WRITEBINCVERTEX2DISK_YAML

!***********************************************************************
!
! Write an input file needed by Cc4s to disk:
! CoulombVertexSingularVectors
!
!***********************************************************************

      SUBROUTINE WRITE_TSVEC_1D_GVEC_YAML(WDES)
         use mkpoints
         IMPLICIT NONE
         TYPE(wavedes) WDES
         integer :: TMPMYCOL,MYNAUX
         INTEGER :: NG,NAUX, NSTRIP, NSTRIPAUX, NNAUX
         COMPLEX(q), ALLOCATABLE :: TSVEC_TMP(:,:)

         NSTRIP=30
         ALLOCATE(TSVEC_TMP(NGVECTOR,NSTRIP))

         IF (ME==0) write(*,*)'- Writing CoulombVertexSingularVectors .'



         IF (ME==0)  THEN

            OPEN(unit = 7,file = "CoulombVertexSingularVectors.yaml",FORM='FORMATTED',access='stream',STATUS='REPLACE')
            WRITE(7,'(A)')    'version: 100'
            WRITE(7,'(A)')    'type: Tensor'
            WRITE(7,'(A)')    'scalarType: Complex64'
            WRITE(7,'(A)')    'dimensions:'
#ifdef gammareal
            WRITE(7,'(A,I6)') '  - length: ',NGVECTOR*2-1
#else
            WRITE(7,'(A,I6)') '  - length: ',NGVECTOR
#endif
            WRITE(7,'(A)')    '    type: Momentum'
            WRITE(7,'(A,I6)') '  - length: ',NOPTAUX
            WRITE(7,'(A)')    '    type: AuxiliaryField'
            WRITE(7,'(A)')    'elements:'
            WRITE(7,'(A)')    '  type: IeeeBinaryFile'
            WRITE(7,'(A)')    'unit: 1.0       #'
            CLOSE(7)

         ENDIF


         IF (ME==0) OPEN(unit = 7,file = "CoulombVertexSingularVectors.elements",FORM='UNFORMATTED',access='stream',STATUS='REPLACE')

         call BLACS_GRIDINFO(CONTXT_1D, NPROW_1D, NPCOL_1D, MYROW_1D, MYCOL_1D)

         DO NAUX=1,NOPTAUX,NSTRIP
  
            NSTRIPAUX=MIN(NOPTAUX + 1 - NAUX, NSTRIP)
            TSVEC_TMP(:,:)=(0.0_q,0.0_q)
            DO NNAUX=1,NSTRIPAUX

               CALL GLOB2LOC(NAUX+NNAUX-1,DESC_TSVEC_1D(4),NPCOL_1D,DESC_TSVEC_1D(6),TMPMYCOL,MYNAUX)
               IF (TMPMYCOL==MYCOL_1D) THEN 
                  TSVEC_TMP(1:NGVECTOR,NNAUX)=TSVEC_1D(1:NGVECTOR,MYNAUX)
                  CALL ZGEBS2D(CONTXT_1D, 'All', 'i-ring', NGVECTOR, 1, TSVEC_1D(1,MYNAUX), NGVECTOR)
               
               ELSE
                  CALL ZGEBR2D(CONTXT_1D, 'All', 'i-ring', NGVECTOR, 1, TSVEC_TMP(1,NNAUX), NGVECTOR, 0, TMPMYCOL)
               ENDIF

            ENDDO


            DO NNAUX=1,NSTRIPAUX

               IF (ME==0) THEN
                  WRITE(7) TSVEC_TMP(1,NNAUX)
               ENDIF


#ifdef gammareal
               TSVEC_TMP(2:NGVECTOR,NNAUX)=TSVEC_TMP(2:NGVECTOR,NNAUX)*sqrt(0.5_q)
               IF (ME==0) THEN
                  DO NG=2,NGVECTOR
                     WRITE(7) TSVEC_TMP(NG,NNAUX)
                  ENDDO 
               ENDIF
               IF (ME==0) THEN
                  DO NG=2,NGVECTOR
                     WRITE(7) CONJG(TSVEC_TMP(NG,NNAUX))
                  ENDDO 
               ENDIF
#else
               IF (ME==0) THEN
                  DO NG=2,NGVECTOR
                     WRITE(7) TSVEC_TMP(NG,NNAUX)
                  ENDDO 
               ENDIF
#endif

            ENDDO
         ENDDO

         IF (ME==0) CLOSE(7)

      END SUBROUTINE WRITE_TSVEC_1D_GVEC_YAML

!***********************************************************************
!
! Write several input files need by Cc4s to disk including:
! GridVectors
! CoulombPotential
! Momentum
!
!***********************************************************************

      SUBROUTINE WRITE_GVEC_POTFAK_YAML(WDES,LATT_CUR)
         use mkpoints
         USE constant
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(latt) LATT_CUR
         INTEGER :: NG,NI,KQ


         ! Write G-vector Mesh

         IF (ME==0)  THEN

            OPEN(unit = 7,file = "GridVectors.yaml",FORM='FORMATTED',access='stream',STATUS='REPLACE')
            WRITE(7,'(A)')    'version: 100'
            WRITE(7,'(A)')    'type: Tensor'
            WRITE(7,'(A)')    'scalarType: Real64'
            WRITE(7,'(A)')    'dimensions:'
            WRITE(7,'(A)')    '  - length: 3'
            WRITE(7,'(A)')    '    type: Vector'
#ifdef gammareal
            WRITE(7,'(A,I6)') '  - length: ',NGVECTOR*2-1
#else
            WRITE(7,'(A,I6)') '  - length: ',NGVECTOR
#endif
            WRITE(7,'(A)')    '    type: Momentum'
            WRITE(7,'(A)')    'elements:'
            WRITE(7,'(A)')    '  type: TextFile'
            WRITE(7,'(A)')    'unit: 0.529177249       # =(Bohr^-1/Angstrom^-1)'
            WRITE(7,'(A)')    'metaData:'
            WRITE(7,'("  Gi: [",E22.15,",",E22.15,",",E22.15,"]")')  TPI*LATT_CUR%B(1, 1),TPI*LATT_CUR%B(2,1),TPI*LATT_CUR%B(3, 1)
            WRITE(7,'("  Gj: [",E22.15,",",E22.15,",",E22.15,"]")')  TPI*LATT_CUR%B(1, 2),TPI*LATT_CUR%B(2,2),TPI*LATT_CUR%B(3, 2)
            WRITE(7,'("  Gk: [",E22.15,",",E22.15,",",E22.15,"]")')  TPI*LATT_CUR%B(1, 3),TPI*LATT_CUR%B(2,3),TPI*LATT_CUR%B(3, 3)
            CLOSE(7)

         ENDIF


         IF (ME==0) OPEN(unit = 7,file = "GridVectors.elements",FORM='FORMATTED',access='stream',STATUS='REPLACE')
         
         KQ=1 !k-points to be implemented

         DO NG=1,NGVECTOR
            DO NI=1,3
               IF (ME==0) THEN
                  WRITE(7,*) TPI*GVEC_FULL(NI,NG,KQ)
               ENDIF
            ENDDO
         ENDDO 
#ifdef gammareal
         DO NG=2,NGVECTOR
            DO NI=1,3
               IF (ME==0) THEN
                  WRITE(7,*) -TPI*GVEC_FULL(NI,NG,KQ)
               ENDIF
            ENDDO
         ENDDO 
#else
#endif

         IF (ME==0) CLOSE(7)

         IF (ME==0)  THEN

            OPEN(unit = 7,file = "CoulombPotential.yaml",FORM='FORMATTED',access='stream',STATUS='REPLACE')
            WRITE(7,'(A)')    'version: 100'
            WRITE(7,'(A)')    'type: Tensor'
            WRITE(7,'(A)')    'scalarType: Real64'
            WRITE(7,'(A)')    'dimensions:'
#ifdef gammareal
            WRITE(7,'(A,I6)') '  - length: ',NGVECTOR*2-1
#else
            WRITE(7,'(A,I6)') '  - length: ',NGVECTOR
#endif
            WRITE(7,'(A)')    '    type: Momentum'
            WRITE(7,'(A)')    'elements:'
            WRITE(7,'(A)')    '  type: TextFile'
            WRITE(7,'(A)')    'unit: 0.2479966649373453       # =(Eh/eV*Bohr^3/Angstrom^3)'
            CLOSE(7)

         ENDIF


         IF (ME==0) OPEN(unit = 7,file = "CoulombPotential.elements",FORM='FORMATTED',access='stream',STATUS='REPLACE')
         
         KQ=1 !k-points to be implemented

         DO NG=1,NGVECTOR
            IF (ME==0) THEN
               WRITE(7,*) REAL(POTFAK_FULL(NG,KQ)*CONJG(POTFAK_FULL(NG,KQ)),kind=q)
            ENDIF
         ENDDO 
#ifdef gammareal
         DO NG=2,NGVECTOR
            IF (ME==0) THEN
               WRITE(7,*) REAL(POTFAK_FULL(NG,KQ)*CONJG(POTFAK_FULL(NG,KQ)),kind=q)
            ENDIF
         ENDDO 
#else
#endif

         IF (ME==0) CLOSE(7)


         IF (ME==0) write(*,*)'- Writing Momentum .'
!this needs to be adjusted for k-points in future
         IF (ME==0)  THEN

            OPEN(unit = 7,file = "Momentum.yaml",FORM='FORMATTED',access='stream',STATUS='REPLACE')
            WRITE(7,'(A)')    'version: 100'
            WRITE(7,'(A)')    'dimensionType: Momentum'
            WRITE(7,'(A)')    'properties:'
            WRITE(7,'(A)')    '  CrystalMomentum:'
            WRITE(7,'(A)')    '  - - 0'
            WRITE(7,'(A)')    '    - 0'
            WRITE(7,'(A)')    '    - 0'
            WRITE(7,'(A)')    'propertyIndices:'
            WRITE(7,'(A)')    '  CrystalMomentum:'
#ifdef gammareal
            DO NG=1,NGVECTOR*2-1
                           WRITE(7,'(A)')    '  - 0'
            ENDDO
#else
            DO NG=1,NGVECTOR
                           WRITE(7,'(A)')    '  - 0'
            ENDDO
#endif

         ENDIF




      END SUBROUTINE WRITE_GVEC_POTFAK_YAML

!***********************************************************************
!
! Estimate truncation error of truncated optimized auxiliary field based
! using approximated second-order energy and EDIFF.
!
!***********************************************************************

      SUBROUTINE ESTIMATE_ACC_NOPTAUX(IO,W,WDES,INFO)
         use mkpoints
         IMPLICIT NONE
         TYPE (in_struct) IO
         type(wavespin) :: W
         TYPE(wavedes) WDES
         TYPE (info_struct) INFO
         !local variables
         INTEGER :: NG,NGP
         INTEGER :: I,J
         COMPLEX(q), ALLOCATABLE :: CVERTEX_TMP(:,:,:), CVERTEX_TMP_SINGLE(:), ENERGY(:),TEST_ENERGY(:)
         COMPLEX(q) :: vijab,vijba, VIIAA
         LOGICAL :: ex, CONSISTENT
         INTEGER :: nbstart, nastart
         integer :: TMPMYCOL,MYJ
         INTEGER :: SPI,SPJ,SI,SJ,KP,KQ,ISP,NP

         ALLOCATE(CVERTEX_TMP_SINGLE(NOPTAUX))
         ALLOCATE(ENERGY(NOPTAUX))
         ALLOCATE(TEST_ENERGY(NOPTAUX))
         ALLOCATE(CVERTEX_TMP(NOPTAUX,(NBANDSDUMP),WDES%ISPIN))

         CALL SETUP_TMPCVERTEX_1D(WDES)

         ENERGY=(0.0_q,0.0_q)
         TEST_ENERGY=(0.0_q,0.0_q)

         DO SI=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)

            SPI=SP_ORD(SI)
            I=N_ORD(SI)
            IF (W%FERTOT(I,1,SPI)<0.001_q) CYCLE


            IF (I <= NFREEZE) CYCLE

            CALL PZGEMR2D(NOPTAUX,(NBANDSDUMP),CVERTEX_2D(1,1,I-NFREEZE,1,1,SPI),1,1,&
              DESC_CVERTEX_2D,TMPCVERTEX_1D(1,1),1,1,DESC_TMPCVERTEX_1D,CONTXT_1D)


            CVERTEX_TMP=(0.0_q,0.0_q)

            DO J=1,(NBANDSDUMP+NFREEZE)
               IF (J <= NFREEZE) CYCLE
               IF (PRO_NI(J-NFREEZE)==ME) THEN
                  CALL GLOB2LOC(J-NFREEZE,DESC_TMPCVERTEX_1D(4),NPCOL_1D,DESC_TMPCVERTEX_1D(6),TMPMYCOL,MYJ)
                  CVERTEX_TMP(:,J-NFREEZE,SPI)=TMPCVERTEX_1D(:,MYJ)
               ENDIF
            ENDDO
            CALL M_sum_z(WDES%COMM_INTER, CVERTEX_TMP(1,1,SPI),NOPTAUX*(NBANDSDUMP))

               DO SJ=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)
                  SPJ=SP_ORD(SJ)
                  J=N_ORD(SJ)
                  IF (W%FERTOT(J,1,SPJ)>0.001_q) CYCLE
                  IF (J <= NFREEZE) CYCLE
                  CVERTEX_TMP_SINGLE(:)=(0.0_q,0.0_q)
                  IF (SPI==SPJ) CVERTEX_TMP_SINGLE(:)=CVERTEX_TMP(:,J-NFREEZE,SPI)
                  VIIAA=(0.0_q,0.0_q)
                  DO NG=1,NOPTAUX
#ifdef gammareal
                     VIIAA=VIIAA+(REAL(CVERTEX_TMP_SINGLE(NG),kind=q)*REAL(CVERTEX_TMP_SINGLE(NG),kind=q))
                     VIIAA=VIIAA+(DIMAG(CVERTEX_TMP_SINGLE(NG))*DIMAG(CVERTEX_TMP_SINGLE(NG)))
#else
                     VIIAA=VIIAA+(CONJG(CVERTEX_TMP_SINGLE(NG))*CVERTEX_TMP_SINGLE(NG))
#endif
                     ENERGY(NG)=ENERGY(NG)+2.0_q*VIIAA**2.0_q/(2.0_q*(W%CELTOT(I,1,SPI)-W%CELTOT(J,1,SPJ)))
                  ENDDO
               ENDDO

         ENDDO


         DO NG=NOPTAUX,1,-1
            IF (ABS(ENERGY(NG)-ENERGY(NOPTAUX))>INFO%EDIFF) THEN
               NOPTAUX_EDIFF=NG+1
               IF (ME==0) write(IO%IU6,*)'Truncating auxilliary field using EDIFF=',INFO%EDIFF,'results in ',NOPTAUX_EDIFF,' optimal auxiliary field vectors.'
               EXIT
            ENDIF
         ENDDO
         NOPTAUX=NOPTAUX_EDIFF
         CALL SETUP_OPTCVERTEX_1D(WDES)

         DO SPI=1,WDES%ISPIN
         DO KQ=1,REALNKPTS
         DO KP=1,REALNKPTS
         DO NP=1,NBANDSDUMP

            CALL PZGEMR2D(NGVECTOR,(NBANDSDUMP),CVERTEX_2D(1,1,NP,KP,KQ,SPI),1,1,&
              DESC_CVERTEX_2D,TMPCVERTEX_1D(1,1),1,1,DESC_TMPCVERTEX_1D,CONTXT_1D)

            OPTCVERTEX_1D(1:NOPTAUX_EDIFF,:,NP,KP,KQ,SPI)=TMPCVERTEX_1D(1:NOPTAUX_EDIFF,:)
         ENDDO
         ENDDO
         ENDDO
         ENDDO

         DEALLOCATE(CVERTEX_2D)

      END SUBROUTINE ESTIMATE_ACC_NOPTAUX


      SUBROUTINE SORT_ORBITALINDICES(W,WDES)
         use mkpoints
         IMPLICIT NONE
         type(wavespin) :: W
         TYPE(wavedes) WDES

         LOGICAL, ALLOCATABLE :: INCLUDED(:,:)
         REAL(q) :: CUR_EV, MAX_EV
         INTEGER :: SELECT_I, SELECT_SP,I,SP,SI


         ALLOCATE(N_ORD((NBANDSDUMP+NFREEZE)*WDES%ISPIN))
         ALLOCATE(SP_ORD((NBANDSDUMP+NFREEZE)*WDES%ISPIN))

         DO SI=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)
            N_ORD(SI)=SI !SELECT_I
            SP_ORD(SI)=1 !SELECT_SP
         ENDDO

      END SUBROUTINE SORT_ORBITALINDICES

      SUBROUTINE CALC_GVEC_FULL(WGW,WDES,LATT_CUR)
         use mkpoints
         USE constant
         IMPLICIT NONE
         TYPE (wavedes) WGW
         TYPE (wavedes) WDES
         TYPE (wavedes1) WDES1
         TYPE (latt) :: LATT_CUR
     ! local
         INTEGER    NI,NP,KQ
         REAL(q) :: DKX,DKY,DKZ,GX,GY,GZ,GSQU,GSQUP,GQUAD,SCALE, SCALEFSG, RHOEFF, FOMEGAP2

         DO KQ=1,WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0.00_q)) CYCLE

            CALL SETWDES(WGW,WDES1,KQ)

            DKX=WDES1%VKPT(1)*LATT_CUR%B(1,1)+ &
                WDES1%VKPT(2)*LATT_CUR%B(1,2)+ &
                WDES1%VKPT(3)*LATT_CUR%B(1,3)
            DKY=WDES1%VKPT(1)*LATT_CUR%B(2,1)+ &
                WDES1%VKPT(2)*LATT_CUR%B(2,2)+ &
                WDES1%VKPT(3)*LATT_CUR%B(2,3)
            DKZ=WDES1%VKPT(1)*LATT_CUR%B(3,1)+ &
                WDES1%VKPT(2)*LATT_CUR%B(3,2)+ &
                WDES1%VKPT(3)*LATT_CUR%B(3,3)

            NP=(WGW%NGVECTOR(KQ))

            SCALE=EDEPS/LATT_CUR%OMEGA/TPI**2

            DO NI=1,NP

               GX=(WDES1%IGX(NI)*LATT_CUR%B(1,1)+WDES1%IGY(NI)* &
                    LATT_CUR%B(1,2)+WDES1%IGZ(NI)*LATT_CUR%B(1,3))
               GY=(WDES1%IGX(NI)*LATT_CUR%B(2,1)+WDES1%IGY(NI)* &
                    LATT_CUR%B(2,2)+WDES1%IGZ(NI)*LATT_CUR%B(2,3))
               GZ=(WDES1%IGX(NI)*LATT_CUR%B(3,1)+WDES1%IGY(NI)* &
                    LATT_CUR%B(3,2)+WDES1%IGZ(NI)*LATT_CUR%B(3,3))
               GSQU=(DKX+GX)**2+(DKY+GY)**2+(DKZ+GZ)**2

               GVEC_FULL(1,NI,RKQofKQ(KQ))=DKX+GX
               GVEC_FULL(2,NI,RKQofKQ(KQ))=DKY+GY
               GVEC_FULL(3,NI,RKQofKQ(KQ))=DKZ+GZ

            ENDDO
         ENDDO

      END SUBROUTINE CALC_GVEC_FULL

!***********************************************************************
!
! Write input files needed by Cc4s to disk:
! DeltaIntegralsPPHH
! DeltaIntegralsHH
!
!***********************************************************************

      SUBROUTINE WRITE_D2PAW_YAML(W,WDES)
         use mkpoints
         IMPLICIT NONE
         type(wavespin) :: W
         TYPE(wavedes) WDES
         INTEGER :: NKI,NKJ,NKA,NI,NJ,NA,NKQ,NB,NKB,NG,NKK,NK
         INTEGER :: I,J,A,B,ISP1,ISP2,ISP3,ISP4,SP
         REAL(q) :: FA,FB
         COMPLEX(q), ALLOCATABLE :: CVERTEX_TMP(:,:,:), CVERTEX_TMP_SINGLE(:)
         COMPLEX(q) :: vijab,vijba, emp2
         LOGICAL :: ex, CONSISTENT
         INTEGER :: nbstart, nastart
         INTEGER*8 :: JUNKSIZE,IREC
         CHARACTER*8 :: MAGIC
         integer :: TMPMYCOL,MYJ
         INTEGER :: SPA,SPB,SI,SJ,NOCC,SA,SB,SPI,SPJ
         LOGICAL :: LISEND(WDES%ISPIN*(NBANDSDUMP+NFREEZE))
         INTEGER :: MNB, RMNB
         INTEGER :: LOCB(WDES%ISPIN*(NBANDSDUMP+NFREEZE)), PROCSEND(WDES%ISPIN*(NBANDSDUMP+NFREEZE))
         INTEGER :: tag=202, IERROR
         INTEGER :: NSTRIP, PNSTRIP, SIP, VBMAXP, NBANDSPROC(PROCS), PINDEX, RPINDEX
         GDEF, ALLOCATABLE :: TMPD2PAW(:,:)
         INTEGER :: status(MPI_STATUS_SIZE, PROCS), ROOT
         INTEGER :: req(PROCS), MYCOL

         ALLOCATE(TMPD2PAW(NBANDSDUMP,NBANDSDUMP))

         IF (ME==0)  THEN

            OPEN(unit = 7,file = "DeltaIntegralsPPHH.yaml",FORM='FORMATTED',access='stream',STATUS='REPLACE')
            WRITE(7,'(A)')    'version: 100'
            WRITE(7,'(A)')    'type: Tensor'
#ifdef gammareal
            WRITE(7,'(A)') 'scalarType: Real64'
#else
            WRITE(7,'(A)') 'scalarType: Complex64'
#endif
            WRITE(7,'(A)')    'dimensions:'
            WRITE(7,'(A,I6)') '  - length: ',((NBANDSDUMP-(VBMAX-NFREEZE))*WDES%ISPIN)
            WRITE(7,'(A)')    '    type: State'
            WRITE(7,'(A,I6)') '  - length: ',((NBANDSDUMP-(VBMAX-NFREEZE))*WDES%ISPIN)
            WRITE(7,'(A)')    '    type: State'
            WRITE(7,'(A,I6)') '  - length: ',((VBMAX-NFREEZE)*WDES%ISPIN)
            WRITE(7,'(A)')    '    type: State'
            WRITE(7,'(A,I6)') '  - length: ',((VBMAX-NFREEZE)*WDES%ISPIN)
            WRITE(7,'(A)')    '    type: State'
            WRITE(7,'(A)')    'elements:'
            WRITE(7,'(A)')    '  type: IeeeBinaryFile'
            WRITE(7,'(A)') 'unit: 0.1481847434769  # = (Angstrom^3/Bohr^3)'
            CLOSE(7)

         ENDIF


         IF (ME==0)  THEN
            OPEN(unit = 7,file = "DeltaIntegralsPPHH.elements",FORM='UNFORMATTED',access='stream',STATUS='REPLACE')
         ENDIF

         DO SJ=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)
            SPJ=SP_ORD(SJ)
            J=N_ORD(SJ)
            IF (J <= NFREEZE) CYCLE
            IF (J>VBMAX) CYCLE

         DO SI=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)
            SPI=SP_ORD(SI)
            I=N_ORD(SI)
            IF (I <= NFREEZE) CYCLE
            IF (I>VBMAX) CYCLE


            TMPD2PAW=zero
            DO RMNB=1,NBANDSDUMP+NFREEZE
               IF (RMNB<=NFREEZE) CYCLE
               CALL GLOB2LOC(RMNB,NBANDSDUMP+NFREEZE,NPCOL_1D,1,MYCOL,MNB)
               IF (MYCOL==MYCOL_1D) THEN
                  TMPD2PAW(:,RMNB-NFREEZE)=D2PAW_VVOO(:,MNB-NFREEZELOC,I-NFREEZE,J-NFREEZE,1,1,1)/0.1481847434769
               ENDIF
            ENDDO

            !The communication of the delta integrals can be made much more efficient
#ifdef gammareal
            CALLMPI( M_sum_d(WDES%COMM,TMPD2PAW(1,1),NBANDSDUMP*NBANDSDUMP) )
#else
            CALLMPI( M_sum_z(WDES%COMM,TMPD2PAW(1,1),NBANDSDUMP*NBANDSDUMP) )
#endif

            IF (ME==0) WRITE(7) TMPD2PAW((VBMAX-NFREEZE+1):,(VBMAX-NFREEZE+1):)

         ENDDO
         ENDDO

         IF (ME==0)  THEN
            CLOSE(7)
         ENDIF


         IF (ME==0)  THEN

            OPEN(unit = 7,file = "DeltaIntegralsHH.yaml",FORM='FORMATTED',access='stream',STATUS='REPLACE')
            WRITE(7,'(A)')    'version: 100'
            WRITE(7,'(A)')    'type: Tensor'
#ifdef gammareal
            WRITE(7,'(A)') 'scalarType: Real64'
#else
            WRITE(7,'(A)') 'scalarType: Complex64'
#endif
            WRITE(7,'(A)')    'dimensions:'
            WRITE(7,'(A,I6)') '  - length: ',((VBMAX-NFREEZE)*WDES%ISPIN)
            WRITE(7,'(A)')    '    type: State'
            WRITE(7,'(A,I6)') '  - length: ',((VBMAX-NFREEZE)*WDES%ISPIN)
            WRITE(7,'(A)')    '    type: State'
            WRITE(7,'(A)')    'elements:'
            WRITE(7,'(A)')    '  type: IeeeBinaryFile'
            WRITE(7,'(A)') 'unit: 0.1481847434769  # = (Angstrom^3/Bohr^3)'
            CLOSE(7)

         ENDIF


         IF (ME==0)  THEN
            OPEN(unit = 7,file = "DeltaIntegralsHH.elements",FORM='UNFORMATTED',access='stream',STATUS='REPLACE')
         ENDIF

         DO SJ=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)
            SPJ=SP_ORD(SJ)
            J=N_ORD(SJ)
            IF (J <= NFREEZE) CYCLE

            DO SI=1,WDES%ISPIN*(NBANDSDUMP+NFREEZE)
               I=N_ORD(SI)
               SPI=SP_ORD(SI)
               IF (I <= NFREEZE) CYCLE

               IF ((I<=VBMAX) .and. (J<=VBMAX)) THEN

                  TMPD2PAW=zero
                  DO RMNB=1,NBANDSDUMP+NFREEZE
                     IF (RMNB<=NFREEZE) CYCLE
                     CALL GLOB2LOC(RMNB,NBANDSDUMP+NFREEZE,NPCOL_1D,1,MYCOL,MNB)
                     IF (MYCOL==MYCOL_1D) THEN
                        TMPD2PAW(:,RMNB-NFREEZE)=D2PAW_VVOO(:,MNB-NFREEZELOC,I-NFREEZE,J-NFREEZE,1,1,1)/0.1481847434769
                     ENDIF
                  ENDDO


            !The communication of the delta integrals can be made much more efficient
#ifdef gammareal
                  CALLMPI( M_sum_d(WDES%COMM,TMPD2PAW(1,1),NBANDSDUMP*NBANDSDUMP))
#else
                  CALLMPI( M_sum_z(WDES%COMM,TMPD2PAW(1,1),NBANDSDUMP*NBANDSDUMP))
#endif


                  IF (ME==0) THEN
                     IF (SPI==SPJ) THEN
#ifdef gammareal
                        WRITE(7) TMPD2PAW(I-NFREEZE,J-NFREEZE)
#else
                        WRITE(7) TMPD2PAW(I-NFREEZE,J-NFREEZE)
#endif
                     ELSE
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO

         IF (ME==0)  THEN
            CLOSE(7)
         ENDIF


      END SUBROUTINE WRITE_D2PAW_YAML

!***********************************************************************
!
! Compute two-electron delta integrals in the PAW framework
!
!***********************************************************************

   SUBROUTINE CALC_D2PAW(LATT_CUR, W, WDES, T_INFO, P, IO)
      USE constant
      USE wave
      USE pseudo
      USE lattice
      USE full_kpoints
      USE asa
      USE radial
      IMPLICIT NONE
      TYPE(latt) LATT_CUR
      TYPE(wavespin) W
      TYPE(wavedes) WDES
      TYPE(type_info) T_INFO
      TYPE(potcar) P(T_INFO%NTYP)
      TYPE(in_struct) IO
! local variables
      INTEGER :: NG, KQ, NI, NJ, KI, KJ, NA, NB, KA, KB, ISP, KI_IN_FULL_ORIG
! for D2PAW calc
      INTEGER :: NSTRIP, I
      COMPLEX(q) :: INTC
      TYPE(wavespin) WHF
      TYPE(wavedes1) WDESKI, WDESKJ, WDESKA, WDESKB
      TYPE(wavefun1), ALLOCATABLE :: WI(:), WJ(:), WA(:), WB(:)
! for D2PAW one-center terms
      INTEGER :: CH1, CH2, CH3, CH4, ION, LMMAX, SPECIES, NBI, NBJ, NT
      INTEGER :: POSCH1, POSCH2, POSCH3, POSCH4
      INTEGER :: L0P, L1P, L2P, L3P, M0P, M1P, M2P, M3P, M0PMAX, M1PMAX, M2PMAX, M3PMAX
      COMPLEX(q), ALLOCATABLE :: DELTA_PAW(:,:,:,:,:) ! delta-kernel in PAW channel basis for each ion type 
      COMPLEX(q), ALLOCATABLE :: DELTA_PAW_RAD(:,:,:,:,:) ! delta-kernel in PAW channel basis for radial factor contribition
      COMPLEX(q), ALLOCATABLE :: DELTA_PAW_SPH(:,:,:,:,:) ! delta-kernel in PAW channel basis for angular factor contribition
      REAL(q) :: QRHO
      REAL(q), ALLOCATABLE :: RHO(:)
!For computing YLM3
      INTEGER :: LM01INDX, ISTART01, IEND01, LM23INDX, ISTART23, IEND23, IC01, IC23, NPRO
!  D2PAW_VVOO related variables
      INTEGER :: NSTRIPA, NSTRIPB, NBA, NBB, NBAP, NBBP
      GDEF, ALLOCATABLE :: D2PAW_INT(:,:,:,:) ! intermediate quantity for fast calculation of one-center terms
! intermediate quantities for enhanced BLAS level 3 performance 
      GDEF, ALLOCATABLE :: MAT1(:,:), MAT2(:,:), MAT3(:,:)
      INTEGER :: NTMP1, NTMP2

      ISP=1

!***********************************************************************
! start with D2PAW_OOOO
!***********************************************************************
      IF (IO%IU0 >= 0) WRITE (IO%IU6, *) 'Preparing calculation of delta-kernel integrals.'


!***********************************************************************
! compute four-channel--one-center--delta-kernel for each atomic species
! THIS COULD GO IN A SEPARATE ROUTINE
!***********************************************************************
     !find lmmax of all ions
     LMMAX=0
     DO SPECIES=1,T_INFO%NTYP
        IF (P(SPECIES)%LMMAX > LMMAX) LMMAX=P(SPECIES)%LMMAX
     ENDDO

     !allocate and compute DELTA_PAW
     ALLOCATE(DELTA_PAW(LMMAX,LMMAX,LMMAX,LMMAX,T_INFO%NTYP))
     ALLOCATE(DELTA_PAW_RAD(LMMAX,LMMAX,LMMAX,LMMAX,T_INFO%NTYP))
     ALLOCATE(DELTA_PAW_SPH(LMMAX,LMMAX,LMMAX,LMMAX,T_INFO%NTYP))
     ALLOCATE(D2PAW_INT(LMMAX,LMMAX,W%WDES%NIONS,VBMAX-NFREEZE))

     DELTA_PAW=zero
     DELTA_PAW_RAD=zero
     DELTA_PAW_SPH=zero

     DO SPECIES=1,T_INFO%NTYP

        IF (P(SPECIES)%LMAX==0) EXIT

     !for each ion type compute delta(r1-r2) kernel integral of pseudo and ae partial wave contributions
     !
     !this algorithm is not efficient but will not be a bottle neck

        IF (ALLOCATED(RHO)) DEALLOCATE(RHO)
        ALLOCATE(RHO(P(SPECIES)%R%NMAX))
        RHO=0._q

        POSCH1=1
        DO CH1=1,P(SPECIES)%LMAX
        !get L0P and M0P
        L0P =P(SPECIES)%LPS(CH1)

        POSCH2=1
        DO CH2=1,P(SPECIES)%LMAX
        !get L1P and M1P
        L1P =P(SPECIES)%LPS(CH2)

        POSCH3=1
        DO CH3=1,P(SPECIES)%LMAX
        !get L2P and M2P
        L2P =P(SPECIES)%LPS(CH3)

        POSCH4=1
        DO CH4=1,P(SPECIES)%LMAX
           !get L3P and M3P
           L3P =P(SPECIES)%LPS(CH4)

           !compute factor of delta-kernel contribution coming from radial coordinate integration 
           !
           !(some work is duplicated for all values of m,m',m'' and m''' quantum numbers).

           !N.B. we have to take care of the 4*pi/r^2 factor because the AE and
           !PS radial wavefunctions are transformed with the corresponding
           !pre-factors
           RHO=0._q
           QRHO=0._q
           DO I=1,P(SPECIES)%R%NMAX
              RHO(I)=RHO(I)+(P(SPECIES)%WAE(I,CH1)*P(SPECIES)%WAE(I,CH2)*P(SPECIES)%WAE(I,CH3)*P(SPECIES)%WAE(I,CH4))/P(SPECIES)%R%R(I)**2*4*PI
              RHO(I)=RHO(I)-(P(SPECIES)%WAE(I,CH1)*P(SPECIES)%WAE(I,CH2)*P(SPECIES)%WPS(I,CH3)*P(SPECIES)%WPS(I,CH4))/P(SPECIES)%R%R(I)**2*4*PI
              RHO(I)=RHO(I)-(P(SPECIES)%WPS(I,CH1)*P(SPECIES)%WPS(I,CH2)*P(SPECIES)%WAE(I,CH3)*P(SPECIES)%WAE(I,CH4))/P(SPECIES)%R%R(I)**2*4*PI
              RHO(I)=RHO(I)+(P(SPECIES)%WPS(I,CH1)*P(SPECIES)%WPS(I,CH2)*P(SPECIES)%WPS(I,CH3)*P(SPECIES)%WPS(I,CH4))/P(SPECIES)%R%R(I)**2*4*PI

           ENDDO
           CALL SIMPI(P(SPECIES)%R,RHO,QRHO)

           M0PMAX=2*L0P+1
           M1PMAX=2*L1P+1
           M2PMAX=2*L2P+1
           M3PMAX=2*L3P+1

           DO M0P=0,M0PMAX-1
           DO M1P=0,M1PMAX-1
           DO M2P=0,M2PMAX-1
           DO M3P=0,M3PMAX-1
              DELTA_PAW_RAD(POSCH1+M0P,POSCH2+M1P,POSCH3+M2P,POSCH4+M3P,SPECIES)=QRHO
           ENDDO
           ENDDO
           ENDDO
           ENDDO

           !next we need to compute DELTA_PAW_SPH =  <l,m.l',m'| L,M><L,M | l'',m'',l''',m'''>
           !
           ! For this we ńeed the  YLM3 Coefficients from the asa module
           ! The YLM3 coeffs correspond to the integrals of three spherical harmonics using l,m,l',m' and L,M
           
           ! First look up <l,m.l',m'| L,M>
           !   ... then  <L,M | l'',m'',l''',m'''>
           ! Then contract over L,M and store in DELTA_PAW_SPH

           CALL YLM3LOOKUP(L0P,L1P,LM01INDX)

           DO M0P=0,M0PMAX-1
           DO M1P=0,M1PMAX-1
              LM01INDX=LM01INDX+1

              ISTART01=INDCG(LM01INDX)
              IEND01  =INDCG(LM01INDX+1)

              CALL YLM3LOOKUP(L2P,L3P,LM23INDX)

              DO M2P=0,M2PMAX-1
              DO M3P=0,M3PMAX-1
 
                 LM23INDX=LM23INDX+1

                 ISTART23=INDCG(LM23INDX)
                 IEND23  =INDCG(LM23INDX+1)

                 DO IC01=ISTART01,IEND01-1
                 DO IC23=ISTART23,IEND23-1

                    IF (JL(IC01)==JL(IC23)) THEN
                    IF (JS(IC01)==JS(IC23)) THEN
                       DELTA_PAW_SPH(POSCH1+M0P,POSCH2+M1P,POSCH3+M2P,POSCH4+M3P,SPECIES)= &
                          DELTA_PAW_SPH(POSCH1+M0P,POSCH2+M1P,POSCH3+M2P,POSCH4+M3P,SPECIES) + YLM3(IC01)*YLM3(IC23)
                            
                    ENDIF
                    ENDIF
                 ENDDO
                 ENDDO


              ENDDO
              ENDDO
           ENDDO
           ENDDO

           DO M0P=0,M0PMAX-1
           DO M1P=0,M1PMAX-1
           DO M2P=0,M2PMAX-1
           DO M3P=0,M3PMAX-1
              DELTA_PAW(POSCH1+M0P,POSCH2+M1P,POSCH3+M2P,POSCH4+M3P,SPECIES)=&
                      DELTA_PAW_RAD(POSCH1+M0P,POSCH2+M1P,POSCH3+M2P,POSCH4+M3P,SPECIES) &
                     *DELTA_PAW_SPH(POSCH1+M0P,POSCH2+M1P,POSCH3+M2P,POSCH4+M3P,SPECIES)
           ENDDO
           ENDDO
           ENDDO
           ENDDO




        POSCH4=POSCH4+2*L3P+1
        ENDDO !ch4
        POSCH3=POSCH3+2*L2P+1
        ENDDO !ch3
        POSCH2=POSCH2+2*L1P+1
        ENDDO !ch2
        POSCH1=POSCH1+2*L0P+1
        ENDDO !ch1
     

     ENDDO

     IF (IO%IU0 >= 0) WRITE (IO%IU6, *) 'Transformation of Delta-kernel to one-center basis done.'

! the following work only needs to be done if NBANDSDUMPLOC>0

     ALLOCATE(D2PAW_VVOO(NBANDSDUMP,NBANDSDUMPLOC,VBMAX-NFREEZE,VBMAX-NFREEZE,REALNKPTS,REALNKPTS,REALNKPTS))

!***********************************************************************
!Add pseudo part contribution on real space grid   < ij | delta(r1-r2) | ab >
!and delta-kernel-contributions of PAW terms from all atoms to D2PAW_VVOO for each state
!***********************************************************************
      IF (IO%IU0 >= 0) WRITE (IO%IU6, *) '* Computing DeltaIntegralsPPHH'
      D2PAW_VVOO= zero

      NSTRIP=30

      WHF = W
      WHF%WDES => WDES_FOCK
      CALL SETWDES(WHF%WDES, WDESKI, 0)
      CALL SETWDES(WHF%WDES, WDESKJ, 0)
      CALL SETWDES(WHF%WDES, WDESKA, 0)
      CALL SETWDES(WHF%WDES, WDESKB, 0)


      IF (ALLOCATED(WI)) DEALLOCATE(WI)
      IF (ALLOCATED(WJ)) DEALLOCATE(WJ)
      IF (ALLOCATED(WA)) DEALLOCATE(WA)
      IF (ALLOCATED(WB)) DEALLOCATE(WB)

      ALLOCATE (WI(VBMAX-NFREEZE) )
      ALLOCATE (WJ(VBMAX-NFREEZE) )
      ALLOCATE (WA(NSTRIP) )
      ALLOCATE (WB(NSTRIP) )

      DO NBI = 1, VBMAX-NFREEZE
         CALL NEWWAV(WI(NBI), WDESKI, .TRUE.)
         CALL NEWWAV(WJ(NBI), WDESKJ, .TRUE.)
      END DO
      DO NBI = 1, NSTRIP
         CALL NEWWAV(WA(NBI), WDESKA, .TRUE.)
         CALL NEWWAV(WB(NBI), WDESKB, .TRUE.)
      END DO



      DO KI = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE

         KI_IN_FULL_ORIG = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:, KI), KPOINTS_FULL_ORIG)
         CALL SETWDES(WHF%WDES, WDESKI, KI)
         CALL W1_GATHER_GLB(WHF, NFREEZE+1, VBMAX, ISP, WI)



      DO KJ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE

         CALL SETWDES(WHF%WDES, WDESKJ, KJ)
         CALL W1_GATHER_GLB(WHF, NFREEZE+1, VBMAX, ISP, WJ)

      DO KA = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA) == 0)) CYCLE
         CALL SETWDES(WHF%WDES, WDESKA, KA)
         KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KA) - &
                                     WDES%VKPT(:, KI), KPOINTS_FULL)

         KB = KPOINT_IN_FULL_GRID(-WDES%VKPT(:, KQ) + WDES%VKPT(:, KJ), KPOINTS_FULL)

         CALL SETWDES(WHF%WDES, WDESKB, KB)

      DO NBJ = 1, VBMAX-NFREEZE

         IF (IO%IU0 >= 0) WRITE (IO%IU6, *) 'NBJ=',NBJ
         !CYCLE
         D2PAW_INT= zero

         DO NBI = 1, VBMAX-NFREEZE


         ! loop over atoms
         ! select corresponding species and loop over channels
            DO NI=1,W%WDES%NIONS
               SPECIES=W%WDES%ITYP(NI)
               LMMAX=W%WDES%LMMAX(SPECIES)
               IF (LMMAX==0) CYCLE
               NPRO =W%WDES%LMBASE(NI)

               DO CH1=1,LMMAX
               DO CH2=1,LMMAX
               DO CH3=1,LMMAX
               DO CH4=1,LMMAX

                  D2PAW_INT(CH3,CH4,NI,NBI)= &
                      D2PAW_INT(CH3,CH4,NI,NBI) + &
                      GCONJG(WI(NBI)%CPROJ(CH1+NPRO)*WJ(NBJ)%CPROJ(CH2+NPRO))* &
                      DELTA_PAW(CH1,CH2,CH3,CH4,SPECIES)



               ENDDO !ch4
               ENDDO !ch3
               ENDDO !ch2
               ENDDO !ch1
 
            ENDDO !NI

         ENDDO !NBI
!         IF (IO%IU0 >= 0) WRITE (*, *) 'D2PAW_INT built for NBJ=',NBJ

      ! loop over all bands
         DO NBA=1,NBANDSDUMP,NSTRIP
            NSTRIPA = MIN(NBANDSDUMP - NBA +1, NSTRIP)
            CALL W1_GATHER_GLB(WHF, NBA+NFREEZE, NFREEZE+NBA+NSTRIPA-1, ISP, WA)

         DO NBB=1,NBANDSDUMPLOC,NSTRIP
            NSTRIPB=MIN(NBANDSDUMPLOC+1-NBB,NSTRIP)
            ! FFT{psi_b} to real space
            DO NBBP=1,NSTRIPB !copy and fourier transform NSTRIP wave functions
!               WRITE(*,*)'ME,NFREEZELOC+NBB+NBBP-1',ME,NFREEZELOC+NBB+NBBP-1
               CALL W1_COPY( ELEMENT(WHF,WDESKB,NFREEZELOC+NBB+NBBP-1,ISP),WB(NBBP))
               CALL FFTWAV_W1(WB(NBBP))
            ENDDO


! allocate intermediate quantities for enhanced performance
            NTMP2=WI(VBMAX-NFREEZE)%WDES1%GRID%MPLWV
            IF (ALLOCATED(MAT1)) DEALLOCATE(MAT1)
            IF (ALLOCATED(MAT2)) DEALLOCATE(MAT2)
            IF (ALLOCATED(MAT3)) DEALLOCATE(MAT3)
            ALLOCATE(MAT1(NTMP2,NSTRIPA*NSTRIPB))
            ALLOCATE(MAT2(NTMP2,VBMAX-NFREEZE))
            ALLOCATE(MAT3(NSTRIPA*NSTRIPB,VBMAX-NFREEZE))

            NTMP1=0
            DO NBBP = 1, NSTRIPB
            DO NBAP = 1, NSTRIPA
               NTMP1=NTMP1+1
               DO I=1,NTMP2
                  MAT1(I,NTMP1)=WA(NBAP)%CR(I)*WB(NBBP)%CR(I)
               ENDDO
            ENDDO
            ENDDO
            DO NBI = 1, VBMAX-NFREEZE
               DO I=1,NTMP2
                  MAT2(I,NBI)=GCONJG(WI(NBI)%CR(I)*WJ(NBJ)%CR(I))/NTMP2
               ENDDO
            ENDDO !NBI

#ifdef gammareal
            CALL DGEMM('t', 'n', NSTRIPA*NSTRIPB, VBMAX-NFREEZE, NTMP2,one, &
                          MAT1(1, 1), NTMP2, &
                          MAT2(1, 1), NTMP2, &
                          zero, MAT3(1,1), NSTRIPA*NSTRIPB) 
#else
            CALL ZGEMM('t', 'n', NSTRIPA*NSTRIPB, VBMAX-NFREEZE, NTMP2,one, &
                          MAT1(1, 1), NTMP2, &
                          MAT2(1, 1), NTMP2, &
                          zero, MAT3(1,1), NSTRIPA*NSTRIPB) 
#endif
            DO NBI = 1, VBMAX-NFREEZE
               NTMP1=0
               DO NBBP = 1, NSTRIPB
               DO NBAP = 1, NSTRIPA
                  NTMP1=NTMP1+1
                  D2PAW_VVOO(NBAP + NBA - 1,NBBP + NBB - 1,NBI,NBJ,RKIofKI(KI), RKIofKI(KJ),KA)=MAT3(NTMP1,NBI)
               ENDDO
               ENDDO
            ENDDO !NBI

! slow implementation of above term without blas level 3


!            DO NBI = 1, VBMAX-NFREEZE

!               DO NBAP = 1, NSTRIPA
!               DO NBBP = 1, NSTRIPB


!                  DO NI=1,W%WDES%NIONS
!                     SPECIES=W%WDES%ITYP(NI)
!                     LMMAX=W%WDES%LMMAX(SPECIES)
!                     IF (LMMAX==0) CYCLE
!                     NPRO =W%WDES%LMBASE(NI)

!                     DO CH3=1,LMMAX
!                     DO CH4=1,LMMAX

!                        D2PAW_VVOO(NBA+NBAP-1,NBB+NBBP-1,NBI,NBJ,RKIofKI(KI), RKIofKI(KJ),RKIofKI(KA))= &
!                           D2PAW_VVOO(NBA+NBAP-1,NBB+NBBP-1,NBI,NBJ,RKIofKI(KI),RKIofKI(KJ),RKIofKI(KA)) + &
!                           WA(NBAP)%CPROJ(CH3+NPRO)*WB(NBBP)%CPROJ(CH4+NPRO) * &
!                           D2PAW_INT(CH3,CH4,NI,NBI)

!                     ENDDO !ch4
!                     ENDDO !ch3
    
!                  ENDDO !NI

!               ENDDO !NBBP
!               ENDDO !NBAP

!            ENDDO !NBI

! allocate intermediate quantities for enhanced performance
           IF (ALLOCATED(MAT3)) DEALLOCATE(MAT3)
           ALLOCATE(MAT3(NSTRIPA*NSTRIPB,VBMAX-NFREEZE))
           MAT3=zero

           DO NI=1,W%WDES%NIONS
              SPECIES=W%WDES%ITYP(NI)
              LMMAX=W%WDES%LMMAX(SPECIES)
              IF (LMMAX==0) CYCLE
              NPRO =W%WDES%LMBASE(NI)

              IF (ALLOCATED(MAT1)) DEALLOCATE(MAT1)
              IF (ALLOCATED(MAT2)) DEALLOCATE(MAT2)
              ALLOCATE(MAT1(LMMAX**2,NSTRIPA*NSTRIPB))
              ALLOCATE(MAT2(LMMAX**2,VBMAX-NFREEZE))

 
              NTMP2=0
              DO NBBP = 1, NSTRIPB
              DO NBAP = 1, NSTRIPA
                 NTMP2=NTMP2+1
                 NTMP1=0
                 DO CH3=1,LMMAX
                 DO CH4=1,LMMAX
                    NTMP1=NTMP1+1
                    MAT1(NTMP1,NTMP2)=WA(NBAP)%CPROJ(CH3+NPRO)*WB(NBBP)%CPROJ(CH4+NPRO)
                 ENDDO
                 ENDDO
              ENDDO
              ENDDO
 
              DO NBI = 1, VBMAX-NFREEZE
                 NTMP1=0
                 DO CH3=1,LMMAX
                 DO CH4=1,LMMAX
                    NTMP1=NTMP1+1
                    MAT2(NTMP1,NBI)=D2PAW_INT(CH3,CH4,NI,NBI)
                 ENDDO
                 ENDDO
              ENDDO

! the power of DGEMM comes here
#ifdef gammareal
              CALL DGEMM('t', 'n', NSTRIPA*NSTRIPB, VBMAX-NFREEZE, LMMAX**2,one, &
                          MAT1(1, 1), LMMAX**2, &
                          MAT2(1, 1), LMMAX**2, &
                          one, MAT3(1,1), NSTRIPA*NSTRIPB) 
#else    
              CALL ZGEMM('t', 'n', NSTRIPA*NSTRIPB, VBMAX-NFREEZE, LMMAX**2,one, &
                          MAT1(1, 1), LMMAX**2, &
                          MAT2(1, 1), LMMAX**2, &
                          one, MAT3(1,1), NSTRIPA*NSTRIPB) 
#endif
           ENDDO !NI atoms

           DO NBI = 1, VBMAX-NFREEZE
              NTMP1=0
              DO NBBP = 1, NSTRIPB
              DO NBAP = 1, NSTRIPA
                 NTMP1=NTMP1+1
                 D2PAW_VVOO(NBAP + NBA - 1,NBBP + NBB - 1,NBI,NBJ,RKIofKI(KI), RKIofKI(KJ),KA)= &
                     D2PAW_VVOO(NBAP+NBA-1,NBBP+NBB-1,NBI,NBJ,RKIofKI(KI), RKIofKI(KJ),KA) + MAT3(NTMP1,NBI)
              ENDDO
              ENDDO
           ENDDO !NBI



        ENDDO !NBA
        ENDDO !NBB

     ENDDO !NBJ


     ENDDO !ka loop
     ENDDO !kj loop
     ENDDO !ki loop

   END SUBROUTINE


    SUBROUTINE CHECK_SHIFTED_KPOINTS(WDES,W,KPOINTS)
         USE constant
         USE full_kpoints
         USE mkpoints
         USE wave
         implicit NONE
         TYPE(wavedes) WDES
         TYPE(wavespin) W
         TYPE (kpoints_struct) KPOINTS
         integer :: MKI,KI
         LOGICAL :: GAMMA_FOUND
         LOGICAL :: SECOND_GAMMA
         REAL(q) :: SECOND_GAMMA_WEIGHT
         REAL(q) :: GAMMA_WEIGHT
         INTEGER :: NRKI, NRKQ
         REAL(q) :: TINY=1.D-8

         GAMMA_FOUND=.FALSE.
         SHIFTED_KPOINTS=.TRUE.
         SECOND_GAMMA=.FALSE.
         do MKI=1,WDES%NKPTS
            if ( (ABS(W%WDES%VKPT(1,MKI))<(TINY)) .and. (ABS(W%WDES%VKPT(2,MKI))<(TINY)) .and. (ABS(W%WDES%VKPT(3,MKI))<(TINY))) then
               IF (GAMMA_FOUND) SECOND_GAMMA=.TRUE.
               IF (SECOND_GAMMA) SECOND_GAMMA_WEIGHT=W%WDES%WTKPT(MKI)
               IF (.NOT. GAMMA_FOUND) THEN
                  GAMMA_WEIGHT=W%WDES%WTKPT(MKI)
                  GAMMA_FOUND=.true.
               ENDIF
            endif
         enddo
         if ((GAMMA_FOUND) .and. (GAMMA_WEIGHT/=0.0_q)) THEN
            SHIFTED_KPOINTS=.FALSE.
         endif
         IF (SECOND_GAMMA) SHIFTED_KPOINTS=.TRUE.

         REALNKPTS=WDES%NKPTS
         IF (SHIFTED_KPOINTS) REALNKPTS=WDES%NKPTS/2

         ALLOCATE(RKQofKQ(WDES%NKPTS))
         ALLOCATE(RKIofKI(WDES%NKPTS))
         ALLOCATE(KQofRKQ(REALNKPTS))
         ALLOCATE(KIofRKI(REALNKPTS))
         NRKI=0
         NRKQ=0
         RKQofKQ=0
         RKIofKI=0
         KQofRKQ=0
         KIofRKI=0
         IF (SHIFTED_KPOINTS) THEN
            DO KI=1,WDES%NKPTS
               IF ((WDES%WTKPT(KI)/=0.00_q))  THEN
                  NRKI=NRKI+1
                  RKIofKI(KI)=NRKI
                  KIofRKI(NRKI)=KI
               ENDIF
               IF ((WDES%WTKPT(KI)==0.00_q))  THEN
                  NRKQ=NRKQ+1
                  RKQofKQ(KI)=NRKQ
                  KQofRKQ(NRKQ)=KI
               ENDIF

            ENDDO
         ELSE
            DO KI=1,WDES%NKPTS
               RKIofKI(KI)=KI
               KIofRKI(KI)=KI
               RKQofKQ(KI)=KI
               KQofRKQ(KI)=KI
            ENDDO
         ENDIF

        WTKPT=1.0_q/REALNKPTS

    END SUBROUTINE CHECK_SHIFTED_KPOINTS

!***********************************************************************
! Read some usefule control flags from the INCAR file
!***********************************************************************

    SUBROUTINE CC4S_INCAR_READER(IU5,IU6,IU0)
      USE base
      USE reader_tags
      INTEGER IU5,IU6,IU0
      INTEGER NTYP
      ! local variables
      INTEGER IDUM,N,IERR
      REAL(q) RDUM
      COMPLEX(q) CDUM
      LOGICAL LOPEN,LDUM,LRPA
      CHARACTER (1) CHARAC
      CHARACTER (40) SZNAM
      CHARACTER (40)   STRING
      INTEGER, EXTERNAL :: LENGTH


      LOPEN=.TRUE.

      NFREEZE=0
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NFREEZE', NFREEZE, IERR, WRITEXMLINCAR)
      
      LLFTODDUMP=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LFTODDUMP', LLFTODDUMP, IERR, WRITEXMLINCAR)

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBANDSHIGH', NBANDSDUMP, IERR, WRITEXMLINCAR)

      LSKIPFOCK=.TRUE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSKIPFOCK', LSKIPFOCK, IERR, WRITEXMLINCAR)


      RETURN

    END SUBROUTINE CC4S_INCAR_READER

#endif // scaLAPACK

END MODULE CC4S
EOF
    echo "Adding file  cc4s.F " 
cp .tmp_vaspcc4s_patch  cc4s.F
 
#######################
# new file   ccsd.F
#######################
cat > .tmp_vaspcc4s_patch<<"EOF"
#include "symbol.inc"
!*********************************************************************************
!
! MODULE :: CCSD
! AUTHOR :: Andreas Grüneis ( andreas.grueneis@tuwien.ac.at )
!
! DESCRIPTION: 
!
! THIS MODULE PERFORMS COUPLED CLUSTER THEORY CALCULATIONS AND IS PART OF VASP
! 
! ================================================
! THE FOLLOWING SETTINGS ARE SUPPORTED:
! ================================================
!
!   * RESTRICTED (WDES%ISPIN=1) CANONICAL ORBITALS ONLY
!   * GAMMA-CENTERED OR SHIFTED K-MESHES
!   * FOR THE SPECIAL CASE OF GAMMA-CENTERED 1X1X1 AND 2X2X2 K-MESHES REAL ORBITALS (LORBITALREAL=.TRUE.)
!      CAN BE EMPLOYED. THIS REDUCES THE MEMORY FOOTPRINT SIGNIFICANTLY.
!
!
! ================================================
!  STRUCTURE OF THE MAIN ROUTINES IN THIS  MODULE
! ================================================
!   
!                    VASP
!                     :
!            CALCULATE_CCSD (main routine that iterates until amplitudes are converged or max_it reached)
!                     |      
!                     |      
!          ----------------------------------------------------------------
!          |                     |     |             |     |              |
!   CALC_T1_T2_N -> (T1_N,T2_N)  |     |             |     |              |
!          |                     |     |             |     |         various disk I/O routines
!          |                     |     |             |     |    to read amplitudes, control flags
!          |                     | UPDATE_T1 -> (T1) |     |        and write converged results
!    ...........             UPDATE_T2 -> (T2)       |     |
!  Various routines that                             |   CALC_CCSD_ENERGY -> E_CCSD
! compute intermediates and                        CALC_SFACTOR -> Structure Factor
! contract them with T1 and T2
!
! Other experimental routines that compute (T) and manipulate all sorts of quantities and
! contributions to the coupled cluster many-body perturbtaion expansion are not mentioned here.
!  
!
! ================================================
! SOME TECHNICAL DETAILS AND DESIGN PRINCIPLES:
! ================================================
!
!   * THIS MODULE PARALLELIZES OVER A SINGLE K-POINT INDEX ONLY. THEREFORE THE MAXIMUM NUMBER OF
!      ALLOWED MPI-PROCESSES IS LIMITED BY THE NUMBER OF K-POINTS.
!   * ALL QUANTTIES INVOLVING TWO OR MORE KPOINT INDICES ARE DISTRIBUTED USING MPI.
!      QUANTITIES WITH ONE K-POINT INDEX ARE NOT DISTRIBUTED AND REPLICATED ON ALL MPI PROCESSES
!   * THE REQUIRED COULOMB INTEGRALS ARE COMPUTED ON-THE-FLY AND CONTRACTED WITH CORRESPONDING
!      T1 AND T2 AMPLITUDES OR INTERMEDIATE QUANTITIES
!    
! ================================================
! ABOUT THE THEORY
! ================================================
!
!   * THE IMPLEMENTED EXPRESSIONS HAVE BEEN TAKEN FROM 
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
! ===============================================================================================
! THIS MODULE HAS BEEN CONTINUOUSLY MODIFIED AND DEVELOPED BY ANDREAS GRÜNEIS SINCE 2010.
! UNFORTUNATELY DUE TO A LACK OF TIME MANY PARTS ARE NOT YET DOCUMENTED IN A CONSISTENT MANNER
! AND FULLY OPTIMIZED.
! IN FUTURE WE AIM TO COMPLETELY REPLACE THIS MOUDLE WITH A NEW CODE WRITTEN
! IN C++ BASED ON A TENSOR FRAMEWORK AND INTERFACES TO VASP.
!
! TODOs: -Add more documentation 
!        -Remove all unused variables
!        -Improve code profiling FLOP Counters etc.
!        -.....
! ===============================================================================================
! 
!
!*********************************************************************************


MODULE ccsd
   USE prec
   USE fock
   USE chi_base
   USE lattice
   USE wpot
   IMPLICIT NONE

#ifdef scaLAPACK

!*********************************************************************************
!
! Allocation of Coulomb-Vertex Arrays:
!    FTOD_PW and FTOD_OC refer to plane-wave and PAW-one-center part, respectively.
!    FTOD_PW_AB refers to particle-particle part of Coulomb-Vertex.
!    These quantitites are needed to build the corresponding Coulomb integrals by
!    integrating over the auxiliary index:
!      e.g.:  <ab|cd> = SUM_g FTOD_PW_AB(g,n_a,n_c,k_q,k_c,1) * FTOD_PW_AB(g,n_b,n_d,k_q,k_d,2)
! 
!*********************************************************************************
   COMPLEX(q), ALLOCATABLE, PUBLIC :: FTOD_PW(:, :, :), FTOD_PW_AB(:, :, :, :, :, :), FTOD_PW_AI(:, :, :, :, :, :)
   COMPLEX(q), ALLOCATABLE, PUBLIC :: FTOD_PW_IA(:, :, :, :, :, :), FTOD_PW_S_AB(:, :, :, :, :, :), FTOD_PW_S_AI(:, :, :, :, :, :)
   COMPLEX(q), ALLOCATABLE :: FTOD_PW_AI_NOPOT(:, :, :, :, :, :)
   COMPLEX(q), ALLOCATABLE :: FTOD_PW_IJ_NOPOT(:, :, :, :, :, :)
   COMPLEX(q), ALLOCATABLE, PUBLIC :: FTOD_PW_IJ(:, :, :, :, :, :), FTOD_PW_S_IA(:, :, :, :, :, :), FTOD_PW_S_IJ(:, :, :, :, :, :)
   GDEF, ALLOCATABLE, PUBLIC :: FTOD_OC(:, :, :), FTOD_OC_AB(:, :, :, :, :, :), FTOD_OC_AI(:, :, :, :, :, :)
   GDEF, ALLOCATABLE, PUBLIC :: FTOD_OC_IA(:, :, :, :, :, :)
   GDEF, ALLOCATABLE, PUBLIC :: FTOD_OC_IJ(:, :, :, :, :, :)
   GDEF, ALLOCATABLE, PUBLIC ::  OC_IJ_TMP(:, :, :, :), PW_IJ_TMP(:, :, :, :)
   GDEF, ALLOCATABLE, PUBLIC ::  OC_IA_TMP(:, :, :, :), PW_IA_TMP(:, :, :, :), OC_IA_TMP2(:, :, :, :), PW_IA_TMP2(:, :, :, :)
   GDEF, ALLOCATABLE, PUBLIC ::  OC_AI_TMP(:, :, :, :), PW_AI_TMP(:, :, :, :)
   GDEF, ALLOCATABLE, PUBLIC ::  OC_AB_TMP(:, :, :, :), PW_AB_TMP(:, :, :, :)

!*********************************************************************************
!
! Allocation of Intermediate quantities needed for the amplitude equations and elsewhere
! 
!*********************************************************************************
   ! important intermediates for complex orbitals
   GDEF, ALLOCATABLE, PUBLIC :: TE4O(:, :), VVOO(:, :, :, :), OVOV(:, :, :, :), JLIK(:, :, :, :), &
      KLIJ(:, :, :, :), OOVV(:, :, :, :), OVVO(:, :, :, :), VV(:, :), VOVO(:, :, :, :), OVVV(:, :, :, :), &
      OOOV(:, :, :, :), VV2(:, :), VVVV(:, :, :, :)
   ! important intermediates for real orbitals
   REAL(q), ALLOCATABLE, PUBLIC :: TE4O_R(:, :), VVOO_R(:, :, :, :), OVOV_R(:, :, :, :), JLIK_R(:, :, :, :), &
                   KLIJ_R(:, :, :, :), OOVV_R(:, :, :, :), OVVO_R(:, :, :, :), VV_R(:, :), VOVO_R(:, :, :, :), OVVV_R(:, :, :, :), &
                                   OOOV_R(:, :, :, :), VV2_R(:, :), VVVV_R(:, :, :, :)
   ! important intermediates for complex orbitals
   GDEFS, ALLOCATABLE, PUBLIC :: T2(:, :, :, :, :, :, :), T2_N(:, :, :, :, :, :, :), VVVV_S(:, :, :, :), &
      CHI_KLIJ(:, :, :, :, :, :, :), CHI_KLIJ_TMP(:, :, :, :, :, :), OOVV_S(:, :, :, :), VVOO_S(:, :, :, :), &
      T2_TMP(:, :, :, :, :, :), CHI_CKAI(:, :, :, :, :, :, :), VVOO2_S(:, :, :, :), VOVO2_S(:, :, :, :), &
      VOVO_S(:, :, :, :), VOVO3_S(:, :, :, :), T2_KTMP(:, :, :, :, :), T2_KTMP2(:, :, :, :, :), &
      CHI_CKAI_KTMP2(:, :, :, :, :), CHI_CKAI_KTMP(:, :, :, :, :), K_KI(:, :, :), K_AC(:, :, :), K_KC(:, :, :), &
      L_KI(:, :, :), L_AC(:, :, :), OOVV2_S(:, :, :, :), VV_S(:, :), VV2_S(:, :), OVVO_S(:, :, :, :), &
      T1(:, :, :), T1_N(:, :, :), VO_S(:, :), OV_S(:, :), VOOV_S(:, :, :, :), OVVV_S(:, :, :, :), T1_T(:, :, :), &
      OOOV_S(:, :, :, :), VO2_S(:, :), OV2_S(:, :), F_KI(:, :, :), F_AI(:, :, :), F_BA(:, :, :), F_KC(:, :, :), &
      T2_VEL(:, :, :, :, :, :, :), PRET2(:, :, :, :, :, :, :)
   ! important intermediates for basis set correction
   GDEF, ALLOCATABLE :: D2_OOOO(:, :, :,:) ! delta kernel integrals
   GDEF, ALLOCATABLE :: D2PAW_OOOO(:, :, :,:) ! delta kernel integrals
   GDEF, ALLOCATABLE :: D2PAW_VVOO(:, :, :,:,:,:,:) ! delta kernel integrals
   GDEF, ALLOCATABLE :: EMP2_PAIR_CBS(:,:,:,:)
   GDEF, ALLOCATABLE :: EMP2_PAIR_FBS(:,:,:,:)
   GDEF, ALLOCATABLE :: GMP2_PAIR(:,:,:,:) ! mp2 correlation hole depth scaling factors
   GDEF, ALLOCATABLE :: GCC_PAIR(:,:,:,:)  ! cc correlation hole depth scaling factors
   ! important intermediates for real orbitals
   REAL(qs), ALLOCATABLE, PUBLIC :: T2_R(:, :, :, :, :, :, :), T2_N_R(:, :, :, :, :, :, :), VVVV_S_R(:, :, :, :), &
                    CHI_KLIJ_R(:, :, :, :, :, :, :), CHI_KLIJ_TMP_R(:, :, :, :, :, :), OOVV_S_R(:, :, :, :), VVOO_S_R(:, :, :, :), &
                        T2_TMP_R(:, :, :, :, :, :), CHI_CKAI_R(:, :, :, :, :, :, :), VVOO2_S_R(:, :, :, :), VOVO2_S_R(:, :, :, :), &
                                 VOVO_S_R(:, :, :, :), VOVO3_S_R(:, :, :, :), T2_KTMP_R(:, :, :, :, :), T2_KTMP2_R(:, :, :, :, :), &
               CHI_CKAI_KTMP2_R(:, :, :, :, :), CHI_CKAI_KTMP_R(:, :, :, :, :), K_KI_R(:, :, :), K_AC_R(:, :, :), K_KC_R(:, :, :), &
                       L_KI_R(:, :, :), L_AC_R(:, :, :), OOVV2_S_R(:, :, :, :), VV_S_R(:, :), VV2_S_R(:, :), OVVO_S_R(:, :, :, :), &
          T1_R(:, :, :), T1_N_R(:, :, :), VO_S_R(:, :), OV_S_R(:, :), VOOV_S_R(:, :, :, :), OVVV_S_R(:, :, :, :), T1_T_R(:, :, :), &
           OOOV_S_R(:, :, :, :), VO2_S_R(:, :), OV2_S_R(:, :), F_KI_R(:, :, :), F_AI_R(:, :, :), F_BA_R(:, :, :), F_KC_R(:, :, :), &
                                    T2_VEL_R(:, :, :, :, :, :, :), PRET2_R(:, :, :, :, :, :, :)

   GDEFS :: STEP, MU

!*********************************************************************************
!
! NGVECTOR and NHVECTOR refer to the number of auxilliary basis functions used
! for the expansion of the co-densities and Coulomb-Vertex
! 
!*********************************************************************************
   INTEGER, PUBLIC, SAVE :: NGVECTOR, NHVECTOR
!*********************************************************************************
!
! VBMAX and NUNOCC refer to the number of occupied and unoccupied states per K-point.
! H is a VASP object needed to handle PAW terms
! 
!*********************************************************************************
   INTEGER, PUBLIC, SAVE :: VBMAX, NUNOCC, START_BAND
   TYPE(one_center_handle), POINTER, PUBLIC, SAVE :: H

!*********************************************************************************
!
! Some work arrays needed for the parallelization over K-points to different processors
! 
!*********************************************************************************
   integer, ALLOCATABLE, PUBLIC, SAVE :: KPTS_MKPTS(:)
   integer, ALLOCATABLE, PUBLIC, SAVE :: PROCS_KPTS(:)
   integer, ALLOCATABLE, PUBLIC, SAVE :: MKPTS_KPTS(:)
   integer, PUBLIC, SAVE :: MY_NKPTS

!*********************************************************************************
!
! Some steps in the parallelization employ Scalapack routines. We need some variable for BLACS, etc. . :
!  BLACS related variables and function
!  PROCS.. number of processors, ME... processor number, NPROW... number of rows in process grid
!  NPCOL... number of columns in process grid, myrow,mycol... my coordinates in process grid
!  mb... blocking size of rows for block cyclic distribution
!  nb... blocking size of columns for block cyclic distribution
!
!*********************************************************************************

   INTEGER, EXTERNAL :: NUMROC, BLACS_PNUM
   INTEGER, PUBLIC, SAVE  :: CONTXT_COLS, CONTXT_GRID
   INTEGER, PUBLIC, SAVE :: NPROW_GRID, ncc
   INTEGER, PUBLIC, SAVE :: PROCS, ME
   INTEGER, PRIVATE, SAVE :: NPROW, NPCOL, MYROW, MYCOL, NBLOCKAB
   !Scalapack array descriptor for FTOD_PW__ before redistribution
   integer, dimension(9)   :: desc_FTOD_PW_br, desc_FTOD_OC_br
   !Scalapack array descriptor for FTOD_PW_ after redistribution
   integer, dimension(9)   :: desc_FTOD_PW, desc_FTOD_OC

!*********************************************************************************
!
! Some LOGICALs to control different types of CC approximations
! 
!*********************************************************************************

   LOGICAL :: SINGLES, RING, LCCD, CANONICAL, CCMP2, PRECONDITION, LORBREAL, LPSPPL
   LOGICAL :: LMETAL
   LOGICAL :: DAMPED, EH_SCREENING, ASYM_RING, LDISTING, LBRUECKNER, CONVERGED
   COMPLEX(q), ALLOCATABLE, PUBLIC :: E_QP(:, :, :), E_QP_X(:, :, :)
   INTEGER, ALLOCATABLE :: NI_N(:), NA_N(:)

!*********************************************************************************
!
! Central result variables and control variables for various solvers
!
!*********************************************************************************

   COMPLEX(q) :: E_CCSD, E_CCSD_X, E_CCSD_OLD, CDIAG, RES, DELTA_E_CCSD, E_CCSD_FS
   COMPLEX(q) :: E_MP2
   REAL(q) :: AX, TP, HFSCR
   COMPLEX(qs) :: mix
   REAL(q) :: CCMIX
   CHARACTER(LEN=10) ::  DIR_APP
   LOGICAL :: LWRITET
   INTEGER :: NELM
   INTEGER :: MAX_IT, NTHREADS

!*********************************************************************************
!
! This module also contains a prototype implementation for perturbative triples (T)
! that needs various intermediate arrays.
!
!*********************************************************************************

   COMPLEX(qs), ALLOCATABLE :: W_ijk_abc(:, :, :, :, :, :, :)
   COMPLEX(qs), ALLOCATABLE :: Wint_ijk_abc(:, :, :, :, :, :, :)
   COMPLEX(qs), ALLOCATABLE :: Wtmp_ijk_abc(:, :, :, :, :, :, :)
   INTEGER, ALLOCATABLE :: KINDEX(:, :, :)
   INTEGER :: EMSCALCWINT
   COMPLEX(q) :: ETRIPLES

!*********************************************************************************
!
! For shifted (non-gamma-centered k-meshes) we need a set of the following
! work arrays to handle all k-points in the employed k-meshes
! KI refers to k-point of one-electron state
! KQ refers to momentum transfer of two one-electron states
! N.B.: {KQ} will always be gamma-centered
!       {KI} is defined by the KPOINTS file in VASP
!
!*********************************************************************************

   LOGICAL :: SHIFTED_KPOINTS
   REAL(q) :: WTKPT
   INTEGER :: REALNKPTS
   INTEGER, ALLOCATABLE :: RKQofKQ(:), RKIofKI(:)
   INTEGER, ALLOCATABLE :: KQofRKQ(:), KIofRKI(:)

!*********************************************************************************
!
! For the calculation of the electronic transition structure factor we need
! the following work arrays. Furthermore a number of experimental structure factor
! related quantities are also included
!
!*********************************************************************************

   LOGICAL :: LSFACTOR
   REAL(q), ALLOCATABLE :: POTFAK_FULL(:, :), GVECLEN(:, :), GVECX(:, :), GVECY(:, :), GVECZ(:, :)
   COMPLEX(q), ALLOCATABLE :: SDFACTOR(:, :), SXFACTOR(:, :)
   REAL(q), ALLOCATABLE :: REG_SFACTOR(:,:,:) !, REG_SXFACTOR(:,:,:)
   REAL(q), ALLOCATABLE :: REG_X(:), REG_Y(:), REG_Z(:)
   COMPLEX(q), ALLOCATABLE :: STFACTOR(:, :), SSFACTOR(:, :)
   COMPLEX(q), ALLOCATABLE :: F12TFACTOR(:, :), F12SFACTOR(:, :)
   COMPLEX(q), ALLOCATABLE :: QDMATRIX(:, :), QXMATRIX(:, :)
   COMPLEX(q), ALLOCATABLE :: QSMATRIX(:, :), QTMATRIX(:, :)
   COMPLEX(q), ALLOCATABLE :: QSMATRIX_INV(:, :), QTMATRIX_INV(:, :)

!*********************************************************************************
!
! These are arrays for a number of unfinished experimental features that are
! and not yet tested.
!
!*********************************************************************************
   REAL(8), ALLOCATABLE :: MF_R(:, :), MMU_R(:, :), T2TMP_R(:)
   COMPLEX*16, ALLOCATABLE :: MF(:, :), MMU(:, :), T2TMP(:)
   COMPLEX*16, ALLOCATABLE :: MFS(:, :), MMUS(:, :), T1TMP(:)
   INTEGER :: DIISMAXSD, NVECL, NVECLS, nit_diis
!!!!! DIIS for Brueckner orbitals not used anymore
   COMPLEX*16, ALLOCATABLE :: WF(:, :), WMUS(:, :), WTMP(:)
   INTEGER :: DIISMAXSDW, NVECLW, nit_diisw
!!!!!  FOR BRUECKNER
   INTEGER :: MAX_IT_BRUECKNER
! Some other variables used for testing performance in an experimental version
   INTEGER :: TETS, TET2S, TETSORT
!
   REAL(q) :: NFLOAT4O, tetfpo, tetcom, tetr2d, tetvv, tetcc, tetvc, tetcontv
   COMPLEX(q) :: E_N, E_N_X, DELTA
!
   INTEGER :: CCNINT


CONTAINS

!***********************************************************************
!
! VASP distinguishes between gamma-only (gamma-point only and real orbitals)
! and non-gamma-only (arbitrary k-meshes and complex orbitals)
! implementations using pre-compiler statements. This CCSD module only works in the
! non-gamma-only framework.
!
!***********************************************************************
#ifdef gammareal
      SUBROUTINE CALCULATE_CCSD(SYMM,P,WDES,W,LATT_CUR,T_INFO,INFO,IO,KPOINTS,WGW,ENCUTGW, ENCUTGWSOFT, LMAXMP2,AMIX,BMIX,NONLR_S,NONL_S,GRID,LATT_INI)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE base
      USE mkpoints
      USE mpimy
      USE ini
      USE choleski
      USE hamil_struct_def
      USE fileio
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      TYPE(wavespin) W
      TYPE(type_info) T_INFO
      TYPE(info_struct) INFO
      TYPE(potcar) P(T_INFO%NTYP)
      TYPE(latt) LATT_CUR, LATT_INI
      TYPE(grid_3d) GRID
      RGRID SV(DIMREAL(GRID%MPLWV), WDES%NCDIJ) 
      TYPE(nonlr_struct) NONLR_S
      TYPE(nonl_struct) NONL_S
      TYPE(symmetry) :: SYMM
      INTEGER LMAXMP2, NBANDSGW
      REAL(q) :: ENCUTGW, ENCUTGWSOFT, change, denom, AMIX, BMIX
      TYPE(in_struct) IO
      TYPE(kpoints_struct) KPOINTS

   END SUBROUTINE CALCULATE_CCSD
#else

      SUBROUTINE CALCULATE_CCSD(SYMM,P,WDES,W,LATT_CUR,T_INFO,INFO,IO,KPOINTS,WGW,ENCUTGW, ENCUTGWSOFT, LMAXMP2,AMIX,BMIX,NONLR_S,NONL_S,GRID,LATT_INI)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE base
      USE mkpoints
      USE mpimy
      USE ini
      USE choleski
      USE hamil_struct_def
      USE fileio
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      TYPE(wavespin) W
      TYPE(type_info) T_INFO
      TYPE(info_struct) INFO
      TYPE(potcar) P(T_INFO%NTYP)
      TYPE(latt) LATT_CUR, LATT_INI
      TYPE(grid_3d) GRID       ! grid for wavefunctions
      TYPE(nonlr_struct) NONLR_S
      TYPE(nonl_struct) NONL_S
      TYPE(symmetry) :: SYMM
      INTEGER LMAXMP2, NBANDSGW
      REAL(q) :: ENCUTGW, ENCUTGWSOFT, change, denom, AMIX, BMIX
      TYPE(in_struct) IO
      TYPE(kpoints_struct) KPOINTS
      ! local
      INTEGER :: NI, NJ, NR, NS, KI, KJ, KR, KS, KQ, KQ_, direction, KR_, NA, MKA, MKJ, &
                 GKQ, MKB, KC, ITERATION, ITERATION_BRU, START_IT, ITBRUECKNER, NBRUECKNER
      integer :: TWOE4ORBITAL_COLS, TWOE4ORBITAL_ROWS, niteration, nwfup
      integer :: time_array1(8), time_array2(8), ems, KN, NN, omp_get_num_threads


!***********************************************************************
!
! Initialize some system dependent variables and parameters used for the solver
!
!***********************************************************************

      ncc = 2

      nthreads = 1 ! OMP implementation not used anymore

      LORBREAL = .FALSE.
      IF (WDES%LORBITALREAL) LORBREAL = .TRUE.
      AMIX = 0.25_q

      VBMAX = 0
      DO KI = 1, W%WDES%NKPTS
         VBMAX = MAX(LAST_FILLED_XI_NOMOD(W, KI, 1), VBMAX)
      END DO

      LMETAL = .FALSE.
      IF (VBMAX > (INFO%NELECT/2)) LMETAL = .TRUE.
      IF (IO%IU0 > 0) THEN
         WRITE (*, *) 'VBMAX=', VBMAX, 'NELECT=', INFO%NELECT, "LMETAL=", LMETAL
      END IF

      DIISMAXSD = 4
      DIISMAXSDW = 3
      nit_diis = 1

      MAX_IT = 10
      CCMP2 = .FALSE.
      CANONICAL = .TRUE.
      SINGLES = .TRUE.
      LCCD = .FALSE.  !double-checked with older version using diamond till 3x3x3 kpoints
      RING = .TRUE.  !double-checked with ACFDT version using diamond till 3x3x3 kpoints
      LDISTING = .FALSE.
      LBRUECKNER = .FALSE.
      MAX_IT_BRUECKNER = 1
      IF (LBRUECKNER) SINGLES = .TRUE.
      IF (LBRUECKNER) CANONICAL = .FALSE.
      IF (LBRUECKNER) MAX_IT_BRUECKNER = 20
      E_CCSD = zero
      E_CCSD_X = zero
      CONVERGED = .FALSE.

!***********************************************************************
!
! Read/Overwrite (additional) parameters that can be controlled using the INCAR file of VASP
!
!***********************************************************************

      CALL CCINCAR_READER(IO%IU5, IO%IU6, IO%IU0)
      mix = (1.0_q, 0.0_q)*CCMIX
      IF (IO%IU0 > 0) THEN
         WRITE (*, *) 'Linear mixer is set to ', mix
         WRITE (*, *) 'LDISTING=', LDISTING
      END IF

      ASYM_RING = .FALSE.
      DAMPED = .FALSE.
      PRECONDITION = .FALSE.
      CDIAG = (0.0_q, 0.0_q)
      EH_SCREENING = .FALSE.  !  experimental part that employs screening for some diagrams
      AX = 0.5
      TP = 8._q*ATAN(1._q)
      HFSCR = 0.0
      MU = 1.0
      STEP = 0.5
      DELTA = (0.0_qs, 0.0_qs)
      IF (RING) SINGLES = .FALSE.
      IF (LCCD) SINGLES = .FALSE.

!***********************************************************************
!
! Check which type of k-mesh is used
!
!***********************************************************************

      CALL CHECK_FULL_KPOINTS
      CALL CHECK_SHIFTED_KPOINTS(WDES, W, KPOINTS)
      IF ((IO%IU0 > 0) .and. (SHIFTED_KPOINTS)) write (IO%IU0, *) 'You are using a shifted k-mesh.'
      !Initialize the 1D process grid
      CALL INIT_BLACS_COLS()

      CALL INIT_MY_KPOINTS(WDES, WGW)

!***********************************************************************
!
! Compute Coulomb-Vertex
!
!***********************************************************************

      CALL CALC_FTOD(WDES, WGW, W, P, T_INFO, LATT_CUR, LMAXMP2, ENCUTGW, ENCUTGWSOFT, IO, FSG_STORE(1), 0)

!***********************************************************************
!
! Write out system-specific settings: number of occupied and unoccupied state, k-points 
!
!***********************************************************************

      IF (IO%IU0 > 0) write (*, *) ''
      IF (IO%IU0 > 0) WRITE (*, *) 'VBMAX=', VBMAX
      IF (IO%IU0 > 0) WRITE (*, *) 'NUNOCC=', NUNOCC
      IF (IO%IU0 > 0) WRITE (*, *) 'WDES%NKPTS=', WDES%NKPTS
      IF (IO%IU0 > 0) WRITE (*, *) 'KPOINTS_ORIG%NKPTS=', KPOINTS_ORIG%NKPTS

!***********************************************************************
!
! NBLOCKAB controls the maximum block size for the on-the-fly calculation of the
! <ab|cd> Integrals used to contract with the amplitudes t^cd_ij .
! Increase for higher efficency and larger memory footprint.
!
!***********************************************************************

      NBLOCKAB = MIN(80, NUNOCC) !WDES%NB_TOT-VBMAX)

!***********************************************************************
!
! Initialise the amplitudes.
!
!***********************************************************************
      IF (IO%IU0 > 0) WRITE (*, *) 'setting up T2 amplitudes'
      CALL SETUP_T2AMPLITUDES(WDES)
      IF (.not. CANONICAL) CALL SETUP_FOCKMATRIX(W, WGW, WDES, IO)

      IF (IO%IU0 >= 0) write (*, *) ''
      IF (IO%IU0 >= 0) WRITE (IO%IU0, *) 'Calculating T2AMPLITUDES'
      IF (IO%IU0 >= 0) write (*, *) ''
      IF (.not. WDES%LORBITALREAL) THEN
         T2 = (0.0_qs, 0.0_qs)
         T2_N = (0.0_qs, 0.0_qs)
         T1 = (0.0_qs, 0.0_qs)
         T1_N = (0.0_qs, 0.0_qs)
      ELSE
         T2_R = (0.0_qs)
         T2_N_R = (0.0_qs)
         T1_R = (0.0_qs)
         T1_N_R = (0.0_qs)
      END IF

      START_IT = 1

!***********************************************************************
!
! Check if T1 and T2 amplitudes have already been written to disk from a
! previous calculation. Read them in if they can be found.
!
!***********************************************************************

      CALL IN_T2(IO, WDES, W, START_IT)
      CALL IN_T1(IO, WDES, W, START_IT)

      IF (LPSPPL) CALL IN_EMP2_PAIR_CBS(IO, WDES, W, START_IT)

!***********************************************************************
!
! Main loop for solving the CCSD amplitude equations. Break before MAX_IT
! if desired convergence has been reached. Otherwise break after MAX_IT
! iterations.
!
!***********************************************************************

      DO ITERATION = START_IT, MAX_IT

         CALL START_TIMING("G")

         CALL CALC_T1_T2_N(WDES, WGW, W, ITERATION, KPOINTS)
         IF (IO%IU0 > 0) write (*, *) 'updating T2 amplitudes'
         CALL UPDATE_T2(WDES, WGW, W, IO, ITERATION)
         IF (IO%IU0 > 0) write (*, *) 'updating T1 amplitudes'
         IF (SINGLES) CALL UPDATE_T1(WDES, WGW, W, IO, ITERATION)

         IF (IO%IU0 > 0) write (*, *) 'calculating CC energy'
         E_CCSD_OLD = E_CCSD
         CALL CALC_CCSD_ENERGY(WDES, WGW, W, IO, ITERATION)
         
         IF ((ITERATION==1) .and. (LPSPPL)) THEN
            CALL CALC_D2PAW(LATT_CUR, W, WDES, T_INFO, P, IO, ENCUTGW, ENCUTGWSOFT)
            CALL CALC_GMP2_PAIR(LATT_CUR, WGW, WDES, IO, ENCUTGW, ENCUTGWSOFT)
!            STOP
         ENDIF

         IF (ABS(E_CCSD_OLD - E_CCSD) < INFO%EDIFF) THEN
            CONVERGED = .TRUE.
       IF (IO%IU6 >= 0) WRITE (IO%IU6, *) 'CCSD energy appears to be converged: Delta E=', ABS(E_CCSD_OLD - E_CCSD), '<', INFO%EDIFF
            IF (IO%IU6 >= 0) WRITE (*, *) 'CCSD energy appears to be converged: Delta E=', ABS(E_CCSD_OLD - E_CCSD), '<', INFO%EDIFF
         END IF

!***********************************************************************
!
! Write T1 and T2 amplitudes to disk after every iteration. 
!
!***********************************************************************
         IF (LWRITET) THEN
            IF (.not. CCMP2) CALL OUT_T2(IO, WDES, W)
            CALL OUT_T1(IO, WDES, W)
         END IF

         CALL STOP_TIMING("G", IO%IU6, "LOOP CCSD")

         IF (((ITERATION==MAX_IT) .and. (.not. CONVERGED)) .and. (IO%IU6>=0)) WRITE(IO%IU6,*)'final E_CCSD (WARNING not converged!)=',E_CCSD
         IF (CONVERGED) THEN
!***********************************************************************
!
! Print the converged correlation energy.
!
!***********************************************************************
            IF (IO%IU6 >= 0) WRITE (IO%IU6, *) 'final E_CCSD=', E_CCSD
            EXIT
         END IF

      end do

!***********************************************************************
!
! Compute the electronic transition structure factor and energy (again)
! and write everything to disk and output files.
!
!***********************************************************************

      IF ((LSFACTOR)) CALL CALC_SFACTOR(LATT_CUR, WGW, WDES, IO, ENCUTGW, ENCUTGWSOFT)

      CALL CALC_CCSD_ENERGY(WDES, WGW, W, IO, ITERATION)

      E_N = E_CCSD
      E_N_X = E_CCSD_X

      IF (IO%IU0 >= 0) WRITE (*, *) 'CCSD correlation energy:', E_CCSD
      IF (IO%IU0 >= 0) WRITE (IO%IU6, *) 'CCSD correlation energy:', E_CCSD
      IF (IO%IU0 >= 0) WRITE (*, *) 'CCSD_X correlation energy:', E_CCSD_X
      IF (IO%IU0 >= 0) WRITE (IO%IU6, *) 'CCSD_X correlation energy:', E_CCSD_X

      IF ((LPSPPL)) THEN
!         CALL CALC_GMP2_PAIR(LATT_CUR, WGW, WDES, IO, ENCUTGW, ENCUTGWSOFT)
         CALL CALC_PSPPL(LATT_CUR, WGW, WDES, IO, ENCUTGW, ENCUTGWSOFT)
      ENDIF

!***********************************************************************
!
! Some experimental features for the particle-particle-ladder term
!
!***********************************************************************

      IF (.not. CCMP2) THEN
         IF (IO%IU0 > 0) write (*, *) 'calculating pp ladder-only T2 amplitudes'
         LCCD = .TRUE.
         CALL CALC_T1_T2_N_PPL(WDES, WGW, W, ITERATION, KPOINTS)
         CALL UPDATE_T2_LADDER(WDES, WGW, W, IO, ITERATION)
         CALL CALC_CCSD_ENERGY_PPL(WDES, WGW, W, IO, ITERATION)
      END IF

      IF (LWRITET) THEN
         IF (.not. CCMP2) CALL OUT_T2(IO, WDES, W)
         CALL OUT_T1(IO, WDES, W)
         IF ((.not. CCMP2)) CALL OUT_FTOD_PW(IO, WDES, W)
      END IF

!***********************************************************************
!
! Exit BLACS and bye, bye!
!
!***********************************************************************

      call BLACS_GRIDEXIT(contxt_cols)
      call BLACS_GRIDEXIT(contxt_grid)

   END SUBROUTINE CALCULATE_CCSD

!***********************************************************************
!
!  This routine calculates the fourier-transformed overlap integrals(FTOD) <i|-G|a>
!  and <j|G|b>*4*pi*e^2/(G+q)^2 aka Coulomb-Vertex and calls the redistribution routine.
!  This is done for the plane-wave part, as well as for the one-center terms.
!  Note, that in the gamma-only version it is enough to save <i|-G|a>*SQRT(4*pi*e^2/(G+q)^2)
!  for the plane-wave part, because (<i|-G|a>)*=<i|G|a>.
!  However, this approximation cannot be made for the one-center terms in the gamma-only
!  version, because of practical reasons.
!  The <i|G|a> quantities are most important for the construction of the two electron
!  four orbital integrals <ij|ab>.  <ij|ab>=4pi e^2 \sum_G <i|-G|a> <j|G|b> /(G+q)^2
!
!***********************************************************************

   SUBROUTINE CALC_FTOD(WDES, WGW, W, P, T_INFO, LATT_CUR, LMAXMP2, ENCUTGW, ENCUTGWSOFT, IO, FSG, IT)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      TYPE(wavespin) W
      TYPE(type_info) T_INFO
      TYPE(potcar) P(T_INFO%NTYP)
      TYPE(latt) LATT_CUR
      INTEGER LMAXMP2
      REAL(q) :: ENCUTGW, ENCUTGWSOFT
      TYPE(in_struct) IO
      REAL(q) :: FSG ! singularity correction
      INTEGER IT
      ! local variables
      TYPE(wavespin) WHF
      TYPE(wavedes1) WGWQ
      TYPE(wavedes1) WDESKI, WDESKA, WDESKB
      TYPE(wavefun1), ALLOCATABLE :: WI(:), WA(:), WB(:)
      INTEGER KQ, KI, KA, KB, KI_IN_FULL_ORIG, KQ_, MKQ
      INTEGER NSTRIP, NSTRIPA, ISP, NFFT
      INTEGER NBI, NBA, NBAA, i, nbi_start, nbi_end, rnba
      INTEGER NP, MNP ! number of plane waves for the overlap density i*(r)a(r) GCHGIA(r)
      REAL(q) :: POTFAK(GRIDHF%MPLWV), MSCALE
      REAL(q) :: POTFAK_S(GRIDHF%MPLWV)
      COMPLEX(q) CPHASE(GRIDHF%MPLWV)
      COMPLEX(q) CPHASE2(GRIDHF%MPLWV)
      LOGICAL LPHASE
      LOGICAL LPHASE2
      COMPLEX(q) :: GWORK(MAX(GRIDHF%MPLWV, WGW%GRID%MPLWV)) !work array for calculating GCHGIA
      COMPLEX(q), ALLOCATABLE :: GCHGIA(:, :, :)  ! charge
      COMPLEX(q), ALLOCATABLE :: GCHGIA_S(:, :, :)  ! charge
      GDEF, ALLOCATABLE :: CRHOIA(:, :)    ! one-center charge
      GDEF, ALLOCATABLE :: CRHOIB(:, :)    ! one-center charge
      GDEF, ALLOCATABLE :: CRHOLM(:)      ! augmentation occupancy matrix
      COMPLEX(q), ALLOCATABLE :: tmp_FTOD_PW(:, :, :, :)
      COMPLEX(q), ALLOCATABLE :: tmp_FTOD_PW_S(:, :, :, :)
      GDEF, ALLOCATABLE :: tmp_FTOD_OC(:, :, :, :)
      Real(qs) :: mem_req

      NGVECTOR = MAXVAL(WGW%NGVECTOR(:))

!***********************************************************************
!
! Structure factor related quantities including experimental stuff
!
!***********************************************************************

      IF (LSFACTOR) THEN
         IF (.not. ALLOCATED(POTFAK_FULL)) THEN
            ALLOCATE (POTFAK_FULL(NGVECTOR, REALNKPTS))
            ALLOCATE (SDFACTOR(NGVECTOR, REALNKPTS))
            ALLOCATE (SXFACTOR(NGVECTOR, REALNKPTS))
            ALLOCATE (STFACTOR(NGVECTOR, REALNKPTS))
            ALLOCATE (SSFACTOR(NGVECTOR, REALNKPTS))
            ALLOCATE (F12TFACTOR(NGVECTOR, REALNKPTS))
            ALLOCATE (F12SFACTOR(NGVECTOR, REALNKPTS))
            ALLOCATE (QDMATRIX(NGVECTOR*REALNKPTS, NGVECTOR*REALNKPTS))
            ALLOCATE (QXMATRIX(NGVECTOR*REALNKPTS, NGVECTOR*REALNKPTS))
            ALLOCATE (QSMATRIX(NGVECTOR*REALNKPTS, NGVECTOR*REALNKPTS))
            ALLOCATE (QTMATRIX(NGVECTOR*REALNKPTS, NGVECTOR*REALNKPTS))
            ALLOCATE (QSMATRIX_INV(NGVECTOR*REALNKPTS, NGVECTOR*REALNKPTS))
            ALLOCATE (QTMATRIX_INV(NGVECTOR*REALNKPTS, NGVECTOR*REALNKPTS))
            POTFAK_FULL = 0.0_q
            SDFACTOR = zero
            SXFACTOR = zero
            STFACTOR = zero
            SSFACTOR = zero
            F12TFACTOR = zero
            F12SFACTOR = zero
         END IF
         IF (.not. ALLOCATED(GVECLEN)) THEN
            ALLOCATE (GVECLEN(NGVECTOR, REALNKPTS))
            ALLOCATE (GVECX(NGVECTOR, REALNKPTS))
            ALLOCATE (GVECY(NGVECTOR, REALNKPTS))
            ALLOCATE (GVECZ(NGVECTOR, REALNKPTS))
            GVECLEN = 0.0_q
            GVECX = 0.0_q
            GVECY = 0.0_q
            GVECZ = 0.0_q
         END IF
      END IF

      CALL CHECK_FULL_KPOINTS

      IF (LMAXMP2 >= 0) THEN
         CALL SET_UP_ONE_CENTER_H(WDES, P, T_INFO, LMAXMP2, H)
      END IF
      WHF = W
      WHF%WDES => WDES_FOCK
      NSTRIP = 30
      CALL SETWDES(WHF%WDES, WDESKI, 0)
      CALL SETWDES(WHF%WDES, WDESKA, 0)
      CALL SETWDES(WHF%WDES, WDESKB, 0)

      NUNOCC = WDES%NB_TOT - VBMAX
      IF (LMETAL) NUNOCC = WDES%NB_TOT

!***********************************************************************
!
! We will compute slices of co-densities to reduce memory footprint at this stage.
!
!***********************************************************************

      ALLOCATE (WI(NSTRIP), WA(NSTRIP), WB(NSTRIP))
      DO NBI = 1, NSTRIP
         CALL NEWWAV(WI(NBI), WDESKI, .TRUE.)
      END DO
      DO NBA = 1, NSTRIP
         CALL NEWWAV(WA(NBA), WDESKA, .TRUE.)
      END DO
      DO NBA = 1, NSTRIP
         CALL NEWWAV(WB(NBA), WDESKB, .TRUE.)
      END DO

      NGVECTOR = MAXVAL(WGW%NGVECTOR(:))
      NHVECTOR = 0
      IF (ASSOCIATED(H)) THEN
         NHVECTOR = H%TOTAL_ENTRIES
      END IF
      ! Make sure that NHVECTOR is suited for the BLACS process grid
      IF (MOD(NHVECTOR, NPROW_GRID) /= 0) THEN
         NHVECTOR = NHVECTOR + (NPROW_GRID - MOD(NHVECTOR, NPROW_GRID))
      END IF

      !Initialize the 2D process grid
      CALL INIT_BLACS_GRID(WDES)

      IF (IT == 0) THEN
         IF (ALLOCATED(FTOD_PW)) THEN
            DEALLOCATE (FTOD_PW)
         END IF
         CALL SETUP_FTOD(WDES)
      END IF

      !Setup the descriptors for the distributed Coulomb Vertex
      ALLOCATE (tmp_FTOD_PW(NGVECTOR, WDES%NBANDS, NSTRIP, ncc))
      IF (EH_SCREENING) ALLOCATE (tmp_FTOD_PW_S(NGVECTOR, WDES%NBANDS, NSTRIP, 1))
      IF (ASSOCIATED(H)) THEN
         ALLOCATE (tmp_FTOD_OC(NHVECTOR, WDES%NBANDS, NSTRIP, 2))
      END IF

      IF (IO%IU0 > 0) THEN
         WRITE (IO%IU0, *)
         WRITE (IO%IU0, *) 'Calculating fourier transformed overlap densities:'
      END IF

      call BLACS_GRIDINFO(CONTXT_COLS, NPROW, NPCOL, MYROW, MYCOL)

      IF (ASSOCIATED(H)) THEN
         ALLOCATE (CRHOIA(NHVECTOR, NSTRIP), CRHOIB(NHVECTOR, NSTRIP))
         FTOD_OC = zero
      END IF

      spin: DO ISP = 1, WDES%ISPIN
      kqloop: DO KQ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.00_q)) CYCLE
         IF (IO%IU0 >= 0) WRITE (IO%IU0, *)
         IF (IO%IU0 >= 0) THEN
            IF (WDES%ISPIN == 1) THEN
               WRITE (IO%IU0, '("NQ=",I4,3F10.4,", ")') KQ, WDES%VKPT(:, KQ)
            ELSE
               WRITE (IO%IU0, '("NQ=",I4,3F10.4,", ",A1,A1,", ")') KQ, WDES%VKPT(:, KQ), ISP
            END IF
         END IF
         CALL SETWDES(WGW, WGWQ, KQ)

         NP = WGWQ%NGVECTOR
         IF (NP > NGVECTOR) THEN
            WRITE (*, *) 'Internal error in "Calc_2orbital_response": NP larger than NGVECTOR'
            EXIT
         END IF
         ALLOCATE (GCHGIA(NP, NSTRIP, 2), CRHOLM(AUG_DES%NPRO*WDES%NRSPINORS))
         IF (EH_SCREENING) ALLOCATE (GCHGIA_S(NP, NSTRIP, 2))
         kiloop: DO KI = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
            tmp_FTOD_PW = (0._q, 0._q)
            IF (ASSOCIATED(H)) tmp_FTOD_OC = zero

            CALL GWPROGRESS(IO%IU0, KI, WDES%NKPTS, KQ, WDES%NKPTS)
            KI_IN_FULL_ORIG = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:, KI), KPOINTS_FULL_ORIG)
            ! collect NSTRIP nbi bands at k_i
            CALL SETWDES(WHF%WDES, WDESKI, KI)

            DO NBI_start = 1, PROCS*WDES%NBANDS, NSTRIP
               call BLACS_GRIDINFO(CONTXT_COLS, NPROW, NPCOL, MYROW, MYCOL)
               NBI_end = NBI_start + MIN(PROCS*WDES%NBANDS - NBI_start + 1, NSTRIP) - 1

               CALL W1_GATHER_GLB(WHF, NBI_start, NBI_end, ISP, WI)
               ! k_b = k_i - k_q - G
               KB = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KQ) + WDES%VKPT(:, KI), KPOINTS_FULL)
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) Write (*, *) 'error in calc_ftod with shifted k-mesh'
               ! k_a = k_i + k_q - G
               KA = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - WDES%VKPT(:, KQ), KPOINTS_FULL)
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA) == 0)) Write (*, *) 'error in calc_ftod with shifted k-mesh'

               CALL SETWDES(WHF%WDES, WDESKA, KA)

               ! CPHASE(r) = e^iGr, where G = k_i - k_q - k_b
               CALL SETPHASE(WDES%VKPT(:, KI) - WDES%VKPT(:, KQ) - WDES%VKPT(:, KA), GRIDHF, CPHASE, LPHASE)

               POTFAK(:) = 0.0_q

               CALL SET_GFAC_WITHOUT_WEIGHT(GRIDHF, LATT_CUR, KI, KA, FSG, POTFAK)
               ! 1/(G+q)**2
                  IF (LSFACTOR) CALL SET_GVEC(WGWQ, LATT_CUR, FSG, GVECLEN(1,RKQofKQ(KQ)), GVECX(1,RKQofKQ(KQ)), GVECY(1,RKQofKQ(KQ)), GVECZ(1,RKQofKQ(KQ)))

               IF (ENCUTGW /= ENCUTGWSOFT .AND. ENCUTGW > 0 &
                   .AND. ENCUTGWSOFT > 0) THEN
                  CALL SET_GFAC_WAVEFUN(WGWQ, LATT_CUR, FSG, POTFAK, ENCUTGW, ENCUTGWSOFT)
                  IF (EH_SCREENING) THEN
                     HFSCREEN = HFSCR
                     CALL SET_GFAC_WAVEFUN(WGWQ, LATT_CUR, FSG, POTFAK_S, ENCUTGW, ENCUTGWSOFT)
                     HFSCREEN = 0
                  END IF
               ELSE
                  CALL SET_GFAC_WAVEFUN(WGWQ, LATT_CUR, FSG, POTFAK)
                  IF (EH_SCREENING) THEN
                     HFSCREEN = HFSCR
                     CALL SET_GFAC_WAVEFUN(WGWQ, LATT_CUR, FSG, POTFAK_S)
                     HFSCREEN = 0
                  END IF
               END IF
               IF (LSFACTOR) THEN
                  POTFAK_FULL(1:NGVECTOR, RKQofKQ(KQ)) = POTFAK(1:NGVECTOR)
                  POTFAK = 1.0_q
               END IF


!
! This is for an experimental feature that employs screeened interaction for the PH-ladder diagrams
!
               IF (EH_SCREENING) THEN
                  IF (ME == 0) WRITE (*, *) 'calculating statically screened e-h interaction from POTFAK'
                  IF ((KQ == 1) .and. (KI == 1) .and. (NBI_start == 1)) write (*, *) 'MSCALE', MSCALE, AX, TP, HFSCR, AEXX
               END IF

               ! loop over all bands
               DO NBA = 1, WDES%NBANDS, NSTRIP
                  NSTRIPA = MIN(WDES%NBANDS + 1 - NBA, NSTRIP)
                  ! FFT{psi_a} to real space
                  DO NBAA = 1, NSTRIPA !copy and fourier transform NSTRIP wave functions
                     CALL W1_COPY(ELEMENT(WHF, WDESKA, NBA + NBAA - 1, ISP), WA(NBAA))
                     CALL FFTWAV_W1(WA(NBAA))
                  END DO
                  ! loop over valence bands only
                  DO NBI = 1, MIN(PROCS*WDES%NBANDS - NBI_start + 1, NSTRIP)
                     GCHGIA = 0
                     IF (EH_SCREENING) GCHGIA_S = 0

                     IF (ASSOCIATED(H)) THEN
                        CRHOIA = zero
                        CRHOIB = zero
                        CRHOLM = zero
                     END IF
                     !loop over all bands in NSTRIP
                     DO NBAA = 1, NSTRIPA
                        CALL LOC2GLOB((NBA + NBAA - 1), ME, WDES%NB_TOT, PROCS, 1, rnba)
                        IF ((CCMP2) .and. ((NBI + NBI_START - 1) > VBMAX) .and. ((RNBA) > VBMAX)) CYCLE
                        ! calculate rho(r)=psi_i(r)* psi_a(r) for,
                        ! one center terms and, on the plane wave grid.
                        IF (ASSOCIATED(H)) THEN
                           CALL FOCK_CHARGE_ONE_CENTER_NOINT(WI(NBI), WA(NBAA), &
                                                             GWORK(1), H, CRHOIA(1, NBAA), CRHOLM, SIZE(CRHOLM))
                        ELSE
                           CALL FOCK_CHARGE_NOINT(WI(NBI), WA(NBAA), GWORK(1), &
                                                  CRHOLM, SIZE(CRHOLM))
                        END IF
                        ! Set phase e^iGr, where G = k_i - k_q - k_b
                        IF (LPHASE) THEN
                           CALL APPLY_PHASE(GRIDHF, CPHASE(1), GWORK(1), &
                                            GWORK(1))
                        END IF

                        ! FFT{rho} to reciprocal space
                        CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), &
                                    GWORK(1), GCHGIA(1, NBAA, 1), WGWQ%GRID, .FALSE.)
                        NFFT = NFFT + 1
                        ! multiply with potential factor
                        IF (EH_SCREENING) THEN
                           GCHGIA_S(:, NBAA, 1) = GCHGIA(:, NBAA, 1)
                           CALL APPLY_GFAC_WAVEFUN(WGWQ, GCHGIA_S(1, NBAA, 1), &
                                                   POTFAK_S(1))
                           CALL APPLY_GFAC_WAVEFUN(WGWQ, GCHGIA(1, NBAA, 1), &
                                                   POTFAK(1))
                        ELSE
                           CALL APPLY_GFAC_WAVEFUN(WGWQ, GCHGIA(1, NBAA, 1), &
                                                   POTFAK(1))
                        END IF

                     END DO !NBAA (loop over all bands in NSTRIP)

                     IF (ASSOCIATED(H)) THEN
                        CALL APPLY_PHASE_ONE_CENTER(WHF%WDES, H, CRHOIA(:, :), &
                                                    WHF%WDES%VKPT(:, KI) - WHF%WDES%VKPT(:, KA))
                     END IF

                     IF (ASSOCIATED(H)) THEN
                        ! use CRHOIA as temporary work array
                        CALL APPLY_ONE_CENTER_H(WHF%WDES, H, CRHOIA(:, :), CRHOIB(:, :), NSTRIPA)
                     END IF

                     DO NBAA = 1, NSTRIPA
                        !copy response functions to RESPF_PW and RESPF_OC
                        CALL LOC2GLOB((NBA + NBAA - 1), ME, WDES%NB_TOT, PROCS, 1, rnba)
                        IF ((CCMP2) .and. ((NBI + NBI_START - 1) > VBMAX) .and. ((RNBA) > VBMAX)) CYCLE
                        tmp_FTOD_PW(1:NP, NBAA + NBA - 1, NBI, 1) = (GCHGIA(1:NP, NBAA, 1))
                        IF (EH_SCREENING) tmp_FTOD_PW_S(1:NP, NBAA + NBA - 1, NBI, 1) = (GCHGIA_S(1:NP, NBAA, 1))
                        IF (ASSOCIATED(H)) THEN
                           tmp_FTOD_OC(1:H%TOTAL_ENTRIES, NBAA + NBA - 1, NBI, 1) = (CRHOIA(1:H%TOTAL_ENTRIES, NBAA))
                        END IF

                     END DO

                  END DO !NBI (loop over nstrip nbi bands only)
               END DO !NBA (loop over all bands)


               CALL SETWDES(WHF%WDES, WDESKB, KB)

               ! k_i - k_a = k_q + G
               CALL PHASER_HF(GRIDHF, LATT_CUR, FAST_AUG_FOCK, WDES%VKPT(:, KB) - WDES%VKPT(:, KI))

               ! CPHASE(r) = e^iGr, where G = k_i - k_q - k_a
               CALL SETPHASE(WDES%VKPT(:, KB) - WDES%VKPT(:, KQ) - WDES%VKPT(:, KI), GRIDHF, CPHASE, LPHASE)

               ! loop over all bands
               DO NBA = 1, WDES%NBANDS, NSTRIP
                  NSTRIPA = MIN(WDES%NBANDS + 1 - NBA, NSTRIP)
                  ! FFT{psi_a} to real space
                  DO NBAA = 1, NSTRIPA !copy and fourier transform NSTRIP wave functions
                     CALL W1_COPY(ELEMENT(WHF, WDESKB, NBA + NBAA - 1, ISP), WB(NBAA))
                     CALL FFTWAV_W1(WB(NBAA))
                  END DO
                  ! loop over valence bands only
                  DO NBI = 1, MIN(PROCS*WDES%NBANDS - NBI_start + 1, NSTRIP)
                     GCHGIA = 0

                     IF (ASSOCIATED(H)) THEN
                        CRHOIB = 0
                        CRHOLM = 0
                     END IF
                     !loop over all bands in NSTRIP
                     DO NBAA = 1, NSTRIPA
                        CALL LOC2GLOB((NBA + NBAA - 1), ME, WDES%NB_TOT, PROCS, 1, rnba)
                        IF ((CCMP2) .and. ((NBI + NBI_START - 1) > VBMAX) .and. ((RNBA) > VBMAX)) CYCLE
                        ! GCHGIA number 2 for X-changed waves
                        ! calculate rho(r)=psi_i(r)* psi_a(r) for,
                        ! one center terms and, on the plane wave grid.
                        IF (ASSOCIATED(H)) THEN
                           CALL FOCK_CHARGE_ONE_CENTER_NOINT(WB(NBAA), WI(NBI), &
                                                             GWORK(1), H, CRHOIB(1, NBAA), CRHOLM, SIZE(CRHOLM))
                        ELSE
                           CALL FOCK_CHARGE_NOINT(WB(NBAA), WI(NBI), GWORK(1), &
                                                  CRHOLM, SIZE(CRHOLM))
                        END IF
                        ! Set phase e^iGr, where G = k_i - k_q - k_a
                        IF (LPHASE) THEN
                           CALL APPLY_PHASE(GRIDHF, CPHASE(1), GWORK(1), &
                                            GWORK(1))
                           !IF (KQ==1) WRITE(*,*)'error: no apply_phase need for kq=1'
                        END IF

                        ! FFT{rho} to reciprocal space
                        CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), &
                                    GWORK(1), GCHGIA(1, NBAA, 2), WGWQ%GRID, .FALSE.)
                        NFFT = NFFT + 1
                        ! multiply with potential factor
                     END DO !NBAA (loop over all bands in NSTRIP)

                     IF (ASSOCIATED(H)) THEN
                        CALL APPLY_PHASE_ONE_CENTER(WHF%WDES, H, CRHOIB(:, :), &
                                                    WHF%WDES%VKPT(:, KB) - WHF%WDES%VKPT(:, KI))
                     END IF

                     DO NBAA = 1, NSTRIPA
                        !copy co-densities to FTOD_PW and FTOD_OC

                        tmp_FTOD_PW(1:NP, NBAA + NBA - 1, NBI, 2) = (GCHGIA(1:NP, NBAA, 2))*(1.0_q/GRIDHF%NPLWV)
                        IF (ASSOCIATED(H)) THEN
                           tmp_FTOD_OC(1:H%TOTAL_ENTRIES, NBAA + NBA - 1, NBI, 2) = CRHOIB(1:H%TOTAL_ENTRIES, NBAA)
                        END IF

                        CALL LOC2GLOB(NBAA + NBA - 1, MYCOL, PROCS*WDES%NBANDS, PROCS, 1, rnba)
                     END DO
                  END DO !NBI (loop over nstrip bands only)
               END DO !NBA loop over all bands

               IF (EH_SCREENING) CALL REDISTRIBUTE_FTOD_S_GRID(W, WDES, KI, KQ, ISP, tmp_FTOD_PW_S, tmp_FTOD_OC, nbi_start, nbi_end)
               CALL REDISTRIBUTE_FTOD_GRID(W, WDES, KI, KQ, ISP, tmp_FTOD_PW, tmp_FTOD_OC, nbi_start, nbi_end)


            END DO !NBI_start loop over all bands
         END DO kiloop

         DEALLOCATE (GCHGIA, CRHOLM)
         IF (EH_SCREENING) DEALLOCATE (GCHGIA_S)
      END DO kqloop
      END DO spin

      IF (ALLOCATED(CRHOIA)) DEALLOCATE (CRHOIA)
      IF (ALLOCATED(CRHOIB)) DEALLOCATE (CRHOIB)
      IF (ALLOCATED(tmp_FTOD_OC)) DEALLOCATE (tmp_FTOD_OC)
      IF (ALLOCATED(tmp_FTOD_PW)) DEALLOCATE (tmp_FTOD_PW)
      IF (ALLOCATED(tmp_FTOD_PW_S)) DEALLOCATE (tmp_FTOD_PW_S)

      DO NBI = 1, NSTRIP
         CALL DELWAV(WI(NBI), .TRUE.)
      END DO
      DO NBA = 1, NSTRIP
         CALL DELWAV(WA(NBA), .TRUE.)
      END DO
      DO NBA = 1, NSTRIP
         CALL DELWAV(WB(NBA), .TRUE.)
      END DO
      DEALLOCATE (WI, WA, WB)

      IF (LSFACTOR) THEN
         FTOD_PW_AI_NOPOT = FTOD_PW_AI
         FTOD_PW_IJ_NOPOT = FTOD_PW_IJ
         CALL APPLY_POTFAK_ALL(WDES)
      END IF
      RETURN
   END SUBROUTINE CALC_FTOD

!***********************************************************************
!
! This subroutine redistributes the fourier-transformed overlap integrals <i|G|a>
! from the column-only process grid (CONTXT) to the quadratic
! process grid (CONTXT_GRID)
!
!***********************************************************************

   SUBROUTINE REDISTRIBUTE_FTOD_GRID(W, WDES, KI, KQ, ISP, tmp_FTOD_PW, tmp_FTOD_OC, nbi_start, nbi_end)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      INTEGER :: KA, KI, KQ, ISP, nbi_start, nbi_end
      GDEF, TARGET :: tmp_FTOD_OC(:, :, :, :)
      COMPLEX(q), TARGET :: tmp_FTOD_PW(:, :, :, :)
      INTEGER :: FTOD_PW_rows, FTOD_PW_cols, cc, FTOD_OC_rows, FTOD_OC_cols
      INTEGER :: NBI, FTOD_PW_rows_br, NBA

      KA = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - WDES%VKPT(:, KQ), KPOINTS_FULL)

      ! Prepare array descriptors for ScaLAPACK
      call BLACS_GRIDINFO(CONTXT_GRID, NPROW, NPCOL, MYROW, MYCOL)

      desc_FTOD_PW_br(1) = 1              ! descriptor type
      desc_FTOD_PW_br(2) = contxt_cols         ! blacs context
      desc_FTOD_PW_br(3) = NGVECTOR    ! global number of rows
      desc_FTOD_PW_br(4) = (PROCS*WDES%NBANDS) ! global number of cols
      desc_FTOD_PW_br(5) = NGVECTOR     ! row block size
      desc_FTOD_PW_br(6) = 1             ! col block size
      desc_FTOD_PW_br(7) = 0              ! initial process row
      desc_FTOD_PW_br(8) = 0              ! initial process col
      desc_FTOD_PW_br(9) = MAX(1, (NGVECTOR)) ! leading dimension of local array

      !Distribute the fourier-transformed overlap integrals
      desc_FTOD_PW = desc_FTOD_PW_br
      desc_FTOD_PW(2) = contxt_cols
      desc_FTOD_PW(6) = (PROCS*WDES%NBANDS)
      desc_FTOD_PW(8) = PROCS_KPTS(KI)


      DO NBI = nbi_start, nbi_end
         DO cc = 1, ncc
            CALL PZGEMR2D(NGVECTOR, (PROCS*WDES%NBANDS), tmp_FTOD_PW(1, 1, (NBI - nbi_start + 1), cc), 1, 1, &
                          desc_FTOD_PW_br, FTOD_PW(1, 1, cc), 1, 1, desc_FTOD_PW, contxt_grid)
            IF (ME == PROCS_KPTS(KI)) THEN
               DO NBA = 1, (PROCS*WDES%NBANDS)
                  IF (.not. LMETAL) THEN
                     IF (FILLED_MP2_ORBITAL(W%FERTOT(NBI, KI, 1))) THEN
                       IF (FILLED_MP2_ORBITAL(W%FERTOT(NBA, KA, 1))) FTOD_PW_IJ(:, NBA, NBI, RKQofKQ(KQ), MKPTS_KPTS(KI), cc) = &
                           (FTOD_PW(:, NBA, cc))
                IF (EMPTY_MP2_ORBITAL(W%FERTOT(NBA, KA, 1))) FTOD_PW_IA(:, NBA - VBMAX, NBI, RKQofKQ(KQ), MKPTS_KPTS(KI), cc) = &
                           (FTOD_PW(:, NBA, cc))
                     ELSE
               IF (FILLED_MP2_ORBITAL(W%FERTOT(NBA, KA, 1))) FTOD_PW_AI(:, NBA, NBI - VBMAX, RKQofKQ(KQ), MKPTS_KPTS(KI), cc) = &
                           (FTOD_PW(:, NBA, cc))
                        IF (CCMP2) CYCLE
        IF (EMPTY_MP2_ORBITAL(W%FERTOT(NBA, KA, 1))) FTOD_PW_AB(:, NBA - VBMAX, NBI - VBMAX, RKQofKQ(KQ), MKPTS_KPTS(KI), cc) = &
                           (FTOD_PW(:, NBA, cc))
                     END IF
                  ELSE
                     IF (VBMAX >= (NBI)) THEN
                        IF (VBMAX >= NBA) FTOD_PW_IJ(:, NBA, NBI, RKQofKQ(KQ), MKPTS_KPTS(KI), cc) = (FTOD_PW(:, NBA, cc))
                     END IF
                     IF (VBMAX >= NBI) FTOD_PW_IA(:, NBA, NBI, RKQofKQ(KQ), MKPTS_KPTS(KI), cc) = (FTOD_PW(:, NBA, cc))
                     IF (VBMAX >= NBA) FTOD_PW_AI(:, NBA, NBI, RKQofKQ(KQ), MKPTS_KPTS(KI), cc) = (FTOD_PW(:, NBA, cc))
                     IF (CCMP2) CYCLE
                     FTOD_PW_AB(:, NBA, NBI, RKQofKQ(KQ), MKPTS_KPTS(KI), cc) = (FTOD_PW(:, NBA, cc))
                  END IF
               END DO
            END IF
         END DO
      END DO

      !ONE-center part FTOD

      IF (ASSOCIATED(H)) THEN

         ! Prepare array descriptors for ScaLAPACK

         desc_FTOD_OC_br(1) = 1            ! descriptor type
         desc_FTOD_OC_br(2) = contxt_cols       ! blacs context
         desc_FTOD_OC_br(3) = NHVECTOR     ! global number of rows
         desc_FTOD_OC_br(4) = (PROCS*WDES%NBANDS) ! global number of cols
         desc_FTOD_OC_br(5) = NHVECTOR     ! row block size
         desc_FTOD_OC_br(6) = 1            ! col block size
         desc_FTOD_OC_br(7) = 0            ! initial process row
         desc_FTOD_OC_br(8) = 0            ! initial process col
         desc_FTOD_OC_br(9) = MAX(1, (NHVECTOR)) ! leading dimension of local array

         desc_FTOD_OC = desc_FTOD_OC_br
         desc_FTOD_OC(2) = contxt_cols
         desc_FTOD_OC(6) = (PROCS*WDES%NBANDS)
         desc_FTOD_OC(8) = PROCS_KPTS(KI)

         DO NBI = nbi_start, nbi_end
            DO cc = 1, 2
               CALL PZGEMR2D(NHVECTOR, (PROCS*WDES%NBANDS), tmp_FTOD_OC(1, 1, NBI - nbi_start + 1, cc), 1, 1, &
                             desc_FTOD_OC_br, FTOD_OC(1, 1, cc), 1, 1, desc_FTOD_OC, contxt_grid)

               IF (ME == PROCS_KPTS(KI)) THEN
                  DO NBA = 1, (PROCS*WDES%NBANDS)
                     IF (FILLED_MP2_ORBITAL(W%FERTOT(NBI, KI, 1))) THEN
                        IF (FILLED_MP2_ORBITAL(W%FERTOT(NBA, KA, 1))) FTOD_OC_IJ(:, NBA, NBI, RKQofKQ(KQ), MKPTS_KPTS(KI), cc) = &
                           FTOD_OC(:, NBA, cc)
                   IF (EMPTY_MP2_ORBITAL(W%FERTOT(NBA, KA, 1))) FTOD_OC_IA(:, NBA - VBMAX, NBI, RKQofKQ(KQ), MKPTS_KPTS(KI), cc) = &
                           FTOD_OC(:, NBA, cc)
                     ELSE
                  IF (FILLED_MP2_ORBITAL(W%FERTOT(NBA, KA, 1))) FTOD_OC_AI(:, NBA, NBI - VBMAX, RKQofKQ(KQ), MKPTS_KPTS(KI), cc) = &
                           FTOD_OC(:, NBA, cc)
           IF (EMPTY_MP2_ORBITAL(W%FERTOT(NBA, KA, 1))) FTOD_OC_AB(:, NBA - VBMAX, NBI - VBMAX, RKQofKQ(KQ), MKPTS_KPTS(KI), cc) = &
                           FTOD_OC(:, NBA, cc)
                     END IF
                  END DO
               END IF

            END DO
         END DO

      END IF

   END SUBROUTINE REDISTRIBUTE_FTOD_GRID

!***********************************************************************
!
! This subroutine redistributes the fourier-transformed overlap integrals <i|G|a>
! from the column-only process grid (CONTXT) to the quadratic
! process grid (CONTXT_GRID) using screened interactions (experimental feature)
!
!***********************************************************************

   SUBROUTINE REDISTRIBUTE_FTOD_S_GRID(W, WDES, KI, KQ, ISP, tmp_FTOD_PW_S, tmp_FTOD_OC_S, nbi_start, nbi_end)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      INTEGER :: KA, KI, KQ, ISP, nbi_start, nbi_end
      GDEF, TARGET :: tmp_FTOD_OC_S(:, :, :, :)
      COMPLEX(q), TARGET :: tmp_FTOD_PW_S(:, :, :, :)
      INTEGER :: FTOD_PW_rows, FTOD_PW_cols, cc, FTOD_OC_rows, FTOD_OC_cols
      INTEGER :: NBI, FTOD_PW_rows_br, NBA

      KA = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - WDES%VKPT(:, KQ), KPOINTS_FULL)

      ! Prepare array descriptors for ScaLAPACK
      call BLACS_GRIDINFO(CONTXT_GRID, NPROW, NPCOL, MYROW, MYCOL)

      desc_FTOD_PW_br(1) = 1              ! descriptor type
      desc_FTOD_PW_br(2) = contxt_cols         ! blacs context
      desc_FTOD_PW_br(3) = NGVECTOR    ! global number of rows
      desc_FTOD_PW_br(4) = (PROCS*WDES%NBANDS) ! global number of cols
      desc_FTOD_PW_br(5) = NGVECTOR     ! row block size
      desc_FTOD_PW_br(6) = 1             ! col block size
      desc_FTOD_PW_br(7) = 0              ! initial process row
      desc_FTOD_PW_br(8) = 0              ! initial process col
      desc_FTOD_PW_br(9) = MAX(1, (NGVECTOR)) ! leading dimension of local array
      !Distribute the fourier-transformed overlap integrals

      desc_FTOD_PW = desc_FTOD_PW_br
      desc_FTOD_PW(2) = contxt_cols
      desc_FTOD_PW(6) = (PROCS*WDES%NBANDS)
      desc_FTOD_PW(8) = PROCS_KPTS(KI)

      DO NBI = nbi_start, nbi_end
         DO cc = 1, 1
            CALL PZGEMR2D(NGVECTOR, (PROCS*WDES%NBANDS), tmp_FTOD_PW_S(1, 1, (NBI - nbi_start + 1), cc), 1, 1, &
                          desc_FTOD_PW_br, FTOD_PW(1, 1, cc), 1, 1, desc_FTOD_PW, contxt_grid)
            IF (ME == PROCS_KPTS(KI)) THEN
               DO NBA = 1, (PROCS*WDES%NBANDS)
                  IF (FILLED_MP2_ORBITAL(W%FERTOT(NBI, KI, 1))) THEN
                     IF (FILLED_MP2_ORBITAL(W%FERTOT(NBA, KA, 1))) FTOD_PW_S_IJ(:, NBA, NBI, KQ, MKPTS_KPTS(KI), cc) = &
                        (FTOD_PW(:, NBA, cc))
                     IF (EMPTY_MP2_ORBITAL(W%FERTOT(NBA, KA, 1))) FTOD_PW_S_IA(:, NBA - VBMAX, NBI, KQ, MKPTS_KPTS(KI), cc) = &
                        (FTOD_PW(:, NBA, cc))
                  ELSE
                     IF (FILLED_MP2_ORBITAL(W%FERTOT(NBA, KA, 1))) FTOD_PW_S_AI(:, NBA, NBI - VBMAX, KQ, MKPTS_KPTS(KI), cc) = &
                        (FTOD_PW(:, NBA, cc))
                     IF (CCMP2) CYCLE
               IF (EMPTY_MP2_ORBITAL(W%FERTOT(NBA, KA, 1))) FTOD_PW_S_AB(:, NBA - VBMAX, NBI - VBMAX, KQ, MKPTS_KPTS(KI), cc) = &
                        (FTOD_PW(:, NBA, cc))
                  END IF
               END DO
            END IF
         END DO
      END DO

   END SUBROUTINE REDISTRIBUTE_FTOD_S_GRID

!***********************************************************************
!
!  This routine initializes a process grid which is used for the
!  pre-calculation of the fourier-transformed overlap integrals <i|G|a>
!  This process grid has NPROCS(number of processors used) columns and 1 row.
!  The assigned context variable is called CONTXT_COLS.
!
!***********************************************************************

   SUBROUTINE INIT_BLACS_COLS()
      implicit none
      INTEGER :: a_PRCS, i
      REAL(q) :: NPCOL_TMP

      NPROW_GRID = 1
      !first we create a column-only process grid in order to
      !calculate the <i|-G|a> and <j|G|b> quantities
      call BLACS_PINFO(ME, PROCS)
      nprow = 1
      npcol = PROCS
      call BLACS_PINFO(ME, PROCS)
      call BLACS_GET(0, 0, CONTXT_COLS)

      call BLACS_GRIDINIT(CONTXT_COLS, 'R', NPROW, NPCOL)
      call BLACS_GRIDINFO(CONTXT_COLS, NPROW, NPCOL, MYROW, MYCOL)
      !since, a quadratic process grid is going to be needed for
      !the matrix-matrix multiplications, its optimal row and column
      !numbers are estimated for the given number of processors
      a_PRCS = CEILING(SQRT(Real(PROCS, kind=q)))
      IF (a_PRCS == SQRT(Real(PROCS, kind=q))) THEN
         NPROW_GRID = a_PRCS
         NPCOL_TMP = a_PRCS
      END IF
      IF (a_PRCS /= SQRT(Real(PROCS, kind=q))) THEN
         DO i = 1, CEILING(SQRT(Real(PROCS, kind=q)))
            NPCOL_TMP = PROCS/Real(i, kind=q)
            IF ((NPCOL_TMP - INT(NPCOL_TMP)) == 0) THEN
               NPROW_GRID = i
            END IF
         END DO
         NPCOL_TMP = PROCS/NPROW_GRID
      END IF

      IF (ABS(NPROW_GRID - a_PRCS) > (a_PRCS/2.0)) THEN
         WRITE (*, *) 'The allocated number of CPUs does not allow for the use of an efficient process grid.'
         WRITE (*, *) 'Suggested number of CPUs: 4, 16, 32, ....'
      END IF
   END SUBROUTINE INIT_BLACS_COLS

!***********************************************************************
!
! This routine initializes the process grid.
! The routine tries to create the most quadratic grid possible for the given
! number of processors.
!
!***********************************************************************

   SUBROUTINE INIT_BLACS_GRID(WDES)
      implicit none
      TYPE(wavedes) WDES

      call BLACS_PINFO(ME, PROCS)

      nprow = NPROW_GRID
      npcol = PROCS/NPROW
      IF (MOD(PROCS/NPROW, 1) /= 0) THEN
         WRITE (*, *) 'internal error in INIT_BLACS_GRID: Bad process grid'
      END IF

      IF (NGVECTOR > (PROCS*WDES%NBANDS)) THEN
         IF (NPROW_GRID < (PROCS/NPROW)) THEN
            NPCOL = NPROW_GRID
            NPROW = PROCS/NPROW_GRID
            NPROW_GRID = NPROW
         END IF
      END IF
      IF (NGVECTOR < (PROCS*WDES%NBANDS)) THEN
         IF (NPROW_GRID > (PROCS/NPROW)) THEN
            NPCOL = NPROW_GRID
            NPROW = PROCS/NPROW_GRID
            NPROW_GRID = NPROW
         END IF
      END IF
      IF ((MYROW == 0) .AND. (MYCOL == 0)) THEN
         WRITE (*, '(A,I3,A,I3,A)') 'The allocated processors form a', NPROW_GRID, 'x', NPCOL, ' grid.'
      END IF

      call BLACS_PINFO(ME, PROCS)
      call BLACS_GET(0, 0, CONTXT_GRID)
      call BLACS_GRIDINIT(CONTXT_GRID, 'R', NPROW, NPCOL)
      call BLACS_GRIDINFO(CONTXT_GRID, NPROW, NPCOL, MYROW, MYCOL)

   END SUBROUTINE INIT_BLACS_GRID

!***********************************************************************
!
!  Initialize the k-point distribution over MPI processes. Every process handles at least one k-point
!
!***********************************************************************

   SUBROUTINE INIT_MY_KPOINTS(WDES, WGW)
      implicit none
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      INTEGER :: MNK, ATPROCS

      call BLACS_PINFO(ME, PROCS)

      MY_NKPTS = 0
      ATPROCS = -1
      DO MNK = 1, WDES%NKPTS
         IF ((WDES%WTKPT(MNK) /= 0.0_q)) THEN
            ATPROCS = MOD((ATPROCS + 1), PROCS)
            IF (ATPROCS == ME) MY_NKPTS = MY_NKPTS + 1
         END IF
      END DO
      ALLOCATE (KPTS_MKPTS(MY_NKPTS))
      ALLOCATE (PROCS_KPTS(WDES%NKPTS))
      ALLOCATE (MKPTS_KPTS(WDES%NKPTS))

      MKPTS_KPTS(:) = 0
      MY_NKPTS = 0
      ATPROCS = -1
      DO MNK = 1, WDES%NKPTS
         IF ((WDES%WTKPT(MNK) /= 0.0_q)) THEN
            ATPROCS = MOD((ATPROCS + 1), PROCS)
            IF (ATPROCS == ME) MY_NKPTS = MY_NKPTS + 1
            IF (ATPROCS == ME) KPTS_MKPTS(MY_NKPTS) = MNK
            IF (ATPROCS == ME) MKPTS_KPTS(MNK) = MY_NKPTS
            PROCS_KPTS(MNK) = ATPROCS
         END IF
      END DO

   END SUBROUTINE INIT_MY_KPOINTS

!***********************************************************************
!
!  Here we allocate all important and required arrays for the amplitudes
!  and amplitude equations.
!
!***********************************************************************

   SUBROUTINE SETUP_T2AMPLITUDES(WDES)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER MB, NB, T2AMPLITUDES_ROWS, T2AMPLITUDES_COLS, NI, NA
      CALL BLACS_GRIDINFO(CONTXT_GRID, NPROW, NPCOL, MYROW, MYCOL)
      NPROW = NPROW_GRID
      NPCOL = PROCS/NPROW

      IF (ALLOCATED(T2)) THEN
         DEALLOCATE (T2)
      END IF
      IF (ALLOCATED(T2_N)) THEN
         DEALLOCATE (T2_N)
      END IF

      IF (.not. WDES%LORBITALREAL) THEN

         ALLOCATE (T2((NUNOCC), (NUNOCC), VBMAX, VBMAX, REALNKPTS, REALNKPTS, MY_NKPTS))
         ALLOCATE (CHI_KLIJ((VBMAX), (VBMAX), VBMAX, VBMAX, REALNKPTS, REALNKPTS, MY_NKPTS))
         ALLOCATE (CHI_KLIJ_TMP((VBMAX), (VBMAX), VBMAX, VBMAX, REALNKPTS, REALNKPTS))
         ALLOCATE (CHI_CKAI((NUNOCC), VBMAX, (NUNOCC), VBMAX, REALNKPTS, REALNKPTS, MY_NKPTS))
         ALLOCATE (T2_N((NUNOCC), (NUNOCC), VBMAX, VBMAX, REALNKPTS, REALNKPTS, MY_NKPTS))
         ALLOCATE (T2_TMP((NUNOCC), (NUNOCC), VBMAX, VBMAX, REALNKPTS, REALNKPTS))
         ALLOCATE (T2_KTMP((NUNOCC), (NUNOCC), VBMAX, VBMAX, REALNKPTS))
         ALLOCATE (T2_KTMP2((NUNOCC), VBMAX, (NUNOCC), VBMAX, REALNKPTS))
         IF (DAMPED) ALLOCATE (T2_VEL((NUNOCC), (NUNOCC), VBMAX, VBMAX, REALNKPTS, REALNKPTS, MY_NKPTS))
         IF (PRECONDITION) ALLOCATE (PRET2((NUNOCC), (NUNOCC), VBMAX, VBMAX, REALNKPTS, REALNKPTS, MY_NKPTS))
         IF (DAMPED) T2_VEL = (0.0_qs, 0.0_qs)
         IF (PRECONDITION) PRET2 = (0.0_qs, 0.0_qs)
         ALLOCATE (CHI_CKAI_KTMP((NUNOCC), VBMAX, (NUNOCC), VBMAX, REALNKPTS))
         ALLOCATE (CHI_CKAI_KTMP2((NUNOCC), VBMAX, (NUNOCC), VBMAX, REALNKPTS))
         ALLOCATE (K_KI(VBMAX, VBMAX, REALNKPTS))
         ALLOCATE (K_AC((NUNOCC), (NUNOCC), REALNKPTS))
         ALLOCATE (K_KC((VBMAX), (NUNOCC), REALNKPTS))
         ALLOCATE (L_KI(VBMAX, VBMAX, REALNKPTS))
         ALLOCATE (L_AC((NUNOCC), (NUNOCC), REALNKPTS))
         ALLOCATE (VV_S((NUNOCC), (NUNOCC)))
         ALLOCATE (VV2_S((NUNOCC), (NUNOCC)))
         ALLOCATE (VV((NUNOCC), (NUNOCC)))
         ALLOCATE (VV2((NUNOCC), (NUNOCC)))
         ALLOCATE (TE4O(VBMAX, VBMAX))
         ALLOCATE (T1((NUNOCC), VBMAX, REALNKPTS))
         ALLOCATE (T1_N((NUNOCC), VBMAX, REALNKPTS))
         ALLOCATE (T1_T((NUNOCC), VBMAX, REALNKPTS))
         ALLOCATE (VO_S((NUNOCC), VBMAX))
         ALLOCATE (OV_S(VBMAX, (NUNOCC)))
         ALLOCATE (VO2_S((NUNOCC), VBMAX))
         ALLOCATE (OV2_S(VBMAX, (NUNOCC)))
         IF (.not. CANONICAL) ALLOCATE (F_KI(VBMAX, VBMAX, REALNKPTS))
         IF (.not. CANONICAL) ALLOCATE (F_AI((NUNOCC), VBMAX, REALNKPTS))
         IF (.not. CANONICAL) ALLOCATE (F_KC(VBMAX, (NUNOCC), REALNKPTS))
         IF (.not. CANONICAL) ALLOCATE (F_BA((NUNOCC), (NUNOCC), REALNKPTS))

         ALLOCATE (VVOO((NUNOCC), (NUNOCC), VBMAX, VBMAX))
         ALLOCATE (OOVV(VBMAX, VBMAX, (NUNOCC), (NUNOCC)))
         ALLOCATE (OOVV_S(VBMAX, VBMAX, (NUNOCC), (NUNOCC)))
         ALLOCATE (OOVV2_S(VBMAX, VBMAX, (NUNOCC), (NUNOCC)))
         ALLOCATE (VVOO_S((NUNOCC), (NUNOCC), VBMAX, VBMAX))
         ALLOCATE (VVOO2_S((NUNOCC), (NUNOCC), VBMAX, VBMAX))
         ALLOCATE (VOVO_S((NUNOCC), VBMAX, (NUNOCC), VBMAX))
         ALLOCATE (VOVO((NUNOCC), VBMAX, (NUNOCC), VBMAX))
         ALLOCATE (VOVO2_S((NUNOCC), VBMAX, (NUNOCC), VBMAX))
         ALLOCATE (VOVO3_S((NUNOCC), VBMAX, (NUNOCC), VBMAX))
         ALLOCATE (OVOV(VBMAX, (NUNOCC), VBMAX, (NUNOCC)))
         ALLOCATE (OVVO(VBMAX, (NUNOCC), (NUNOCC), VBMAX))
         ALLOCATE (OVVV(VBMAX, (NUNOCC), (NUNOCC), (NUNOCC)))
         ALLOCATE (OVVV_S(VBMAX, (NUNOCC), (NUNOCC), (NUNOCC)))
         ALLOCATE (OOOV(VBMAX, VBMAX, VBMAX, (NUNOCC)))
         ALLOCATE (OOOV_S(VBMAX, VBMAX, VBMAX, (NUNOCC)))
         ALLOCATE (VOOV_S((NUNOCC), VBMAX, VBMAX, (NUNOCC)))
         ALLOCATE (OVVO_S(VBMAX, (NUNOCC), (NUNOCC), VBMAX))
         ALLOCATE (KLIJ(VBMAX, VBMAX, VBMAX, VBMAX))
         ALLOCATE (JLIK(VBMAX, VBMAX, VBMAX, VBMAX))
         IF (.not. CCMP2) ALLOCATE (VVVV((NUNOCC), (NUNOCC), NBLOCKAB, NBLOCKAB))
         IF (.not. CCMP2) ALLOCATE (VVVV_S((NUNOCC), (NUNOCC), NBLOCKAB, NBLOCKAB))

      ELSE

         ALLOCATE (T2_R((NUNOCC), (NUNOCC), VBMAX, VBMAX, WDES%NKPTS, WDES%NKPTS, MY_NKPTS))
         ALLOCATE (CHI_KLIJ_R((VBMAX), (VBMAX), VBMAX, VBMAX, WDES%NKPTS, WDES%NKPTS, MY_NKPTS))
         ALLOCATE (CHI_KLIJ_TMP_R((VBMAX), (VBMAX), VBMAX, VBMAX, WDES%NKPTS, WDES%NKPTS))
         ALLOCATE (CHI_CKAI_R((NUNOCC), VBMAX, (NUNOCC), VBMAX, WDES%NKPTS, WDES%NKPTS, MY_NKPTS))
         ALLOCATE (T2_N_R((NUNOCC), (NUNOCC), VBMAX, VBMAX, WDES%NKPTS, WDES%NKPTS, MY_NKPTS))
         ALLOCATE (T2_TMP_R((NUNOCC), (NUNOCC), VBMAX, VBMAX, WDES%NKPTS, WDES%NKPTS))
         ALLOCATE (T2_KTMP_R((NUNOCC), (NUNOCC), VBMAX, VBMAX, WDES%NKPTS))
         ALLOCATE (T2_KTMP2_R((NUNOCC), VBMAX, (NUNOCC), VBMAX, WDES%NKPTS))
         IF (DAMPED) ALLOCATE (T2_VEL_R((NUNOCC), (NUNOCC), VBMAX, VBMAX, WDES%NKPTS, WDES%NKPTS, MY_NKPTS))
         IF (PRECONDITION) ALLOCATE (PRET2_R((NUNOCC), (NUNOCC), VBMAX, VBMAX, WDES%NKPTS, WDES%NKPTS, MY_NKPTS))
         IF (DAMPED) T2_VEL = (0.0_qs)
         IF (PRECONDITION) PRET2 = (0.0_qs)
         ALLOCATE (CHI_CKAI_KTMP_R((NUNOCC), VBMAX, (NUNOCC), VBMAX, WDES%NKPTS))
         ALLOCATE (CHI_CKAI_KTMP2_R((NUNOCC), VBMAX, (NUNOCC), VBMAX, WDES%NKPTS))
         ALLOCATE (K_KI_R(VBMAX, VBMAX, WDES%NKPTS))
         ALLOCATE (K_AC_R((NUNOCC), (NUNOCC), WDES%NKPTS))
         ALLOCATE (K_KC_R((VBMAX), (NUNOCC), WDES%NKPTS))
         ALLOCATE (L_KI_R(VBMAX, VBMAX, WDES%NKPTS))
         ALLOCATE (L_AC_R((NUNOCC), (NUNOCC), WDES%NKPTS))
         ALLOCATE (VV_S_R((NUNOCC), (NUNOCC)))
         ALLOCATE (VV2_S_R((NUNOCC), (NUNOCC)))
         ALLOCATE (VV_R((NUNOCC), (NUNOCC)))
         ALLOCATE (VV2_R((NUNOCC), (NUNOCC)))
         ALLOCATE (TE4O_R(VBMAX, VBMAX))
         ALLOCATE (T1_R((NUNOCC), VBMAX, WDES%NKPTS))
         ALLOCATE (T1_N_R((NUNOCC), VBMAX, WDES%NKPTS))
         ALLOCATE (T1_T_R((NUNOCC), VBMAX, WDES%NKPTS))
         ALLOCATE (VO_S_R((NUNOCC), VBMAX))
         ALLOCATE (OV_S_R(VBMAX, (NUNOCC)))
         ALLOCATE (VO2_S_R((NUNOCC), VBMAX))
         ALLOCATE (OV2_S_R(VBMAX, (NUNOCC)))
         IF (.not. CANONICAL) ALLOCATE (F_KI_R(VBMAX, VBMAX, WDES%NKPTS))
         IF (.not. CANONICAL) ALLOCATE (F_AI_R((NUNOCC), VBMAX, WDES%NKPTS))
         IF (.not. CANONICAL) ALLOCATE (F_KC_R(VBMAX, (NUNOCC), WDES%NKPTS))
         IF (.not. CANONICAL) ALLOCATE (F_BA_R((NUNOCC), (NUNOCC), WDES%NKPTS))

         ALLOCATE (VVOO_R((NUNOCC), (NUNOCC), VBMAX, VBMAX))
         ALLOCATE (OOVV_R(VBMAX, VBMAX, (NUNOCC), (NUNOCC)))
         ALLOCATE (OOVV_S_R(VBMAX, VBMAX, (NUNOCC), (NUNOCC)))
         ALLOCATE (OOVV2_S_R(VBMAX, VBMAX, (NUNOCC), (NUNOCC)))
         ALLOCATE (VVOO_S_R((NUNOCC), (NUNOCC), VBMAX, VBMAX))
         ALLOCATE (VVOO2_S_R((NUNOCC), (NUNOCC), VBMAX, VBMAX))
         ALLOCATE (VOVO_S_R((NUNOCC), VBMAX, (NUNOCC), VBMAX))
         ALLOCATE (VOVO_R((NUNOCC), VBMAX, (NUNOCC), VBMAX))
         ALLOCATE (VOVO2_S_R((NUNOCC), VBMAX, (NUNOCC), VBMAX))
         ALLOCATE (VOVO3_S_R((NUNOCC), VBMAX, (NUNOCC), VBMAX))
         ALLOCATE (OVOV_R(VBMAX, (NUNOCC), VBMAX, (NUNOCC)))
         ALLOCATE (OVVO_R(VBMAX, (NUNOCC), (NUNOCC), VBMAX))
         ALLOCATE (OVVV_R(VBMAX, (NUNOCC), (NUNOCC), (NUNOCC)))
         ALLOCATE (OVVV_S_R(VBMAX, (NUNOCC), (NUNOCC), (NUNOCC)))
         ALLOCATE (OOOV_R(VBMAX, VBMAX, VBMAX, (NUNOCC)))
         ALLOCATE (OOOV_S_R(VBMAX, VBMAX, VBMAX, (NUNOCC)))
         ALLOCATE (VOOV_S_R((NUNOCC), VBMAX, VBMAX, (NUNOCC)))
         ALLOCATE (OVVO_S_R(VBMAX, (NUNOCC), (NUNOCC), VBMAX))
         ALLOCATE (KLIJ_R(VBMAX, VBMAX, VBMAX, VBMAX))
         ALLOCATE (JLIK_R(VBMAX, VBMAX, VBMAX, VBMAX))
         IF (.not. CCMP2) ALLOCATE (VVVV_R((NUNOCC), (NUNOCC), NBLOCKAB, NBLOCKAB))
         IF (.not. CCMP2) ALLOCATE (VVVV_S_R((NUNOCC), (NUNOCC), NBLOCKAB, NBLOCKAB))

      END IF

      IF (LPSPPL) THEN
         ALLOCATE (D2_OOOO(VBMAX, VBMAX, REALNKPTS, REALNKPTS))
         ALLOCATE (D2PAW_OOOO(VBMAX, VBMAX, REALNKPTS, REALNKPTS))
         ALLOCATE (D2PAW_VVOO(NUNOCC,NUNOCC,VBMAX, VBMAX, REALNKPTS, REALNKPTS, MY_NKPTS))
         ALLOCATE (EMP2_PAIR_CBS(VBMAX, VBMAX, REALNKPTS, REALNKPTS))
         ALLOCATE (EMP2_PAIR_FBS(VBMAX, VBMAX, REALNKPTS, REALNKPTS))
         ALLOCATE (GMP2_PAIR(VBMAX, VBMAX, REALNKPTS, REALNKPTS))
         ALLOCATE (GCC_PAIR(VBMAX, VBMAX, REALNKPTS, REALNKPTS))
      ENDIF

      ALLOCATE (NI_N(VBMAX), NA_N(NUNOCC))
      DO NI = 1, VBMAX
         NI_N(NI) = NI
      END DO

      DO NA = 1, NUNOCC
         NA_N(NA) = NA + VBMAX
         IF (LMETAL) NA_N(NA) = NA
      END DO

   END SUBROUTINE SETUP_T2AMPLITUDES

!***********************************************************************
!
!  Allocate the arrays needed for the Coulomb-Vertex representations
!
!***********************************************************************

   SUBROUTINE SETUP_FTOD(WDES)
      IMPLICIT NONE
      TYPE(wavedes) WDES

      call BLACS_GRIDINFO(CONTXT_GRID, NPROW, NPCOL, MYROW, MYCOL)
      !Blocking size for block-cyclic distribution of RESPF_PW

      desc_FTOD_PW(1) = 1              ! descriptor type
      desc_FTOD_PW(2) = CONTXT_GRID      ! blacs context
      desc_FTOD_PW(3) = NGVECTOR       ! global number of rows
      desc_FTOD_PW(4) = (PROCS*WDES%NBANDS)   ! global number of cols
      desc_FTOD_PW(5) = NGVECTOR       ! row block size
      desc_FTOD_PW(6) = (PROCS*WDES%NBANDS)       ! column block size
      desc_FTOD_PW(7) = 0              ! initial process row
      desc_FTOD_PW(8) = 0              ! initial process column
      desc_FTOD_PW(9) = MAX(1, NGVECTOR) ! leading dimension of local array

      ALLOCATE (FTOD_PW(NGVECTOR, (PROCS*WDES%NBANDS), ncc))
      ALLOCATE (FTOD_PW_IA(NGVECTOR, (NUNOCC), VBMAX, REALNKPTS, MY_NKPTS, ncc))
      IF (EH_SCREENING) ALLOCATE (FTOD_PW_S_IA(NGVECTOR, (NUNOCC), VBMAX, REALNKPTS, MY_NKPTS, 1))
      ALLOCATE (FTOD_PW_AI(NGVECTOR, VBMAX, (NUNOCC), REALNKPTS, MY_NKPTS, ncc))
      IF (LSFACTOR) ALLOCATE (FTOD_PW_AI_NOPOT(NGVECTOR, VBMAX, (NUNOCC), REALNKPTS, MY_NKPTS, ncc))
      IF (LSFACTOR) ALLOCATE (FTOD_PW_IJ_NOPOT(NGVECTOR, VBMAX, VBMAX, REALNKPTS, MY_NKPTS, ncc))
      IF (EH_SCREENING) ALLOCATE (FTOD_PW_S_AI(NGVECTOR, VBMAX, (NUNOCC), REALNKPTS, MY_NKPTS, 1))
      IF (.not. CCMP2) ALLOCATE (FTOD_PW_AB(NGVECTOR, (NUNOCC), (NUNOCC), REALNKPTS, MY_NKPTS, ncc))
      IF (EH_SCREENING) ALLOCATE (FTOD_PW_S_AB(NGVECTOR, (NUNOCC), (NUNOCC), REALNKPTS, MY_NKPTS, 1))
      ALLOCATE (FTOD_PW_IJ(NGVECTOR, VBMAX, VBMAX, REALNKPTS, MY_NKPTS, ncc))
      IF (EH_SCREENING) ALLOCATE (FTOD_PW_S_IJ(NGVECTOR, VBMAX, VBMAX, REALNKPTS, MY_NKPTS, 1))

      ALLOCATE (PW_AI_TMP(NGVECTOR, VBMAX, (NUNOCC), REALNKPTS))
      ALLOCATE (PW_IA_TMP(NGVECTOR, (NUNOCC), VBMAX, REALNKPTS))
      ALLOCATE (PW_IA_TMP2(NGVECTOR, (NUNOCC), VBMAX, REALNKPTS))
      ALLOCATE (PW_IJ_TMP(NGVECTOR, VBMAX, VBMAX, REALNKPTS))
      IF (.not. CCMP2) ALLOCATE (PW_AB_TMP(NGVECTOR, (NUNOCC), (NUNOCC), REALNKPTS))

      IF (ASSOCIATED(H)) THEN
         !Blocking size for block-cyclic distribution of RESPF_PW

         desc_FTOD_OC(1) = 1              ! descriptor type
         desc_FTOD_OC(2) = CONTXT_GRID        ! blacs context
         desc_FTOD_OC(3) = NHVECTOR       ! global number of rows
         desc_FTOD_OC(4) = (PROCS*WDES%NBANDS)   ! global number of cols
         desc_FTOD_OC(5) = NHVECTOR             ! row block size
         desc_FTOD_OC(6) = (PROCS*WDES%NBANDS)             ! column block size
         desc_FTOD_OC(7) = 0              ! initial process row
         desc_FTOD_OC(8) = 0              ! initial process column
         desc_FTOD_OC(9) = MAX(1, NHVECTOR) ! leading dimension of local array

         ALLOCATE (FTOD_OC(NHVECTOR, (PROCS*WDES%NBANDS), 2))
         ALLOCATE (FTOD_OC_IA(NHVECTOR, (NUNOCC), VBMAX, REALNKPTS, MY_NKPTS, 2))
         ALLOCATE (FTOD_OC_AI(NHVECTOR, VBMAX, (NUNOCC), REALNKPTS, MY_NKPTS, 2))
         IF (.not. CCMP2) ALLOCATE (FTOD_OC_AB(NHVECTOR, (NUNOCC), (NUNOCC), REALNKPTS, MY_NKPTS, 2))
         ALLOCATE (FTOD_OC_IJ(NHVECTOR, VBMAX, VBMAX, REALNKPTS, MY_NKPTS, 2))
         ALLOCATE (OC_IJ_TMP(NHVECTOR, VBMAX, VBMAX, REALNKPTS))
         ALLOCATE (OC_IA_TMP(NHVECTOR, (NUNOCC), VBMAX, REALNKPTS), OC_IA_TMP2(NHVECTOR, (NUNOCC), VBMAX, REALNKPTS))
         ALLOCATE (OC_AI_TMP(NHVECTOR, VBMAX, (NUNOCC), REALNKPTS))
         IF (.not. CCMP2) ALLOCATE (OC_AB_TMP(NHVECTOR, (NUNOCC), (NUNOCC), REALNKPTS))
      END IF

   END SUBROUTINE SETUP_FTOD

!***********************************************************************
!
! TRANSFORM LOCAL TO GLOBAL ARRAY INDEX FOR BLOCK-CYCLIC ARRAY
! DISTRIBUTION
!
!***********************************************************************

   SUBROUTINE LOC2GLOB(li, p, n, np, nb, gi)
      IMPLICIT NONE
      INTEGER :: li   ! local index
      INTEGER :: p    ! index in processor grid (either MYROW or MYCOL)
      INTEGER :: n    ! global array dimension
      INTEGER :: np   ! processor array dimension
      INTEGER :: nb   ! blocking size
      INTEGER :: gi   ! global index
      INTEGER :: litmp

      litmp = li - 1
      gi = (((litmp/nb)*np) + p)*nb + mod(litmp, nb) + 1
      RETURN
   END SUBROUTINE LOC2GLOB

!***********************************************************************
!
! Unsuccessfull experimental feature for finite size corrections
!
!***********************************************************************

   subroutine APPEND_CDER_TO_FTOD(WDES, W, LATT_CUR, dir)
      USE constant
      USE full_kpoints
      USE mkpoints
      USE wave
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      TYPE(latt) LATT_CUR
      integer :: dir !xyz direction of the dipolemoments
      integer :: mni, mna, mki, mkq, mng, mrng, mrna, mrni
      INTEGER :: FTOD_ROWS, TWOE4ORBITAL_COLS
      GDEF :: CDER_BETWEEN_STATE_IA(3), CDER_BETWEEN_STATE_AI(3)
      real(qs) :: ediff_th

      do mki = 1, MY_NKPTS

         MKQ = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:, KPTS_MKPTS(mki)) - W%WDES%VKPT(:, KPTS_MKPTS(mki)), KPOINTS_FULL) !just to find the gamma-point
         do mni = 1, VBMAX
            do mna = 1, (NUNOCC)

               CDER_BETWEEN_STATE_IA(:) = zero
               CDER_BETWEEN_STATE_AI(:) = zero

               CALL CDER_BETWEEN_STATES_ROTATED( &
                  CDER_BETWEEN_STATE_IA, LATT_CUR, KPTS_MKPTS(mki), 1, mna + VBMAX, mni)

               FTOD_PW_IA(ngvector, mna, mni, mkq, mki, 1) = (CDER_BETWEEN_STATE_IA(dir)*EDEPS/LATT_CUR%OMEGA)
               FTOD_PW_IA(ngvector, mna, mni, mkq, mki, 2) = conjg(CDER_BETWEEN_STATE_IA(dir))

               FTOD_PW_AI(ngvector, mni, mna, mkq, mki, 1) = conjg((CDER_BETWEEN_STATE_IA(dir)*EDEPS/LATT_CUR%OMEGA))
               FTOD_PW_AI(ngvector, mni, mna, mkq, mki, 2) = (CDER_BETWEEN_STATE_IA(dir))

            end do
         end do
      end do

   end subroutine APPEND_CDER_TO_FTOD

!***********************************************************************
!
!  This routine computes the Fock exchange energy.
!  This is a good sanity check of many quantities used by this module.
!
!***********************************************************************

   Subroutine CALC_EFOCK(WDES, WGW, IO)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      TYPE(in_struct) IO
      INTEGER :: KQ, KI, KJ, KQ_, NI, NJ
      write (*, *) 'me,my_kpts', me
      IF (WDES%LORBITALREAL) THEN
         TE4O = zero
      ELSE
         TE4O_R = (0.0_q)
      END IF
      call BLACS_GRIDINFO(CONTXT_COLS, NPROW, NPCOL, MYROW, MYCOL)
      DO KI = 1, WDES%NKPTS
         CALL BCAST2ALL_FTOD_IJ(WDES, KI, 1, PW_IJ_TMP, OC_IJ_TMP)

         DO KJ = 1, MY_NKPTS
         DO NI = 1, VBMAX
         DO NJ = 1, VBMAX
            KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - WDES%VKPT(:, KPTS_MKPTS(KJ)), KPOINTS_FULL)
            CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                            PW_IJ_TMP(1, NJ, NI, RKQofKQ(KQ)), OC_IJ_TMP(1, NJ, NI, RKQofKQ(KQ)), &
                            (1)*1, &
                           FTOD_PW_IJ(1, NI, NJ, RKQofKQ(KQ), RKIofKI(KJ), 2), FTOD_OC_IJ(1, NI, NJ, RKQofKQ(KQ), RKIofKI(KJ), 2), &
                            (1)*1, TE4O(1, 1), TE4O_R(1, 1), (1.0_q, 0.0_q))

         END DO
         END DO
         END DO
      END DO

      IF (.not. WDES%LORBITALREAL) THEN
         TE4O(1, 1) = GCONJG(TE4O(1, 1)*WTKPT*WTKPT)
         CALLMPI(M_sum_z(WGW%COMM_INTER, TE4O(1, 1), 1))
         IF (IO%IU0 >= 0) write (*, *) 'EFOCK', TE4O(1, 1)
         IF (IO%IU0 >= 0) write (*, *) 'VBMAX=', VBMAX
      ELSE
         TE4O_R(1, 1) = (TE4O_R(1, 1)*WTKPT*WTKPT)
         CALLMPI(M_sum_d(WGW%COMM_INTER, TE4O_R(1, 1), 1))
         IF (IO%IU0 >= 0) write (*, *) 'EFOCK', TE4O_R(1, 1)
         IF (IO%IU0 >= 0) write (*, *) 'VBMAX=', VBMAX
      END IF

   END SUBROUTINE CALC_EFOCK

!***********************************************************************
!
!  This routine computes the CCSD correlation energy
!
!***********************************************************************

   SUBROUTINE CALC_CCSD_ENERGY(WDES, WGW, W, IO, ITERATION)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      TYPE(wavespin) W
      TYPE(in_struct) IO
      INTEGER :: KQ, KI, KJ, KQ_, NI, NJ, KA, KB, NA, NB, ITERATION
      GDEFS :: ETMP
      GDEF :: EPAIR

      E_CCSD = zero
      E_CCSD_X = zero
      ETMP = zero

      IF ((ITERATION==1) .and. (LPSPPL)) THEN
         EMP2_PAIR_FBS(:,:,:,:)=zero
      ENDIF



      DO KB = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE

         IF (.not. CANONICAL) THEN
            IF (.not. CCMP2) THEN
            IF (WDES%LORBITALREAL) THEN

               OV2_S_R = F_KC_R(:, :, RKIofKI(KB))
               CALL SORT_O1V1_V1O1(WDES, OV2_S, VO_S, OV2_S_R, VO_S_R)
               CALL SGEMM('T', 'n', 1, 1, &
                          (VBMAX)*(NUNOCC), (2._qs, 0._qs), VO_S_R(1, 1), &
                          (VBMAX)*(NUNOCC), &
                          T1_R(1, 1, RKIofKI(KB)), (VBMAX)*(NUNOCC), &
                          (1._qs, 0._qs), ETMP, 1)
            ELSE
               OV2_S = (F_KC(:, :, RKIofKI(KB)))
               CALL SORT_O1V1_V1O1(WDES, OV2_S, VO_S, OV2_S_R, VO_S_R)
               CALL CGEMM('t', 'n', 1, 1, &
                          (VBMAX)*(NUNOCC), (2._qs, 0._qs), VO_S(1, 1), &
                          (VBMAX)*(NUNOCC), &
                          T1(1, 1, RKIofKI(KB)), (VBMAX)*(NUNOCC), &
                          (1._qs, 0._qs), ETMP, 1)
            END IF
            END IF
         END IF

         CALL BCAST2ALL_FTOD_AI(WDES, KB, 2, PW_AI_TMP, OC_AI_TMP)
         DO KA = 1, MY_NKPTS
         DO KQ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.00_q)) CYCLE
            IF (WDES%LORBITALREAL) THEN
               OOVV_R = (0.0_q)
            ELSE
               OOVV = zero
            END IF
            KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
            KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) + &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
            KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                      WDES%VKPT(:, KJ), KPOINTS_FULL)
            DO NA = 1, (NUNOCC)
            DO NB = 1, (NUNOCC)
               CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                               PW_AI_TMP(1, 1, NB, RKQofKQ(KQ_)), OC_AI_TMP(1, 1, NB, RKQofKQ(KQ_)), &
                               VBMAX, &
                               FTOD_PW_AI(1, 1, NA, RKQofKQ(KQ_), KA, 1), FTOD_OC_AI(1, 1, NA, RKQofKQ(KQ_), KA, 1), &
                               VBMAX, OOVV(1, 1, NA, NB), OOVV_R(1, 1, NA, NB), zero)
            END DO
            END DO

            DO NI = 1, VBMAX
            DO NJ = 1, VBMAX
               DO NB = 1, (NUNOCC)
               DO NA = 1, (NUNOCC)
                  IF (WDES%LORBITALREAL) THEN
                     E_CCSD_X = E_CCSD_X + (OOVV_R(NI, NJ, NA, NB))* &
                                (T2_R(NA, NB, NI, NJ, KI, KJ, KA))*KPOINTS_FULL%WTKPT(1) &
                                *KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1)
                  ELSE
                     E_CCSD_X = E_CCSD_X + (OOVV(NI, NJ, NA, NB))* &
                                (T2(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA))*WTKPT &
                                *WTKPT*WTKPT
                  END IF

               END DO
               END DO
            END DO
            END DO

            IF (SINGLES .and. (KI == KPTS_MKPTS(KA))) THEN
               DO NI = 1, VBMAX
               DO NJ = 1, VBMAX
                  DO NB = 1, (NUNOCC)
                  DO NA = 1, (NUNOCC)
                  IF (WDES%LORBITALREAL) THEN
                     E_CCSD_X = E_CCSD_X + (OOVV_R(NI, NJ, NA, NB))* &
                                (T1_R(NA, NI, KI)*T1_R(NB, NJ, KJ))*WTKPT &
                                *WTKPT
                  ELSE
                     E_CCSD_X = E_CCSD_X + (OOVV(NI, NJ, NA, NB))* &
                                (T1(NA, NI, RKIofKI(KI))*T1(NB, NJ, RKIofKI(KJ)))*WTKPT &
                                *WTKPT
                  END IF
                  END DO
                  END DO
               END DO
               END DO
            END IF

            IF (WDES%LORBITALREAL) THEN
               OOVV_R = -(0.5_qs)*OOVV_R
            ELSE
               OOVV = -CONJG((0.5_qs, 0.0_qs)*OOVV)
            END IF

            DO NA = 1, (NUNOCC)
            DO NB = 1, (NUNOCC)
               CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                               FTOD_PW_AI(1, 1, NA, RKQofKQ(KQ), KA, 1), FTOD_OC_AI(1, 1, NA, RKQofKQ(KQ), KA, 1), &
                               VBMAX, &
                               PW_AI_TMP(1, 1, NB, RKQofKQ(KQ)), OC_AI_TMP(1, 1, NB, RKQofKQ(KQ)), &
                               VBMAX, OOVV(1, 1, NA, NB), OOVV_R(1, 1, NA, NB), one)
            END DO
            END DO

            IF (WDES%LORBITALREAL) THEN
               OOVV_R = ((2.0_qs)*OOVV_R)
            ELSE
               OOVV = CONJG((2.0_qs, 0.0_qs)*OOVV)
            END IF


            DO NI = 1, VBMAX
            DO NJ = 1, VBMAX
               DO NB = 1, (NUNOCC)
               DO NA = 1, (NUNOCC)

                  IF (WDES%LORBITALREAL) THEN
                     E_CCSD = E_CCSD + OOVV_R(NI, NJ, NA, NB)* &
                              (T2_R(NA, NB, NI, NJ, KI, KJ, KA))*KPOINTS_FULL%WTKPT(1) &
                              *KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1)
                  ELSE
                     E_CCSD = E_CCSD + OOVV(NI, NJ, NA, NB)* &
                              (T2(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA))*WTKPT &
                              *WTKPT*WTKPT
                     IF ((ITERATION==1) .and. (LPSPPL)) THEN
                        EMP2_PAIR_FBS(NI,NJ,RKIofKI(KI), RKIofKI(KJ))=&
                           EMP2_PAIR_FBS(NI,NJ,RKIofKI(KI), RKIofKI(KJ))+&
                           OOVV(NI, NJ, NA, NB)* &
                           (T2(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA))*WTKPT &
                           *WTKPT*WTKPT
                     ENDIF

                  END IF

               END DO
               END DO
            END DO
            END DO

            IF (SINGLES .and. (KI == KPTS_MKPTS(KA))) THEN
               DO NI = 1, VBMAX
               DO NJ = 1, VBMAX
                  DO NB = 1, (NUNOCC)
                  DO NA = 1, (NUNOCC)
                     IF (WDES%LORBITALREAL) THEN
                        E_CCSD = E_CCSD + OOVV_R(NI, NJ, NA, NB)* &
                                 (T1_R(NA, NI, KI)*T1_R(NB, NJ, KJ))*KPOINTS_FULL%WTKPT(1) &
                                 *KPOINTS_FULL%WTKPT(1)
                     ELSE
                        E_CCSD = E_CCSD + OOVV(NI, NJ, NA, NB)* &
                                 (T1(NA, NI, RKIofKI(KI))*T1(NB, NJ, RKIofKI(KJ)))*WTKPT &
                                 *WTKPT
                     END IF
                  END DO
                  END DO
               END DO
               END DO
            END IF
         END DO
         END DO
      END DO
      CALL M_sum_z(WGW%COMM_INTER, E_CCSD, 1)
      CALL M_sum_z(WGW%COMM_INTER, E_CCSD_X, 1)
      E_CCSD = E_CCSD + ETMP*WTKPT
      E_CCSD_X = E_CCSD_X + ETMP*WTKPT
      IF (IO%IU0 >= 0) WRITE (IO%IU0, *) 'E_CCSD=', E_CCSD
      IF (IO%IU0 >= 0) WRITE (IO%IU6, *) 'E CCSD=', E_CCSD
      IF (IO%IU0 >= 0) WRITE (IO%IU0, *) 'E_CCSD_X=', E_CCSD_X
      IF (IO%IU0 >= 0) WRITE (IO%IU6, *) 'E CCSD_X=', E_CCSD_X
      IF (ITERATION == 1) E_MP2 = E_CCSD

      IF ((ITERATION==1) .and. (LPSPPL)) THEN
         CALLMPI(M_sum_z(WDES%COMM, EMP2_PAIR_FBS(1,1,1,1), VBMAX**2*REALNKPTS**2 ))
      ENDIF

   END SUBROUTINE CALC_CCSD_ENERGY

!***********************************************************************
!
!  This routine computes the PPL contribution to the CCSD correlation energy
!
!***********************************************************************

   SUBROUTINE CALC_CCSD_ENERGY_PPL(WDES, WGW, W, IO, ITERATION)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      TYPE(wavespin) W
      TYPE(in_struct) IO
      INTEGER :: KQ, KI, KJ, KQ_, NI, NJ, KA, KB, NA, NB, ITERATION
      GDEFS :: ETMP

      E_CCSD = zero
      E_CCSD_X = zero
      ETMP = zero

      DO KB = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE

         IF (.not. CANONICAL) THEN
            IF (.not. CCMP2) THEN
            IF (WDES%LORBITALREAL) THEN

               OV2_S_R = F_KC_R(:, :, RKIofKI(KB))
               CALL SORT_O1V1_V1O1(WDES, OV2_S, VO_S, OV2_S_R, VO_S_R)
               CALL SGEMM('T', 'n', 1, 1, &
                          (VBMAX)*(NUNOCC), (2._qs, 0._qs), VO_S_R(1, 1), &
                          (VBMAX)*(NUNOCC), &
                          T1_R(1, 1, RKIofKI(KB)), (VBMAX)*(NUNOCC), &
                          (1._qs, 0._qs), ETMP, 1)
            ELSE
               OV2_S = (F_KC(:, :, RKIofKI(KB)))
               CALL SORT_O1V1_V1O1(WDES, OV2_S, VO_S, OV2_S_R, VO_S_R)
               CALL CGEMM('t', 'n', 1, 1, &
                          (VBMAX)*(NUNOCC), (2._qs, 0._qs), VO_S(1, 1), &
                          (VBMAX)*(NUNOCC), &
                          T1(1, 1, RKIofKI(KB)), (VBMAX)*(NUNOCC), &
                          (1._qs, 0._qs), ETMP, 1)
            END IF
            END IF
         END IF

         CALL BCAST2ALL_FTOD_AI(WDES, KB, 2, PW_AI_TMP, OC_AI_TMP)
         DO KA = 1, MY_NKPTS
         DO KQ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.00_q)) CYCLE
            IF (WDES%LORBITALREAL) THEN
               OOVV_R = (0.0_q)
            ELSE
               OOVV = zero
            END IF
            KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
            KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) + &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
            KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                      WDES%VKPT(:, KJ), KPOINTS_FULL)
            DO NA = 1, (NUNOCC)
            DO NB = 1, (NUNOCC)
               CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                               PW_AI_TMP(1, 1, NB, RKQofKQ(KQ_)), OC_AI_TMP(1, 1, NB, RKQofKQ(KQ_)), &
                               VBMAX, &
                               FTOD_PW_AI(1, 1, NA, RKQofKQ(KQ_), KA, 1), FTOD_OC_AI(1, 1, NA, RKQofKQ(KQ_), KA, 1), &
                               VBMAX, OOVV(1, 1, NA, NB), OOVV_R(1, 1, NA, NB), zero)
            END DO
            END DO

            DO NI = 1, VBMAX
            DO NJ = 1, VBMAX
               DO NB = 1, (NUNOCC)
               DO NA = 1, (NUNOCC)
                  IF (WDES%LORBITALREAL) THEN
                     E_CCSD_X = E_CCSD_X + (OOVV_R(NI, NJ, NA, NB))* &
                                (T2_N_R(NA, NB, NI, NJ, KI, KJ, KA))*KPOINTS_FULL%WTKPT(1) &
                                *KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1)
                  ELSE
                     E_CCSD_X = E_CCSD_X + (OOVV(NI, NJ, NA, NB))* &
                                (T2_N(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA))*WTKPT &
                                *WTKPT*WTKPT
                  END IF

               END DO
               END DO
            END DO
            END DO

            IF (WDES%LORBITALREAL) THEN
               OOVV_R = -(0.5_qs)*OOVV_R
            ELSE
               OOVV = -CONJG((0.5_qs, 0.0_qs)*OOVV)
            END IF

            DO NA = 1, (NUNOCC)
            DO NB = 1, (NUNOCC)
               CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                               FTOD_PW_AI(1, 1, NA, RKQofKQ(KQ), KA, 1), FTOD_OC_AI(1, 1, NA, RKQofKQ(KQ), KA, 1), &
                               VBMAX, &
                               PW_AI_TMP(1, 1, NB, RKQofKQ(KQ)), OC_AI_TMP(1, 1, NB, RKQofKQ(KQ)), &
                               VBMAX, OOVV(1, 1, NA, NB), OOVV_R(1, 1, NA, NB), one)
            END DO
            END DO

            IF (WDES%LORBITALREAL) THEN
               OOVV_R = ((2.0_qs)*OOVV_R)
            ELSE
               OOVV = CONJG((2.0_qs, 0.0_qs)*OOVV)
            END IF

            DO NI = 1, VBMAX
            DO NJ = 1, VBMAX
               DO NB = 1, (NUNOCC)
               DO NA = 1, (NUNOCC)

                  IF (WDES%LORBITALREAL) THEN
                     E_CCSD = E_CCSD + OOVV_R(NI, NJ, NA, NB)* &
                              (T2_N_R(NA, NB, NI, NJ, KI, KJ, KA))*KPOINTS_FULL%WTKPT(1) &
                              *KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1)
                  ELSE
                     E_CCSD = E_CCSD + OOVV(NI, NJ, NA, NB)* &
                              (T2_N(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA))*WTKPT &
                              *WTKPT*WTKPT
                  END IF

               END DO
               END DO
            END DO
            END DO

         END DO
         END DO
      END DO
      CALL M_sum_z(WGW%COMM_INTER, E_CCSD, 1)
      CALL M_sum_z(WGW%COMM_INTER, E_CCSD_X, 1)
      E_CCSD = E_CCSD + ETMP*WTKPT
      E_CCSD_X = E_CCSD_X + ETMP*WTKPT
      IF (IO%IU0 >= 0) WRITE (IO%IU0, *) 'E_PPL=', E_CCSD
      IF (IO%IU0 >= 0) WRITE (IO%IU6, *) 'E PPL=', E_CCSD
      IF (IO%IU0 >= 0) WRITE (IO%IU0, *) 'E_PPL_X=', E_CCSD_X
      IF (IO%IU0 >= 0) WRITE (IO%IU6, *) 'E_PPL_X=', E_CCSD_X
      IF (ITERATION == 1) E_MP2 = E_CCSD
      IF ((ITERATION == MAX_IT) .and. (IO%IU6 >= 0)) WRITE (IO%IU6, *) 'final E_PPL=', E_CCSD

   END SUBROUTINE CALC_CCSD_ENERGY_PPL

!***********************************************************************
!
!  This routine broadcasts a set of hole-hole Coulomb-Vertices to all MPI processes
!
!***********************************************************************

   Subroutine BCAST2ALL_FTOD_IJ(WDES, KI, MNCC, PW_IJ_TMP, OC_IJ_TMP)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KI, MNCC
      GDEF :: PW_IJ_TMP(:, :, :, :), OC_IJ_TMP(:, :, :, :)

      call BLACS_PINFO(ME, PROCS)

      IF (PROCS_KPTS(KI) == ME) THEN
               CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, VBMAX*VBMAX*REALNKPTS, FTOD_PW_IJ(1,1,1,1,MKPTS_KPTS(KI),MIN(NCC,MNCC)), NGVECTOR)
         PW_IJ_TMP(:, :, :, :) = FTOD_PW_IJ(:, :, :, :, MKPTS_KPTS(KI), MIN(NCC, MNCC))
      ELSE
     CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, VBMAX*VBMAX*REALNKPTS, PW_IJ_TMP(1, 1, 1, 1), NGVECTOR, 0, PROCS_KPTS(KI))
      END IF

      IF (ASSOCIATED(H)) THEN

         IF (PROCS_KPTS(KI) == ME) THEN
 CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, VBMAX*VBMAX*REALNKPTS, FTOD_OC_IJ(1, 1, 1, 1, MKPTS_KPTS(KI), MNCC), NHVECTOR)
            OC_IJ_TMP(:, :, :, :) = FTOD_OC_IJ(:, :, :, :, MKPTS_KPTS(KI), MNCC)
         ELSE
     CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, VBMAX*VBMAX*REALNKPTS, OC_IJ_TMP(1, 1, 1, 1), NHVECTOR, 0, PROCS_KPTS(KI))
         END IF
      END IF
   END Subroutine BCAST2ALL_FTOD_IJ

!***********************************************************************
!
!  This routine broadcasts a set of particle-hole Coulomb-Vertices to all MPI processes
!
!***********************************************************************

   Subroutine BCAST2ALL_FTOD_IA(WDES, KI, MNCC, MPW_IA_TMP, MOC_IA_TMP)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KI, MNCC
      GDEF :: MPW_IA_TMP(:, :, :, :), MOC_IA_TMP(:, :, :, :)

      MPW_IA_TMP(:, :, :, :) = zero
      IF (ASSOCIATED(H)) MOC_IA_TMP(:, :, :, :) = zero
      call BLACS_PINFO(ME, PROCS)

      IF (PROCS_KPTS(KI) == ME) THEN
               CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, (NUNOCC)*VBMAX*REALNKPTS, FTOD_PW_IA(1,1,1,1,MKPTS_KPTS(KI),MIN(NCC,MNCC)), NGVECTOR)
         MPW_IA_TMP(:, :, :, :) = FTOD_PW_IA(:, :, :, :, MKPTS_KPTS(KI), MIN(NCC, MNCC))
      ELSE
 CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, (NUNOCC)*VBMAX*REALNKPTS, MPW_IA_TMP(1, 1, 1, 1), NGVECTOR, 0, PROCS_KPTS(KI))
      END IF

      IF (ASSOCIATED(H)) THEN

         IF (PROCS_KPTS(KI) == ME) THEN
   CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, (NUNOCC)*VBMAX*REALNKPTS, FTOD_OC_IA(1,1,1,1,MKPTS_KPTS(KI),MNCC), NHVECTOR)
            MOC_IA_TMP(:, :, :, :) = FTOD_OC_IA(:, :, :, :, MKPTS_KPTS(KI), MNCC)
         ELSE
 CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, (NUNOCC)*VBMAX*REALNKPTS, MOC_IA_TMP(1, 1, 1, 1), NHVECTOR, 0, PROCS_KPTS(KI))
         END IF

      END IF

   END Subroutine BCAST2ALL_FTOD_IA

!***********************************************************************
!
!  This routine broadcasts a set of particle-hole Coulomb-Vertices to all MPI processes
!
!***********************************************************************

   Subroutine BCAST2ALL_FTOD_AI(WDES, KA, MNCC, PW_AI_TMP, OC_AI_TMP)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KA, MNCC
      GDEF :: PW_AI_TMP(:, :, :, :), OC_AI_TMP(:, :, :, :)

      PW_AI_TMP(:, :, :, :) = zero
      IF (ASSOCIATED(H)) OC_AI_TMP(:, :, :, :) = zero
      call BLACS_PINFO(ME, PROCS)

      IF (PROCS_KPTS(KA) == ME) THEN
               CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, VBMAX*(NUNOCC)*REALNKPTS, FTOD_PW_AI(1,1,1,1,MKPTS_KPTS(KA),MIN(NCC,MNCC)), NGVECTOR)
         PW_AI_TMP(:, :, :, :) = FTOD_PW_AI(:, :, :, :, MKPTS_KPTS(KA), MIN(NCC, MNCC))
      ELSE
  CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, VBMAX*(NUNOCC)*REALNKPTS, PW_AI_TMP(1, 1, 1, 1), NGVECTOR, 0, PROCS_KPTS(KA))
      END IF

      IF (ASSOCIATED(H)) THEN

         IF (PROCS_KPTS(KA) == ME) THEN
   CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, VBMAX*(NUNOCC)*REALNKPTS, FTOD_OC_AI(1,1,1,1,MKPTS_KPTS(KA),MNCC), NHVECTOR)
            OC_AI_TMP(:, :, :, :) = FTOD_OC_AI(:, :, :, :, MKPTS_KPTS(KA), MNCC)
         ELSE
  CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, VBMAX*(NUNOCC)*REALNKPTS, OC_AI_TMP(1, 1, 1, 1), NHVECTOR, 0, PROCS_KPTS(KA))
         END IF

      END IF

   END Subroutine BCAST2ALL_FTOD_AI

!***********************************************************************
!
!  This routine broadcasts a set of particle-particle Coulomb-Vertices to all MPI processes
!
!***********************************************************************

   Subroutine BCAST2ALL_FTOD_AB(WDES, KA, MNCC, PW_AB_TMP, OC_AB_TMP)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KA, MNCC, MKA
      GDEF :: PW_AB_TMP(:, :, :, :), OC_AB_TMP(:, :, :, :)
      integer :: time_array1(8), time_array2(8), ems

      call date_and_time(values=time_array1)

      call BLACS_PINFO(ME, PROCS)
      DO MKA = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(MKA) /= 0)) CYCLE
         IF (PROCS_KPTS(KA) == ME) THEN
               CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, (NUNOCC)*(NUNOCC), FTOD_PW_AB(1,1,1,RKQofKQ(MKA),MKPTS_KPTS(KA),MIN(NCC,MNCC)), NGVECTOR)
            PW_AB_TMP(:, :, :, RKQofKQ(MKA)) = FTOD_PW_AB(:, :, :, RKQofKQ(MKA), MKPTS_KPTS(KA), MIN(NCC, MNCC))
         ELSE
   CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, (NUNOCC)*(NUNOCC), PW_AB_TMP(1,1,1,RKQofKQ(MKA)), NGVECTOR,0,PROCS_KPTS(KA))
         END IF

         IF (ASSOCIATED(H)) THEN

            IF (PROCS_KPTS(KA) == ME) THEN
                  CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, (NUNOCC)*(NUNOCC), FTOD_OC_AB(1,1,1,RKQofKQ(MKA),MKPTS_KPTS(KA),MNCC), NHVECTOR)
               OC_AB_TMP(:, :, :, RKQofKQ(MKA)) = FTOD_OC_AB(:, :, :, RKQofKQ(MKA), MKPTS_KPTS(KA), MNCC)
            ELSE
   CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, (NUNOCC)*(NUNOCC), OC_AB_TMP(1,1,1,RKQofKQ(MKA)), NHVECTOR,0,PROCS_KPTS(KA))
            END IF

         END IF
      END DO
      call date_and_time(values=time_array2)
      ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
             time_array2(8) - time_array1(8))
      TETS = TETS + ems

   END Subroutine BCAST2ALL_FTOD_AB

!***********************************************************************
!
!  This routine broadcasts the \chi_klij intermediate to all MPI processes
!
!***********************************************************************

   Subroutine BCAST2ALL_CHI_KLIJ(WDES, KK, CHI_TMP, CHI_TMP_R)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KK, MNCC, CHIL, CHIM
      GDEFS :: CHI_TMP(:, :, :, :, :, :)
      REAL(qs) :: CHI_TMP_R(:, :, :, :, :, :)

      call BLACS_PINFO(ME, PROCS)

      IF (LORBREAL) THEN
         IF (PROCS_KPTS(KK) == ME) THEN
            CALL SGEBS2D(CONTXT_COLS, 'All', 'i-ring', VBMAX*VBMAX, (VBMAX*VBMAX*REALNKPTS*REALNKPTS), &
                         CHI_KLIJ_R(1, 1, 1, 1, 1, 1, MKPTS_KPTS(KK)), VBMAX*VBMAX)
            CHI_TMP_R(:, :, :, :, :, :) = CHI_KLIJ_R(:, :, :, :, :, :, MKPTS_KPTS(KK))
         ELSE
               CALL SGEBR2D(CONTXT_COLS, 'All', 'i-ring', VBMAX*VBMAX, (VBMAX*VBMAX*REALNKPTS*REALNKPTS), CHI_TMP_R(1,1,1,1,1,1), VBMAX*VBMAX,0,PROCS_KPTS(KK))
         END IF
      ELSE
         IF (PROCS_KPTS(KK) == ME) THEN
            CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', VBMAX*VBMAX, (VBMAX*VBMAX*REALNKPTS*REALNKPTS), &
                         CHI_KLIJ(1, 1, 1, 1, 1, 1, MKPTS_KPTS(KK)), VBMAX*VBMAX)
            CHI_TMP(:, :, :, :, :, :) = CHI_KLIJ(:, :, :, :, :, :, MKPTS_KPTS(KK))
         ELSE
               CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', VBMAX*VBMAX, (VBMAX*VBMAX*REALNKPTS*REALNKPTS), CHI_TMP(1,1,1,1,1,1), VBMAX*VBMAX,0,PROCS_KPTS(KK))
         END IF
      END IF

   END Subroutine BCAST2ALL_CHI_KLIJ

!***********************************************************************
!
!  This routine broadcasts a T2 amplitude to all MPI processes
!
!***********************************************************************

   Subroutine BCAST2ALL_T2(WDES, KA, T2_MTMP, T2_MTMP_R)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KA, CHIL, CHIM
      GDEFS :: T2_MTMP(:, :, :, :, :, :)
      REAL(qs) :: T2_MTMP_R(:, :, :, :, :, :)
      integer :: time_array1(8), time_array2(8), ems

      call BLACS_PINFO(ME, PROCS)
      call date_and_time(values=time_array1)

      IF (LORBREAL) THEN
         IF (PROCS_KPTS(KA) == ME) THEN
            CALL SGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC), (VBMAX*VBMAX*REALNKPTS*REALNKPTS), &
                         T2_R(1, 1, 1, 1, 1, 1, MKPTS_KPTS(KA)), (NUNOCC*NUNOCC))
            T2_MTMP_R(:, :, :, :, :, :) = T2_R(:, :, :, :, :, :, MKPTS_KPTS(KA))
         ELSE
            CALL SGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC), (VBMAX*VBMAX*REALNKPTS*REALNKPTS), &
                         T2_MTMP_R(1, 1, 1, 1, 1, 1), (NUNOCC*NUNOCC), 0, PROCS_KPTS(KA))
         END IF
      ELSE
         IF (PROCS_KPTS(KA) == ME) THEN
            CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC), (VBMAX*VBMAX*REALNKPTS*REALNKPTS), &
                         T2(1, 1, 1, 1, 1, 1, MKPTS_KPTS(KA)), (NUNOCC*NUNOCC))
            T2_MTMP(:, :, :, :, :, :) = T2(:, :, :, :, :, :, MKPTS_KPTS(KA))
         ELSE
            CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC), (VBMAX*VBMAX*REALNKPTS*REALNKPTS), &
                         T2_MTMP(1, 1, 1, 1, 1, 1), (NUNOCC*NUNOCC), 0, PROCS_KPTS(KA))
         END IF
      END IF

      call date_and_time(values=time_array2)
      ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
             time_array2(8) - time_array1(8))
      TET2S = TET2S + ems

   END Subroutine BCAST2ALL_T2

   Subroutine BCAST2ALL_T2_N(WDES, KA, T2_MTMP, T2_MTMP_R)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KA, CHIL, CHIM
      GDEFS :: T2_MTMP(:, :, :, :, :, :)
      REAL(qs) :: T2_MTMP_R(:, :, :, :, :, :)
      integer :: time_array1(8), time_array2(8), ems

      call BLACS_PINFO(ME, PROCS)
      call date_and_time(values=time_array1)

      IF (LORBREAL) THEN
         IF (PROCS_KPTS(KA) == ME) THEN
            CALL SGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX*REALNKPTS*REALNKPTS), &
                         T2_N_R(1, 1, 1, 1, 1, 1, MKPTS_KPTS(KA)), (NUNOCC))
            T2_MTMP_R(:, :, :, :, :, :) = T2_N_R(:, :, :, :, :, :, MKPTS_KPTS(KA))
         ELSE
            CALL SGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX*REALNKPTS*REALNKPTS), &
                         T2_MTMP_R(1, 1, 1, 1, 1, 1), (NUNOCC), 0, PROCS_KPTS(KA))
         END IF
      ELSE
         IF (PROCS_KPTS(KA) == ME) THEN
            CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX*REALNKPTS*REALNKPTS), &
                         T2_N(1, 1, 1, 1, 1, 1, MKPTS_KPTS(KA)), (NUNOCC))
            T2_MTMP(:, :, :, :, :, :) = T2_N(:, :, :, :, :, :, MKPTS_KPTS(KA))
         ELSE
            CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX*REALNKPTS*REALNKPTS), &
                         T2_MTMP(1, 1, 1, 1, 1, 1), (NUNOCC), 0, PROCS_KPTS(KA))
         END IF
      END IF

      call date_and_time(values=time_array2)
      ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
             time_array2(8) - time_array1(8))
      TET2S = TET2S + ems

   END Subroutine BCAST2ALL_T2_N

   Subroutine BCAST2ALL_T2_KJ_KA(WDES, KJ, KA, T2_MTMP, T2_MTMP_R)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KA, CHIL, CHIM, KJ
      GDEFS :: T2_MTMP(:, :, :, :, :)
      REAL(qs) :: T2_MTMP_R(:, :, :, :, :)
      integer :: time_array1(8), time_array2(8), ems

      call BLACS_PINFO(ME, PROCS)
      call date_and_time(values=time_array1)
      IF (LORBREAL) THEN
         IF (PROCS_KPTS(KA) == ME) THEN
            CALL SGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX*REALNKPTS), &
                         T2_R(1, 1, 1, 1, 1, KJ, MKPTS_KPTS(KA)), (NUNOCC))
            T2_MTMP_R(:, :, :, :, :) = T2_R(:, :, :, :, :, KJ, MKPTS_KPTS(KA))
         ELSE
            CALL SGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX*REALNKPTS), &
                         T2_MTMP_R(1, 1, 1, 1, 1), (NUNOCC), 0, PROCS_KPTS(KA))
         END IF
      ELSE
         IF (PROCS_KPTS(KA) == ME) THEN
            CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX*REALNKPTS), &
                         T2(1, 1, 1, 1, 1, RKIofKI(KJ), MKPTS_KPTS(KA)), (NUNOCC))
            T2_MTMP(:, :, :, :, :) = T2(:, :, :, :, :, RKIofKI(KJ), MKPTS_KPTS(KA))
         ELSE
            CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX*REALNKPTS), &
                         T2_MTMP(1, 1, 1, 1, 1), (NUNOCC), 0, PROCS_KPTS(KA))
         END IF
      END IF
      call date_and_time(values=time_array2)
      ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
             time_array2(8) - time_array1(8))
      TET2S = TET2S + ems

   END Subroutine BCAST2ALL_T2_KJ_KA

   Subroutine BCAST2ALL_T2_KI_KJ_KA(WDES, KI, KJ, KA, T2_MTMP, T2_MTMP_R)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KA, CHIL, CHIM, KJ, KI
      GDEFS :: T2_MTMP(:, :, :, :)
      REAL(qs) :: T2_MTMP_R(:, :, :, :)
      integer :: time_array1(8), time_array2(8), ems

      call BLACS_PINFO(ME, PROCS)
      call date_and_time(values=time_array1)
      IF (LORBREAL) THEN
         IF (PROCS_KPTS(KA) == ME) THEN
            CALL SGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX), &
                         T2_R(1, 1, 1, 1, KI, KJ, MKPTS_KPTS(KA)), (NUNOCC))
            T2_MTMP_R(:, :, :, :) = T2_R(:, :, :, :, KI, KJ, MKPTS_KPTS(KA))
         ELSE
            CALL SGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX), &
                         T2_MTMP_R(1, 1, 1, 1), (NUNOCC), 0, PROCS_KPTS(KA))
         END IF
      ELSE
         IF (PROCS_KPTS(KA) == ME) THEN
            CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX), &
                         T2(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KJ), MKPTS_KPTS(KA)), (NUNOCC))
            T2_MTMP(:, :, :, :) = T2(:, :, :, :, RKIofKI(KI), RKIofKI(KJ), MKPTS_KPTS(KA))
         ELSE
            CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX), &
                         T2_MTMP(1, 1, 1, 1), (NUNOCC), 0, PROCS_KPTS(KA))
         END IF
      END IF
      call date_and_time(values=time_array2)
      ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
             time_array2(8) - time_array1(8))
      TET2S = TET2S + ems

   END Subroutine BCAST2ALL_T2_KI_KJ_KA

   Subroutine BCAST2ALL_AS_T2_KJ_KA(WDES, KJ, KA, T2_MTMP, T2_MTMP_R)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KA, CHIL, CHIM, KJ, KK, KB, NJ, NK
      GDEFS :: T2_MTMP(:, :, :, :, :)
      REAL(qs) :: T2_MTMP_R(:, :, :, :, :)
      integer :: time_array1(8), time_array2(8), ems

      call BLACS_PINFO(ME, PROCS)
      call date_and_time(values=time_array1)
      IF (LORBREAL) THEN
         IF (PROCS_KPTS(KA) == ME) THEN
            DO KK = 1, WDES%NKPTS
            DO NK = 1, VBMAX
            DO NJ = 1, VBMAX
                  IF ((.not. LDISTING) .and. ((.NOT. (RING)) .or. (ASYM_RING))) T2_MTMP_R(:,:,NK,NJ,KK)=-T2_R(:,:,NK,NJ,KK,KJ,MKPTS_KPTS(KA))+(2.0_qs)*T2_R(:,:,NJ,NK,KJ,KK,MKPTS_KPTS(KA))
               IF (RING .and. (.not. ASYM_RING)) T2_MTMP_R(:, :, NK, NJ, KK) = (2.0_qs)*T2_R(:, :, NJ, NK, KJ, KK, MKPTS_KPTS(KA))
         IF ((LDISTING)) T2_MTMP_R(:,:,NK,NJ,KK)=-T2_R(:,:,NK,NJ,KK,KJ,MKPTS_KPTS(KA))+(2.0_qs)*T2_R(:,:,NJ,NK,KJ,KK,MKPTS_KPTS(KA))
            END DO
            END DO
            END DO

            CALL SGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX*WDES%NKPTS), &
                         T2_MTMP_R(1, 1, 1, 1, 1), (NUNOCC))
         ELSE
            CALL SGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX*WDES%NKPTS), &
                         T2_MTMP_R(1, 1, 1, 1, 1), (NUNOCC), 0, PROCS_KPTS(KA))
         END IF
      ELSE
         IF (PROCS_KPTS(KA) == ME) THEN
            DO KK = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KK) == 0)) CYCLE
               DO NK = 1, VBMAX
               DO NJ = 1, VBMAX
                  IF ((.not. LDISTING) .and. ((.NOT. (RING)) .or. (ASYM_RING))) T2_MTMP(:,:,NK,NJ,RKIofKI(KK))=-T2(:,:,NK,NJ,RKIofKI(KK),RKIofKI(KJ),MKPTS_KPTS(KA))+(2.0_qs,0.0_qs)*T2(:,:,NJ,NK,RKIofKI(KJ),RKIofKI(KK),MKPTS_KPTS(KA))
                  IF (RING .and. (.not. ASYM_RING)) T2_MTMP(:,:,NK,NJ,RKIofKI(KK))=(2.0_qs,0.0_qs)*T2(:,:,NJ,NK,RKIofKI(KJ),RKIofKI(KK),MKPTS_KPTS(KA))
                  IF ((LDISTING)) T2_MTMP(:,:,NK,NJ,RKIofKI(KK))=-T2(:,:,NK,NJ,RKIofKI(KK),RKIofKI(KJ),MKPTS_KPTS(KA))+(2.0_qs,0.0_qs)*T2(:,:,NJ,NK,RKIofKI(KJ),RKIofKI(KK),MKPTS_KPTS(KA))
               END DO
               END DO
            END DO

            CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX*REALNKPTS), &
                         T2_MTMP(1, 1, 1, 1, 1), (NUNOCC))
         ELSE
            CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX*REALNKPTS), &
                         T2_MTMP(1, 1, 1, 1, 1), (NUNOCC), 0, PROCS_KPTS(KA))
         END IF
      END IF

      call date_and_time(values=time_array2)
      ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
             time_array2(8) - time_array1(8))
      TET2S = TET2S + ems

   END Subroutine BCAST2ALL_AS_T2_KJ_KA

   Subroutine BCAST2ALL_ASonly_T2_KJ_KA(WDES, KJ, KA, T2_MTMP, T2_MTMP_R)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KA, CHIL, CHIM, KJ, KK, KB, NJ, NK
      GDEFS :: T2_MTMP(:, :, :, :, :)
      REAL(qs) :: T2_MTMP_R(:, :, :, :, :)
      integer :: time_array1(8), time_array2(8), ems

      call BLACS_PINFO(ME, PROCS)
      call date_and_time(values=time_array1)
      IF (LORBREAL) THEN
         IF (PROCS_KPTS(KA) == ME) THEN
            DO KK = 1, WDES%NKPTS
            DO NK = 1, VBMAX
            DO NJ = 1, VBMAX
                  IF ((.not. LDISTING) .and. ((.NOT. (RING)) .or. (ASYM_RING))) T2_MTMP_R(:,:,NK,NJ,KK)=-T2_R(:,:,NK,NJ,KK,KJ,MKPTS_KPTS(KA))+(2.0_qs)*T2_R(:,:,NJ,NK,KJ,KK,MKPTS_KPTS(KA))
               IF (RING .and. (.not. ASYM_RING)) T2_MTMP_R(:, :, NK, NJ, KK) = (2.0_qs)*T2_R(:, :, NJ, NK, KJ, KK, MKPTS_KPTS(KA))
               IF (LDISTING) T2_MTMP_R(:, :, NK, NJ, KK) = -T2_R(:, :, NK, NJ, KK, KJ, MKPTS_KPTS(KA))
            END DO
            END DO
            END DO

            CALL SGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX*WDES%NKPTS), &
                         T2_MTMP_R(1, 1, 1, 1, 1), (NUNOCC))
         ELSE
            CALL SGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX*WDES%NKPTS), &
                         T2_MTMP_R(1, 1, 1, 1, 1), (NUNOCC), 0, PROCS_KPTS(KA))
         END IF
      ELSE
         IF (PROCS_KPTS(KA) == ME) THEN
            DO KK = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KK) == 0)) CYCLE
               DO NK = 1, VBMAX
               DO NJ = 1, VBMAX
                  IF ((.not. LDISTING) .and. ((.NOT. (RING)) .or. (ASYM_RING))) T2_MTMP(:,:,NK,NJ,RKIofKI(KK))=-T2(:,:,NK,NJ,RKIofKI(KK),RKIofKI(KJ),MKPTS_KPTS(KA))+(2.0_qs,0.0_qs)*T2(:,:,NJ,NK,RKIofKI(KJ),RKIofKI(KK),MKPTS_KPTS(KA))
                  IF (RING .and. (.not. ASYM_RING)) T2_MTMP(:,:,NK,NJ,RKIofKI(KK))=(2.0_qs,0.0_qs)*T2(:,:,NJ,NK,RKIofKI(KJ),RKIofKI(KK),MKPTS_KPTS(KA))
                  IF (LDISTING) T2_MTMP(:, :, NK, NJ, RKIofKI(KK)) = -T2(:, :, NK, NJ, RKIofKI(KK), RKIofKI(KJ), MKPTS_KPTS(KA))
               END DO
               END DO
            END DO

            CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX*REALNKPTS), &
                         T2_MTMP(1, 1, 1, 1, 1), (NUNOCC))
         ELSE
            CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX*REALNKPTS), &
                         T2_MTMP(1, 1, 1, 1, 1), (NUNOCC), 0, PROCS_KPTS(KA))
         END IF
      END IF

      call date_and_time(values=time_array2)
      ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
             time_array2(8) - time_array1(8))
      TET2S = TET2S + ems

   END Subroutine BCAST2ALL_ASonly_T2_KJ_KA

   Subroutine BCAST2ALL_CHI_CKAI_KJ_KK_KA(WDES, KJ, KK, KA, T2_MTMP, T2_MTMP_R)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KA, CHIL, CHIM, KK, KJ
      GDEFS :: T2_MTMP(:, :, :, :)
      REAL(qs) :: T2_MTMP_R(:, :, :, :)
      integer :: time_array1(8), time_array2(8), ems

      call BLACS_PINFO(ME, PROCS)
      call date_and_time(values=time_array1)
      IF (LORBREAL) THEN
         IF (PROCS_KPTS(KA) == ME) THEN
            CALL SGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX), &
                         CHI_CKAI_R(1, 1, 1, 1, KJ, KK, MKPTS_KPTS(KA)), (NUNOCC))
            T2_MTMP_R(:, :, :, :) = CHI_CKAI_R(:, :, :, :, KJ, KK, MKPTS_KPTS(KA))
         ELSE
            CALL SGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX), &
                         T2_MTMP_R(1, 1, 1, 1), (NUNOCC), 0, PROCS_KPTS(KA))
         END IF
      ELSE
         IF (PROCS_KPTS(KA) == ME) THEN
            CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX), &
                         CHI_CKAI(1, 1, 1, 1, RKIofKI(KJ), RKIofKI(KK), MKPTS_KPTS(KA)), (NUNOCC))
            T2_MTMP(:, :, :, :) = CHI_CKAI(:, :, :, :, RKIofKI(KJ), RKIofKI(KK), MKPTS_KPTS(KA))
         ELSE
            CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (NUNOCC)*(VBMAX*VBMAX), &
                         T2_MTMP(1, 1, 1, 1), (NUNOCC), 0, PROCS_KPTS(KA))
         END IF
      END IF

      call date_and_time(values=time_array2)
      ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
             time_array2(8) - time_array1(8))
      TET2S = TET2S + ems

   END Subroutine BCAST2ALL_CHI_CKAI_KJ_KK_KA

!***********************************************************************
!
!  This routine does most of the BLAS level 3 work needed to compute
!  all kinds of Coulomb integrals on-the-fly from the Coulomb-Vertices.
!
!***********************************************************************

   SUBROUTINE CONTR_FTOD(NGV, NHV, PW_XY, OC_XY, LNXY, PW_XY_, OC_XY_, &
                         LNXY_, XY2E4ORB, XY2E4ORB_R, beta)
      IMPLICIT NONE
      INTEGER :: NGV, NHV, LNXY, LNXY_
      GDEF :: XY2E4ORB, PW_XY, OC_XY, &
         PW_XY_, OC_XY_, beta
      REAL(q) :: XY2E4ORB_R

      IF (.not. LORBREAL) THEN
         CALL ZGEMM(trans, 'n', (LNXY), (LNXY_), &
                    (NGV), (1._q, 0._q), PW_XY, (NGV), &
                    PW_XY_, (NGV), &
                    beta, XY2E4ORB, (LNXY))

         IF (ASSOCIATED(H)) THEN
            CALL ZGEMM(trans, 'n', (LNXY), (LNXY_), &
                       (NHV), (1._q, 0._q), OC_XY, (NHV), &
                       OC_XY_, (NHV), &
                       (1._q, 0._q), XY2E4ORB, (LNXY))
         END IF

      ELSE

         CALL DGEMM(trans, 'n', (LNXY), (LNXY_), &
                    (NGV*2), (1._q), PW_XY, (NGV*2), &
                    PW_XY_, (NGV*2), &
                    REAL(beta, kind=q), XY2E4ORB_R, (LNXY))

         IF (ASSOCIATED(H)) THEN

            CALL DGEMM(trans, 'n', (LNXY), (LNXY_), &
                       (NHV*2), (1._q), OC_XY, (NHV*2), &
                       OC_XY_, (NHV*2), &
                       (1._q), XY2E4ORB_R, (LNXY))
         END IF

      END IF

   END SUBROUTINE

!***********************************************************************
!
!  A bunch of sorting routines that resort different arrays in
!  a representation of particle/hole states with different dimensions.
!  Resorting is necessary for efficient contraction afterwards.
!
!***********************************************************************

   SUBROUTINE SORT_O1V1O2V2_V2V1O2O1(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEF :: XY(:, :, :, :), XY2(:, :, :, :)
      REAL(q) :: XY_R(:, :, :, :), XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2

      IF (LORBREAL) THEN
         DO NV1 = 1, (NUNOCC)
         DO NO1 = 1, VBMAX
         DO NV2 = 1, (NUNOCC)
         DO NO2 = 1, VBMAX
            XY2_R(NV1, NV2, NO1, NO2) = XY_R(NO2, NV2, NO1, NV1)
         END DO
         END DO
         END DO
         END DO
      ELSE
         DO NV1 = 1, (NUNOCC)
         DO NO1 = 1, VBMAX
         DO NV2 = 1, (NUNOCC)
         DO NO2 = 1, VBMAX
            XY2(NV1, NV2, NO1, NO2) = XY(NO2, NV2, NO1, NV1)
         END DO
         END DO
         END DO
         END DO
      END IF

   END SUBROUTINE SORT_O1V1O2V2_V2V1O2O1

   SUBROUTINE SORT_V1V2O1O2_O1O2V1V2(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEF :: XY(:, :, :, :), XY2(:, :, :, :)
      REAL(q) :: XY_R(:, :, :, :), XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2

      IF (LORBREAL) THEN
         DO NO2 = 1, VBMAX
         DO NO1 = 1, VBMAX
         DO NV2 = 1, (NUNOCC)
         DO NV1 = 1, (NUNOCC)
            XY2_R(NO1, NO2, NV1, NV2) = XY_R(NV1, NV2, NO1, NO2)
         END DO
         END DO
         END DO
         END DO
      ELSE
         DO NO2 = 1, VBMAX
         DO NO1 = 1, VBMAX
         DO NV2 = 1, (NUNOCC)
         DO NV1 = 1, (NUNOCC)
            XY2(NO1, NO2, NV1, NV2) = XY(NV1, NV2, NO1, NO2)
         END DO
         END DO
         END DO
         END DO
      END IF

   END SUBROUTINE SORT_V1V2O1O2_O1O2V1V2

   SUBROUTINE SORT_V1V2O1O2_O1O2V2V1(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY(:, :, :, :), XY2(:, :, :, :)
      REAL(qs) :: XY_R(:, :, :, :), XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2

      IF (LORBREAL) THEN
      DO NO2 = 1, VBMAX
      DO NO1 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2_R(NO1, NO2, NV2, NV1) = XY_R(NV1, NV2, NO1, NO2)
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO NO2 = 1, VBMAX
      DO NO1 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2(NO1, NO2, NV2, NV1) = XY(NV1, NV2, NO1, NO2)
      END DO
      END DO
      END DO
      END DO
      END IF
   END SUBROUTINE SORT_V1V2O1O2_O1O2V2V1

   SUBROUTINE SORT_V1V2O1O2_O2O1V2V1(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY(:, :, :, :), XY2(:, :, :, :)
      REAL(q) :: XY_R(:, :, :, :), XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2

      IF (LORBREAL) THEN
      DO NO2 = 1, VBMAX
      DO NO1 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2_R(NO2, NO1, NV2, NV1) = XY_R(NV1, NV2, NO1, NO2)
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO NO2 = 1, VBMAX
      DO NO1 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2(NO2, NO1, NV2, NV1) = XY(NV1, NV2, NO1, NO2)
      END DO
      END DO
      END DO
      END DO
      END IF
   END SUBROUTINE SORT_V1V2O1O2_O2O1V2V1

   SUBROUTINE SORT_O1O2V1V2_V2O2O1V1(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY(:, :, :, :), XY2(:, :, :, :)
      REAL(qs) :: XY_R(:, :, :, :), XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2

      IF (LORBREAL) THEN
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
      DO NO2 = 1, VBMAX
      DO NO1 = 1, VBMAX
         XY2_R(NV2, NO2, NO1, NV1) = XY_R(NO1, NO2, NV1, NV2)
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
      DO NO2 = 1, VBMAX
      DO NO1 = 1, VBMAX
         XY2(NV2, NO2, NO1, NV1) = XY(NO1, NO2, NV1, NV2)
      END DO
      END DO
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_O1O2V1V2_V2O2O1V1

   SUBROUTINE SORT_O1O2V1V2_V1O1V2O2(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY(:, :, :, :), XY2(:, :, :, :)
      REAL(qs) :: XY_R(:, :, :, :), XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2

      IF (LORBREAL) THEN
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
      DO NO2 = 1, VBMAX
      DO NO1 = 1, VBMAX
         XY2_R(NV1, NO1, NV2, NO2) = XY_R(NO1, NO2, NV1, NV2)
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
      DO NO2 = 1, VBMAX
      DO NO1 = 1, VBMAX
         XY2(NV1, NO1, NV2, NO2) = XY(NO1, NO2, NV1, NV2)
      END DO
      END DO
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_O1O2V1V2_V1O1V2O2

   SUBROUTINE SORT_V1V2_V2V1(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY(:, :), XY2(:, :)
      REAL(qs) :: XY_R(:, :), XY2_R(:, :)
      INTEGER :: NV1, NV2

      IF (LORBREAL) THEN
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2_R(NV2, NV1) = XY_R(NV1, NV2)
      END DO
      END DO
      ELSE
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2(NV2, NV1) = XY(NV1, NV2)
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_V1V2_V2V1

   SUBROUTINE SORT_V1O1_O1V1(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY(:, :), XY2(:, :)
      REAL(qs) :: XY_R(:, :), XY2_R(:, :)
      INTEGER :: NV1, NO1

      IF (LORBREAL) THEN
      DO NO1 = 1, (VBMAX)
      DO NV1 = 1, (NUNOCC)
         XY2_R(NO1, NV1) = XY_R(NV1, NO1)
      END DO
      END DO
      ELSE
      DO NO1 = 1, (VBMAX)
      DO NV1 = 1, (NUNOCC)
         XY2(NO1, NV1) = XY(NV1, NO1)
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_V1O1_O1V1

   SUBROUTINE SORT_O1V1_V1O1(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY(:, :), XY2(:, :)
      REAL(qs) :: XY_R(:, :), XY2_R(:, :)
      INTEGER :: NV1, NO1

      IF (LORBREAL) THEN
      DO NO1 = 1, (VBMAX)
      DO NV1 = 1, (NUNOCC)
         XY2_R(NV1, NO1) = XY_R(NO1, NV1)
      END DO
      END DO
      ELSE
      DO NO1 = 1, (VBMAX)
      DO NV1 = 1, (NUNOCC)
         XY2(NV1, NO1) = XY(NO1, NV1)
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_O1V1_V1O1

   SUBROUTINE SORT_O1V1V2O2_V1V2O1O2(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEF :: XY(:, :, :, :)
      REAL(q) :: XY_R(:, :, :, :)
      GDEFS :: XY2(:, :, :, :)
      REAL(qs) :: XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2

      IF (LORBREAL) THEN
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
      DO NO1 = 1, VBMAX
         XY2_R(NV1, NV2, NO1, NO2) = XY_R(NO1, NV1, NV2, NO2)
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
      DO NO1 = 1, VBMAX
         XY2(NV1, NV2, NO1, NO2) = XY(NO1, NV1, NV2, NO2)
      END DO
      END DO
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_O1V1V2O2_V1V2O1O2

   SUBROUTINE SORT_V2V1O1O2_V1V2O1O2(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEF :: XY(:, :, :, :)
      REAL(q) :: XY_R(:, :, :, :)
      GDEFS :: XY2(:, :, :, :)
      REAL(qs) :: XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2

      IF (LORBREAL) THEN
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV1 = 1, (NUNOCC)
      DO NV2 = 1, (NUNOCC)
         XY2_R(NV1, NV2, NO1, NO2) = XY_R(NV2, NV1, NO1, NO2)
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV1 = 1, (NUNOCC)
      DO NV2 = 1, (NUNOCC)
         XY2(NV1, NV2, NO1, NO2) = XY(NV2, NV1, NO1, NO2)
      END DO
      END DO
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_V2V1O1O2_V1V2O1O2

   SUBROUTINE SORT_V2V1O1O2_V2O2V1O1(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEF :: XY(:, :, :, :)
      REAL(q) :: XY_R(:, :, :, :)
      GDEFS :: XY2(:, :, :, :)
      REAL(qs) :: XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2

      IF (LORBREAL) THEN
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV1 = 1, (NUNOCC)
      DO NV2 = 1, (NUNOCC)
         XY2_R(NV2, NO2, NV1, NO1) = XY_R(NV2, NV1, NO1, NO2)
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV1 = 1, (NUNOCC)
      DO NV2 = 1, (NUNOCC)
         XY2(NV2, NO2, NV1, NO1) = XY(NV2, NV1, NO1, NO2)
      END DO
      END DO
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_V2V1O1O2_V2O2V1O1

   SUBROUTINE SORT_V1V2O1O2_V2O1V1O2(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY(:, :, :, :)
      REAL(qs) :: XY_R(:, :, :, :)
      GDEFS :: XY2(:, :, :, :)
      REAL(qs) :: XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2

      IF (LORBREAL) THEN
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2_R(NV2, NO1, NV1, NO2) = XY_R(NV1, NV2, NO1, NO2)
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2(NV2, NO1, NV1, NO2) = XY(NV1, NV2, NO1, NO2)
      END DO
      END DO
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_V1V2O1O2_V2O1V1O2

   SUBROUTINE SORT_V1V2O1O2_V2O1V1O2_KTMP(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY(:, :, :, :, :)
      REAL(qs) :: XY_R(:, :, :, :, :)
      GDEFS :: XY2(:, :, :, :, :)
      REAL(qs) :: XY2_R(:, :, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2, KK

      IF (LORBREAL) THEN
      DO KK = 1, WDES%NKPTS
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2_R(NV2, NO1, NV1, NO2, KK) = XY_R(NV1, NV2, NO1, NO2, KK)
      END DO
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO KK = 1, REALNKPTS
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2(NV2, NO1, NV1, NO2, KK) = XY(NV1, NV2, NO1, NO2, KK)
      END DO
      END DO
      END DO
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_V1V2O1O2_V2O1V1O2_KTMP

   SUBROUTINE SORT_V1V2O1O2_V2V1O2O1(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY(:, :, :, :)
      REAL(qs) :: XY_R(:, :, :, :)
      GDEFS :: XY2(:, :, :, :)
      REAL(qs) :: XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2, KK

      IF (LORBREAL) THEN
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2_R(NV2, NV1, NO2, NO1) = XY_R(NV1, NV2, NO1, NO2)
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2(NV2, NV1, NO2, NO1) = XY(NV1, NV2, NO1, NO2)
      END DO
      END DO
      END DO
      END DO
      END IF
   END SUBROUTINE SORT_V1V2O1O2_V2V1O2O1

   SUBROUTINE SORT_T2_V1V2O1O2_V2O1V1O2(WDES, KK, KI, KA, XY2, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY2(:, :, :, :)
      REAL(qs) :: XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2, KK, KI, KA

      IF (LORBREAL) THEN
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2_R(NV2, NO1, NV1, NO2) = T2_R(NV1, NV2, NO1, NO2, KK, KI, KA)
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2(NV2, NO1, NV1, NO2) = T2(NV1, NV2, NO1, NO2, KK, KI, KA)
      END DO
      END DO
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_T2_V1V2O1O2_V2O1V1O2

   SUBROUTINE SORT_AS_T2_V1V2O1O2_V2O1V1O2(WDES, KK, KI, KA, XY2, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY2(:, :, :, :)
      REAL(qs) :: XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2, KK, KI, KA
      IF (LORBREAL) THEN
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2_R(NV2, NO1, NV1, NO2) = 2.0_qs*T2_R(NV1, NV2, NO2, NO1, KI, KK, KA) - &
                                     T2_R(NV1, NV2, NO1, NO2, KK, KI, KA)
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2(NV2, NO1, NV1, NO2) = 2.0_qs*T2(NV1, NV2, NO2, NO1, RKIofKI(KI), RKIofKI(KK), RKIofKI(KA)) - &
                                   T2(NV1, NV2, NO1, NO2, RKIofKI(KK), RKIofKI(KI), RKIofKI(KA))
      END DO
      END DO
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_AS_T2_V1V2O1O2_V2O1V1O2

   SUBROUTINE SORT_V1O1V2O2_V1V2O1O2(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY(:, :, :, :)
      REAL(qs) :: XY_R(:, :, :, :)
      GDEFS :: XY2(:, :, :, :)
      REAL(qs) :: XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2

      IF (LORBREAL) THEN
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2_R(NV1, NV2, NO1, NO2) = XY_R(NV1, NO1, NV2, NO2)
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2(NV1, NV2, NO1, NO2) = XY(NV1, NO1, NV2, NO2)
      END DO
      END DO
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_V1O1V2O2_V1V2O1O2

   SUBROUTINE SORT_V1O1V2O2_V2O1V1O2(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY(:, :, :, :)
      REAL(qs) :: XY_R(:, :, :, :)
      GDEFS :: XY2(:, :, :, :)
      REAL(qs) :: XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2
      IF (LORBREAL) THEN
      DO NV1 = 1, (NUNOCC)
      DO NV2 = 1, (NUNOCC)
         XY2_R(NV1, :, NV2, :) = XY_R(NV2, :, NV1, :)
      END DO
      END DO
      ELSE
      DO NV1 = 1, (NUNOCC)
      DO NV2 = 1, (NUNOCC)
         XY2(NV1, :, NV2, :) = XY(NV2, :, NV1, :)
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_V1O1V2O2_V2O1V1O2

   SUBROUTINE SORT_V1V2O1O2_V2O2V1O1(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY(:, :, :, :)
      REAL(qs) :: XY_R(:, :, :, :)
      GDEFS :: XY2(:, :, :, :)
      REAL(qs) :: XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2
      integer :: time_array1(8), time_array2(8), ems

      call date_and_time(values=time_array1)
      IF (LORBREAL) THEN
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2_R(NV2, NO2, NV1, NO1) = XY_R(NV1, NV2, NO1, NO2)
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2(NV2, NO2, NV1, NO1) = XY(NV1, NV2, NO1, NO2)
      END DO
      END DO
      END DO
      END DO

      END IF
      call date_and_time(values=time_array2)
      ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
             time_array2(8) - time_array1(8))
      TETSORT = TETSORT + ems

   END SUBROUTINE SORT_V1V2O1O2_V2O2V1O1

   SUBROUTINE SORT_V1V2O1O2_V1O1V2O2(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY(:, :, :, :)
      REAL(qs) :: XY_R(:, :, :, :)
      GDEFS :: XY2(:, :, :, :)
      REAL(qs) :: XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2
      integer :: time_array1(8), time_array2(8), ems

      call date_and_time(values=time_array1)
      IF (LORBREAL) THEN
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
         XY2_R(:, NO1, NV2, NO2) = XY_R(:, NV2, NO1, NO2)
      END DO
      END DO
      END DO
      ELSE
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
         XY2(:, NO1, NV2, NO2) = XY(:, NV2, NO1, NO2)
      END DO
      END DO
      END DO
      END IF
      call date_and_time(values=time_array2)
      ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
             time_array2(8) - time_array1(8))
      TETSORT = TETSORT + ems

   END SUBROUTINE SORT_V1V2O1O2_V1O1V2O2

   SUBROUTINE SORT_V1V2O1O2_V1V2O2O1(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY(:, :, :, :)
      REAL(qs) :: XY_R(:, :, :, :)
      GDEFS :: XY2(:, :, :, :)
      REAL(qs) :: XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2
      integer :: time_array1(8), time_array2(8), ems
      IF (LORBREAL) THEN
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
         XY2_R(:, :, NO2, NO1) = XY_R(:, :, NO1, NO2)
      END DO
      END DO
      ELSE
      DO NO1 = 1, VBMAX
      DO NO2 = 1, VBMAX
         XY2(:, :, NO2, NO1) = XY(:, :, NO1, NO2)
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_V1V2O1O2_V1V2O2O1

   SUBROUTINE SORT_O1V1V2O2_V2O2V1O1(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEF :: XY(:, :, :, :)
      REAL(q) :: XY_R(:, :, :, :)
      GDEFS :: XY2(:, :, :, :)
      REAL(qs) :: XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2

      IF (LORBREAL) THEN
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
      DO NO1 = 1, VBMAX
         XY2_R(NV2, NO2, NV1, NO1) = XY_R(NO1, NV1, NV2, NO2)
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
      DO NO1 = 1, VBMAX
         XY2(NV2, NO2, NV1, NO1) = XY(NO1, NV1, NV2, NO2)
      END DO
      END DO
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_O1V1V2O2_V2O2V1O1

   SUBROUTINE SORT_O2V1V2O1_V1O1V2O2(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY(:, :, :, :)
      REAL(qs) :: XY_R(:, :, :, :)
      GDEFS :: XY2(:, :, :, :)
      REAL(qs) :: XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2

      IF (LORBREAL) THEN
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
      DO NO1 = 1, VBMAX
         XY2_R(NV1, NO1, NV2, NO2) = XY_R(NO2, NV1, NV2, NO1)
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
      DO NO1 = 1, VBMAX
         XY2(NV1, NO1, NV2, NO2) = XY(NO2, NV1, NV2, NO1)
      END DO
      END DO
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_O2V1V2O1_V1O1V2O2

   SUBROUTINE SORT_V2O1V1O2_V1V2O1O2(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY(:, :, :, :)
      REAL(qs) :: XY_R(:, :, :, :)
      GDEFS :: XY2(:, :, :, :)
      REAL(qs) :: XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2

      IF (LORBREAL) THEN
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
      DO NO1 = 1, VBMAX
         XY2_R(NV1, NV2, NO1, NO2) = XY_R(NV2, NO1, NV1, NO2)
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO NO2 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
      DO NO1 = 1, VBMAX
         XY2(NV1, NV2, NO1, NO2) = XY(NV2, NO1, NV1, NO2)
      END DO
      END DO
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_V2O1V1O2_V1V2O1O2

   SUBROUTINE SORT_V1V2O1O2_O1V1V2O2(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEFS :: XY(:, :, :, :)
      REAL(qs) :: XY_R(:, :, :, :)
      GDEFS :: XY2(:, :, :, :)
      REAL(qs) :: XY2_R(:, :, :, :)
      INTEGER :: NV1, NO1, NV2, NO2

      IF (LORBREAL) THEN
      DO NO2 = 1, VBMAX
      DO NO1 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2_R(NO1, NV1, NV2, NO2) = XY_R(NV1, NV2, NO1, NO2)
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO NO2 = 1, VBMAX
      DO NO1 = 1, VBMAX
      DO NV2 = 1, (NUNOCC)
      DO NV1 = 1, (NUNOCC)
         XY2(NO1, NV1, NV2, NO2) = XY(NV1, NV2, NO1, NO2)
      END DO
      END DO
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_V1V2O1O2_O1V1V2O2

   SUBROUTINE SORT_JLIK_KLIJ(WDES, XY, XY2, XY_R, XY2_R)
      IMPLICIT NONE
      TYPE(wavedes) WDES
      GDEF :: XY(:, :, :, :), XY2(:, :, :, :)
      REAL(q) :: XY_R(:, :, :, :)
      REAL(q) :: XY2_R(:, :, :, :)
      INTEGER :: NK, NL, NI, NJ
      IF (LORBREAL) THEN
      DO NK = 1, VBMAX
      DO NL = 1, VBMAX
      DO NI = 1, VBMAX
      DO NJ = 1, VBMAX
         XY2_R(NK, NL, NI, NJ) = XY_R(NJ, NL, NI, NK)
      END DO
      END DO
      END DO
      END DO
      ELSE
      DO NK = 1, VBMAX
      DO NL = 1, VBMAX
      DO NI = 1, VBMAX
      DO NJ = 1, VBMAX
         XY2(NK, NL, NI, NJ) = XY(NJ, NL, NI, NK)
      END DO
      END DO
      END DO
      END DO
      END IF

   END SUBROUTINE SORT_JLIK_KLIJ


!***********************************************************************
!
!  This routine updates the T2 amplitudes with a new guess depending on the
!  employed solver and solver parameters.
!
!***********************************************************************

   SUBROUTINE UPDATE_T2(WDES, WGW, W, IO, ITERATION)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      TYPE(wavespin) W
      TYPE(in_struct) IO
      INTEGER :: KQ, KI, KJ, KQ_, NI, NJ, KA, KB, NA, NB, ITERATION
      COMPLEX(q) :: E_CCSD, TMP, TMPPRET2
      COMPLEX(qs) :: GRADIENT
      INTEGER :: NTMP
      REAL(8) :: RMSD

      RES = zero
      NTMP = 0

      DO KB = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0.0_q)) CYCLE
         DO KA = 1, MY_NKPTS
         DO KQ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
            KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) + &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
            KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)

            DO NI = 1, VBMAX
            DO NJ = 1, VBMAX
               DO NB = 1, (NUNOCC)
               DO NA = 1, (NUNOCC)
                  IF (ITERATION == 1) THEN
                     IF (LORBREAL) THEN
                        T2_R(NA, NB, NI, NJ, KI, KJ, KA) = (T2_N_R(NA, NB, NI, NJ, KI, KJ, KA)) &
                                                           /(W%CELTOT(NI_N(NI), KI, 1) + &
                                                             W%CELTOT(NI_N(NJ), KJ, 1) - &
                                                             W%CELTOT(NA_N(NA), KPTS_MKPTS(KA), 1) - &
                                                             W%CELTOT(NA_N(NB), KB, 1))
                     ELSE
                        NTMP = NTMP + 1
                        T2(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA) = (T2_N(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA)) &
                                                                           /(W%CELTOT(NI_N(NI), KI, 1) + &
                                                                             W%CELTOT(NI_N(NJ), KJ, 1) - &
                                                                             W%CELTOT(NA_N(NA), KPTS_MKPTS(KA), 1) - &
                                                                             W%CELTOT(NA_N(NB), KB, 1))
                     END IF
                     CYCLE
                  ELSE
                     IF (LORBREAL) THEN
                        T2_R(NA, NB, NI, NJ, KI, KJ, KA) = T2_R(NA, NB, NI, NJ, KI, KJ, KA)*((1.0_qs) - REAL(mix, kind=qs)) + &
                                 (CDIAG*T2_R(NA, NB, NI, NJ, KI, KJ, KA) + T2_N_R(NA, NB, NI, NJ, KI, KJ, KA))*REAL(mix, kind=qs)/ &
                                                           (W%CELTOT(NI_N(NI), KI, 1) + W%CELTOT(NI_N(NJ), KJ, 1) - &
                                                          W%CELTOT(NA_N(NA), KPTS_MKPTS(KA), 1) - W%CELTOT(NA_N(NB), KB, 1) + CDIAG)
                     ELSE
                        NTMP = NTMP + 1
    T2(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA) = T2(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA)*((1.0_qs, 0.0_qs) - mix) + &
                (CDIAG*T2(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA) + T2_N(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA))*mix/ &
                                                                          (W%CELTOT(NI_N(NI), KI, 1) + W%CELTOT(NI_N(NJ), KJ, 1) - &
                                                          W%CELTOT(NA_N(NA), KPTS_MKPTS(KA), 1) - W%CELTOT(NA_N(NB), KB, 1) + CDIAG)
                     END IF

                  END IF

               END DO
               END DO
            END DO
            END DO
         END DO
         END DO
      END DO

      IF (LMETAL) CALL NULLIFY_T2(WDES, WGW, W, IO, ITERATION)

      IF (PRECONDITION) PRET2 = (0.0_qs, 0.0_qs)
      CALLMPI(M_sum_z(WGW%COMM_INTER, RES, 1))


   END SUBROUTINE UPDATE_T2


!***********************************************************************
!
!  This routine updates the T2 amplitudes for the calculation of the
!  particle-particle ladder contribution to the CCSD energy
!
!***********************************************************************

   SUBROUTINE UPDATE_T2_LADDER(WDES, WGW, W, IO, ITERATION)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      TYPE(wavespin) W
      TYPE(in_struct) IO
      INTEGER :: KQ, KI, KJ, KQ_, NI, NJ, KA, KB, NA, NB, ITERATION
      COMPLEX(q) :: E_CCSD, TMP, TMPPRET2
      COMPLEX(qs) :: GRADIENT
      INTEGER :: NTMP
      REAL(8) :: RMSD

      RES = zero
      NTMP = 0

      DO KB = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0.0_q)) CYCLE
         DO KA = 1, MY_NKPTS
         DO KQ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
            KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) + &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
            KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)

            DO NI = 1, VBMAX
            DO NJ = 1, VBMAX
               DO NB = 1, (NUNOCC)
               DO NA = 1, (NUNOCC)
                  IF (LORBREAL) THEN
                     T2_N_R(NA, NB, NI, NJ, KI, KJ, KA) = (T2_N_R(NA, NB, NI, NJ, KI, KJ, KA)) &
                                                          /(W%CELTOT(NI_N(NI), KI, 1) + &
                                                            W%CELTOT(NI_N(NJ), KJ, 1) - &
                                                            W%CELTOT(NA_N(NA), KPTS_MKPTS(KA), 1) - &
                                                            W%CELTOT(NA_N(NB), KB, 1))
                  ELSE
                     NTMP = NTMP + 1
                     T2_N(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA) = (T2_N(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA)) &
                                                                          /(W%CELTOT(NI_N(NI), KI, 1) + &
                                                                            W%CELTOT(NI_N(NJ), KJ, 1) - &
                                                                            W%CELTOT(NA_N(NA), KPTS_MKPTS(KA), 1) - &
                                                                            W%CELTOT(NA_N(NB), KB, 1))
                  END IF
                  CYCLE

               END DO
               END DO
            END DO
            END DO
         END DO
         END DO
      END DO

      IF (LMETAL) CALL NULLIFY_T2(WDES, WGW, W, IO, ITERATION)

   END SUBROUTINE UPDATE_T2_LADDER

!***********************************************************************
!
!  This routine updates the T1 amplitudes with a new guess depending on the
!  employed solver and solver parameters.
!
!***********************************************************************

   SUBROUTINE UPDATE_T1(WDES, WGW, W, IO, iteration)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      TYPE(wavespin) W
      TYPE(in_struct) IO
      INTEGER :: iteration
      INTEGER :: KQ, KI, KJ, KQ_, NI, NJ, KA, KB, NA, NB
      COMPLEX(q) :: E_CCSD
      INTEGER :: NTMP
      REAL(8) :: RMSD

      IF (LORBREAL) THEN
      DO KI = 1, WDES%NKPTS
         DO NI = 1, VBMAX
            DO NA = 1, (NUNOCC)
               T1_R(NA, NI, KI) = (T1_N_R(NA, NI, KI)*REAL(mix, kind=qs)) &
                                  /(W%CELTOT(NI_N(NI), KI, 1) - &
                                    W%CELTOT(NA_N(NA), KI, 1)) + T1_R(NA, NI, KI)*((1.0_qs) - REAL(mix, kind=qs))
            END DO
         END DO
      END DO
      ELSE
      DO KI = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
         DO NI = 1, VBMAX
            DO NA = 1, (NUNOCC)
               T1(NA, NI, RKIofKI(KI)) = (CDIAG*T1(NA, NI, RKIofKI(KI)) + T1_N(NA, NI, RKIofKI(KI))*mix) &
                                         /(W%CELTOT(NI_N(NI), KI, 1) - &
                                           W%CELTOT(NA_N(NA), KI, 1) + CDIAG) + T1(NA, NI, RKIofKI(KI))*((1.0_qs, 0.0_qs) - mix)

            END DO
         END DO
      END DO
      END IF

      IF (LMETAL) CALL NULLIFY_T1(WDES, WGW, W, IO, ITERATION)

   END SUBROUTINE UPDATE_T1



!***********************************************************************
!
! This routine sets all T2 amplitude components to zero that do not represent
! a particle-hole pair excitation amplitude.
! N.B.: For simplicity we also compute T_ij^kl when LMETAL=.TRUE. ,
!       where k and l can be occupied orbitals. Obviously these contributions
!       need to be zero, which is guaranteed by this routine.
!
!***********************************************************************

   SUBROUTINE NULLIFY_T2(WDES, WGW, W, IO, ITERATION)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      TYPE(wavespin) W
      TYPE(in_struct) IO
      INTEGER :: KQ, KI, KJ, KQ_, NI, NJ, KA, KB, NA, NB, ITERATION
      REAL(q) :: WEIGHT, SMALL

      SMALL = 0.1_q

      DO KB = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
         DO KA = 1, MY_NKPTS
         DO KQ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
            KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) + &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
            KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)

            DO NI = 1, VBMAX
            DO NJ = 1, VBMAX
               DO NB = 1, (NUNOCC)
               DO NA = 1, (NUNOCC)

                  WEIGHT = W%FERTOT(NI, KI, 1)*W%FERTOT(NJ, KJ, 1)* &
                           (1.0_q - W%FERTOT(NA, KPTS_MKPTS(KA), 1))*(1.0_q - W%FERTOT(NB, KB, 1))

                  IF (LORBREAL) THEN
                     IF (WEIGHT < SMALL) T2_R(NA, NB, NI, NJ, KI, KJ, KA) = 0.0_qs
                     IF (WEIGHT < SMALL) T2_N_R(NA, NB, NI, NJ, KI, KJ, KA) = 0.0_qs
                  ELSE
                     IF (WEIGHT < SMALL) T2(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA) = (0.0_qs, 0.0_qs)
                     IF (WEIGHT < SMALL) T2_N(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA) = (0.0_qs, 0.0_qs)
                     IF (ISNANC(T2(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA))) THEN
                        WRITE (*, *) 'isnan', NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA, WEIGHT
                     END IF
                  END IF
               END DO
               END DO
            END DO
            END DO
         END DO
         END DO
      END DO

   END SUBROUTINE NULLIFY_T2

!***********************************************************************
!
! This routine sets all T1 amplitude components to zero that do not represent
! a particle-hole excitation amplitude.
! N.B.: For simplicity we also compute T_i^k when LMETAL=.TRUE. ,
!       where k can be an occupied orbital. Obviously these contributions
!       need to be zero, which is guaranteed by this routine.
!
!***********************************************************************

   SUBROUTINE NULLIFY_T1(WDES, WGW, W, IO, ITERATION)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      TYPE(wavespin) W
      TYPE(in_struct) IO
      INTEGER :: KQ, KI, KJ, KQ_, NI, NJ, KA, KB, NA, NB, ITERATION
      REAL(q) :: WEIGHT, SMALL

      SMALL = 0.1_q

      IF (LORBREAL) THEN
      DO KI = 1, WDES%NKPTS
         DO NI = 1, VBMAX
            DO NA = 1, (NUNOCC)
               WEIGHT = W%FERTOT(NI, KI, 1)*(1.0_q - W%FERTOT(NA, KI, 1))
               IF (WEIGHT < SMALL) T1_R(NA, NI, KI) = 0.0_qs
               IF (WEIGHT < SMALL) T1_N_R(NA, NI, KI) = 0.0_qs
            END DO
         END DO
      END DO
      ELSE
      DO KI = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
         DO NI = 1, VBMAX
            DO NA = 1, (NUNOCC)
               WEIGHT = W%FERTOT(NI, KI, 1)*(1.0_q - W%FERTOT(NA, KI, 1))
               IF (WEIGHT < SMALL) T1(NA, NI, RKIofKI(KI)) = (0.0_qs, 0.0_qs)
               IF (WEIGHT < SMALL) T1_N(NA, NI, RKIofKI(KI)) = (0.0_qs, 0.0_qs)
    IF (ISNANC(T1(NA, NI, RKIofKI(KI)))) WRITE (*, *) 'isnan', NA, NI, RKIofKI(KI), WEIGHT, W%FERTOT(NI, KI, 1), W%FERTOT(NA, KI, 1)
            END DO
         END DO
      END DO
      END IF

   END SUBROUTINE NULLIFY_T1

!***********************************************************************
!
! Very useful routine for debuggin that check is something is nan
!
!***********************************************************************

   FUNCTION ISNANC(X)
      IMPLICIT NONE
      LOGICAL :: ISNANC
      COMPLEX(qs) :: X

!!      IF ((ISNAN(REAL(X))) .OR. (ISNAN(AIMAG(X)))) THEN
      IF ((REAL(X)/=REAL(X)).OR.(AIMAG(X)/=AIMAG(X))) THEN
         ISNANC = .TRUE.
      ELSE
         ISNANC = .FALSE.
      END IF

   END FUNCTION ISNANC

!***********************************************************************
!
! The following routine computes an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE FORM_K_KI(WGW, WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      INTEGER :: KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
                 KC, KB_, KD, NC, ND
      COMPLEX(qs) :: KW2

      KW2 = KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1)*(2.0_qs, 0.0_qs)
      IF (LORBREAL) THEN
         K_KI_R = (0.0_qs)
      ELSE
         K_KI = (0.0_qs, 0.0_qs)
      END IF

      IF (.not. CCMP2) THEN

         DO KC = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KC) == 0)) CYCLE
            CALL BCAST2ALL_T2(WDES, KC, T2_TMP, T2_TMP_R)
            DO KL = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KL) == 0)) CYCLE
               CALL BCAST2ALL_FTOD_IA(WDES, KL, 2, PW_IA_TMP, OC_IA_TMP)
               DO KI = 1, MY_NKPTS

                  KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KI)) - &
                                           WDES%VKPT(:, KC), KPOINTS_FULL)
                  KD = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) + &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)
                  KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KI)) - &
                                            WDES%VKPT(:, KD), KPOINTS_FULL)
                  IF (LORBREAL) THEN
                     VVOO_R = (0.0_qs)
                  ELSE
                     VVOO = zero
                  END IF

                  DO NK = 1, VBMAX
                  DO NL = 1, VBMAX

                     CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                     PW_IA_TMP(1, 1, NL, RKQofKQ(KQ_)), OC_IA_TMP(1, 1, NL, RKQofKQ(KQ_)), &
                                     (NUNOCC), &
                                     FTOD_PW_IA(1, 1, NK, RKQofKQ(KQ_), KI, 1), FTOD_OC_IA(1, 1, NK, RKQofKQ(KQ_), KI, 1), &
                                     (NUNOCC), VV(1, 1), VV_R(1, 1), zero)

                     IF (LORBREAL) THEN
                        VVOO_R(:, :, NL, NK) = (VV_R(:, :))
                     ELSE
                        VVOO(:, :, NL, NK) = CONJG(VV(:, :))
                     END IF

                     IF (LDISTING) THEN
                        CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                        FTOD_PW_IA(1, 1, NK, RKQofKQ(KQ), KI, 1), FTOD_OC_IA(1, 1, NK, RKQofKQ(KQ), KI, 1), &
                                        (NUNOCC), &
                                        PW_IA_TMP(1, 1, NL, RKQofKQ(KQ)), OC_IA_TMP(1, 1, NL, RKQofKQ(KQ)), &
                                        (NUNOCC), VVOO(1, 1, NL, NK), VVOO_R(1, 1, NL, NK), (-0.5_q, 0.0_q))
                        IF (LORBREAL) THEN
                           VVOO_R(:, :, NL, NK) = VVOO_R(:, :, NL, NK)*0.5_q
                        ELSE
                           VVOO(:, :, NL, NK) = VVOO(:, :, NL, NK)*(0.5_q, 0.0_q)
                        END IF
                     ELSE
                        CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                        FTOD_PW_IA(1, 1, NK, RKQofKQ(KQ), KI, 1), FTOD_OC_IA(1, 1, NK, RKQofKQ(KQ), KI, 1), &
                                        (NUNOCC), &
                                        PW_IA_TMP(1, 1, NL, RKQofKQ(KQ)), OC_IA_TMP(1, 1, NL, RKQofKQ(KQ)), &
                                        (NUNOCC), VVOO(1, 1, NL, NK), VVOO_R(1, 1, NL, NK), (-0.5_q, 0.0_q))
                     END IF

                  END DO !NL
                  END DO !NK
                  IF (LORBREAL) THEN
                     VVOO_S_R(:, :, :, :) = T2_TMP_R(:, :, :, :, KPTS_MKPTS(KI), KL)
                     IF ((SINGLES)) THEN
                        IF (KC == KPTS_MKPTS(KI)) THEN
                           DO NI = 1, VBMAX
                           DO NL = 1, VBMAX
                           DO ND = 1, (NUNOCC)
                              VVOO_S_R(:, ND, NI, NL) = VVOO_S_R(:, ND, NI, NL) + &
                                                        T1_R(:, NI, KPTS_MKPTS(KI))*T1_R(ND, NL, KL)
                           END DO
                           END DO
                           END DO
                        END IF
                     END IF
                  ELSE
                     VVOO_S(:, :, :, :) = T2_TMP(:, :, :, :, KPTS_MKPTS(KI), KL)
                     IF ((SINGLES)) THEN
                        IF (KC == KPTS_MKPTS(KI)) THEN
                           DO NI = 1, VBMAX
                           DO NL = 1, VBMAX
                           DO ND = 1, (NUNOCC)
                              VVOO_S(:, ND, NI, NL) = VVOO_S(:, ND, NI, NL) + &
                                                      T1(:, NI, KPTS_MKPTS(KI))*T1(ND, NL, KL)
                           END DO
                           END DO
                           END DO
                        END IF
                     END IF
                  END IF
                  CALL SORT_V1V2O1O2_V1V2O2O1(WDES, VVOO_S, VVOO2_S, VVOO_S_R, VVOO2_S_R)
                  IF (LORBREAL) THEN
                     VVOO_S_R = VVOO_R
                     CALL SGEMM('t', 'n', VBMAX, VBMAX, (NUNOCC)*(NUNOCC)*VBMAX, &
                                REAL(KW2, kind=qs), VVOO_S_R(1, 1, 1, 1), (NUNOCC)*(NUNOCC)*VBMAX, &
                                VVOO2_S_R(1, 1, 1, 1), (NUNOCC)*(NUNOCC)*VBMAX, &
                                (1.0_qs), K_KI_R(1, 1, KPTS_MKPTS(KI)), VBMAX)
                  ELSE
                     VVOO_S = VVOO
                     CALL CGEMM('t', 'n', VBMAX, VBMAX, (NUNOCC)*(NUNOCC)*VBMAX, &
                                KW2, VVOO_S(1, 1, 1, 1), (NUNOCC)*(NUNOCC)*VBMAX, &
                                VVOO2_S(1, 1, 1, 1), (NUNOCC)*(NUNOCC)*VBMAX, &
                                (1.0_qs, 0.0_qs), K_KI(1, 1, KPTS_MKPTS(KI)), VBMAX)
                  END IF
               END DO
            END DO
         END DO

      END IF

      IF (LORBREAL) THEN
         CALLMPI(M_sum_single(WDES%COMM, K_KI_R, SIZE(K_KI_R)))
         IF (.not. CANONICAL) K_KI_R = K_KI_R + F_KI_R
      ELSE
         CALLMPI(M_sum_single(WDES%COMM, K_KI, 2*SIZE(K_KI)))
         IF (.not. CANONICAL) K_KI = K_KI + F_KI
      END IF

   END SUBROUTINE FORM_K_KI

!***********************************************************************
!
! The following routine computes an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE FORM_K_AC(WGW, WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      INTEGER :: KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
                 KC, KB_, NC, ND, KD
      COMPLEX(qs) :: KW

      KW = KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1)
      IF (LORBREAL) THEN
         K_AC_R = (0.0_qs)
      ELSE
         K_AC = (0.0_qs, 0.0_qs)
      END IF

      IF (.not. CCMP2) THEN

         DO KD = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KD) == 0)) CYCLE
            CALL BCAST2ALL_FTOD_AI(WDES, KD, 2, PW_AI_TMP, OC_AI_TMP)
            DO KL = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KL) == 0)) CYCLE
               DO KC = 1, MY_NKPTS

                  KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) - &
                                           WDES%VKPT(:, KD), KPOINTS_FULL)
                  KK = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KC)) - &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)
                  KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KC)) - &
                                            WDES%VKPT(:, KL), KPOINTS_FULL)
                  IF (LORBREAL) THEN
                     OOVV_R = (0.0_qs)
                  ELSE
                     OOVV = zero
                  END IF

                  DO NC = 1, (NUNOCC)
                  DO ND = 1, (NUNOCC)

                     CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                     FTOD_PW_AI(1, 1, NC, RKQofKQ(KQ), KC, 1), FTOD_OC_AI(1, 1, NC, RKQofKQ(KQ), KC, 1), &
                                     (VBMAX), &
                                     PW_AI_TMP(1, 1, ND, RKQofKQ(KQ)), OC_AI_TMP(1, 1, ND, RKQofKQ(KQ)), &
                                     (VBMAX), OOVV(1, 1, ND, NC), OOVV_R(1, 1, ND, NC), zero)

                     CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                     PW_AI_TMP(1, 1, ND, RKQofKQ(KQ_)), OC_AI_TMP(1, 1, ND, RKQofKQ(KQ_)), &
                                     (VBMAX), &
                                     FTOD_PW_AI(1, 1, NC, RKQofKQ(KQ_), KC, 1), FTOD_OC_AI(1, 1, NC, RKQofKQ(KQ_), KC, 1), &
                                     (VBMAX), TE4O(1, 1), TE4O_R(1, 1), zero)

                     IF (LORBREAL) THEN
                        IF (LDISTING) THEN
                           OOVV_R(:, :, ND, NC) = (-(OOVV_R(:, :, ND, NC))*(2.0_q) + (TE4O(:, :)))*0.5_q
                        ELSE
                           OOVV_R(:, :, ND, NC) = (-(OOVV_R(:, :, ND, NC))*(2.0_q) + (TE4O_R(:, :)))
                        END IF
                     ELSE
                        IF (LDISTING) THEN
                           OOVV(:, :, ND, NC) = (-CONJG(OOVV(:, :, ND, NC))*(2.0_q, 0.0_q) + (TE4O(:, :)))*0.5_q
                        ELSE
                           OOVV(:, :, ND, NC) = (-CONJG(OOVV(:, :, ND, NC))*(2.0_q, 0.0_q) + (TE4O(:, :)))
                        END IF
                     END IF

                  END DO !ND
                  END DO !NC

                  IF (LORBREAL) THEN
                     VVOO_S_R(:, :, :, :) = T2_R(:, :, :, :, KK, KL, KC)
                     IF ((SINGLES)) THEN
                        IF (KD == KL) THEN
                           DO NK = 1, VBMAX
                           DO NL = 1, VBMAX
                           DO ND = 1, (NUNOCC)
                              VVOO_S_R(:, ND, NK, NL) = VVOO_S_R(:, ND, NK, NL) - &
                                                        T1_R(:, NK, KPTS_MKPTS(KC))*T1_R(ND, NL, KL)
                           END DO
                           END DO
                           END DO
                        END IF
                     END IF
                  ELSE
                     VVOO_S(:, :, :, :) = T2(:, :, :, :, KK, KL, KC)
                     IF ((SINGLES)) THEN
                        IF (KD == KL) THEN
                           DO NK = 1, VBMAX
                           DO NL = 1, VBMAX
                           DO ND = 1, (NUNOCC)
                              VVOO_S(:, ND, NK, NL) = VVOO_S(:, ND, NK, NL) - &
                                                      T1(:, NK, KPTS_MKPTS(KC))*T1(ND, NL, KL)
                           END DO
                           END DO
                           END DO
                        END IF
                     END IF
                  END IF

                  CALL SORT_V1V2O1O2_O1O2V2V1(WDES, VVOO_S, OOVV_S, VVOO_S_R, OOVV_S_R)
                  IF (LORBREAL) THEN
                     OOVV2_S_R = OOVV_R
                     CALL SGEMM('t', 'n', (NUNOCC), (NUNOCC), &
                                (NUNOCC)*VBMAX*VBMAX, KW, OOVV_S_R(1, 1, 1, 1), (NUNOCC)*VBMAX*VBMAX, &
                                OOVV2_S_R(1, 1, 1, 1), (NUNOCC)*VBMAX*VBMAX, &
                                (1.0_qs), K_AC_R(1, 1, KPTS_MKPTS(KC)), (NUNOCC))
                  ELSE
                     OOVV2_S = OOVV
                     CALL CGEMM('t', 'n', (NUNOCC), (NUNOCC), &
                                (NUNOCC)*VBMAX*VBMAX, KW, OOVV_S(1, 1, 1, 1), (NUNOCC)*VBMAX*VBMAX, &
                                OOVV2_S(1, 1, 1, 1), (NUNOCC)*VBMAX*VBMAX, &
                                (1.0_qs, 0.0_qs), K_AC(1, 1, KPTS_MKPTS(KC)), (NUNOCC))
                  END IF

               END DO
            END DO
         END DO

      END IF

      IF (LORBREAL) THEN
         CALL M_sum_single(WDES%COMM, K_AC_R, SIZE(K_AC_R))
         IF (.not. CANONICAL) K_AC_R = K_AC_R + F_BA_R
      ELSE
         CALL M_sum_single(WDES%COMM, K_AC, 2*SIZE(K_AC))
         IF (.not. CANONICAL) K_AC = K_AC + F_BA
      END IF

   END SUBROUTINE FORM_K_AC

!***********************************************************************
!
! The following routine computes an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE FORM_K_KC(WGW, WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      INTEGER :: KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
                 KC, KB_, NC, ND, KD
      COMPLEX(qs) :: KW

      KW = KPOINTS_FULL%WTKPT(1)
      IF (LORBREAL) THEN
         K_KC_R = (0.0_qs)
      ELSE
         K_KC = (0.0_qs, 0.0_qs)
      END IF

      IF (.not. CCMP2) THEN
      DO KL = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KL) == 0)) CYCLE
         CALL BCAST2ALL_FTOD_AI(WDES, KL, 2, PW_AI_TMP, OC_AI_TMP)
         DO KC = 1, MY_NKPTS

            KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KC)) - &
                                     WDES%VKPT(:, KPTS_MKPTS(KC)), KPOINTS_FULL)
            KD = KL
            KK = KPTS_MKPTS(KC)
            KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) - &
                                      WDES%VKPT(:, KPTS_MKPTS(KC)), KPOINTS_FULL)
            IF (LORBREAL) THEN
               OOVV_S = (0.0_qs)
            ELSE
               OOVV = zero
            END IF

            DO NC = 1, (NUNOCC)
            DO ND = 1, (NUNOCC)

               CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                               FTOD_PW_AI(1, 1, NC, RKQofKQ(KQ), KC, 1), FTOD_OC_AI(1, 1, NC, RKQofKQ(KQ), KC, 1), &
                               (VBMAX), &
                               PW_AI_TMP(1, 1, ND, RKQofKQ(KQ)), OC_AI_TMP(1, 1, ND, RKQofKQ(KQ)), &
                               (VBMAX), OOVV(1, 1, ND, NC), OOVV_R(1, 1, ND, NC), zero)

               CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                               PW_AI_TMP(1, 1, ND, RKQofKQ(KQ_)), OC_AI_TMP(1, 1, ND, RKQofKQ(KQ_)), &
                               (VBMAX), &
                               FTOD_PW_AI(1, 1, NC, RKQofKQ(KQ_), KC, 1), FTOD_OC_AI(1, 1, NC, RKQofKQ(KQ_), KC, 1), &
                               (VBMAX), TE4O(1, 1), TE4O_R(1, 1), zero)
               IF (LORBREAL) THEN
                  OOVV_R(:, :, ND, NC) = (OOVV_R(:, :, ND, NC))*(2.0_q) - (TE4O_R(:, :))
               ELSE
                  OOVV(:, :, ND, NC) = CONJG(OOVV(:, :, ND, NC))*(2.0_q, 0.0_q) - (TE4O(:, :))
               END IF

            END DO !ND
            END DO !NC

            !klcd -> dlkc
            CALL SORT_O1O2V1V2_V2O2O1V1(WDES, OOVV_S, VOOV_S, OOVV_S_R, VOOV_S_R)
            IF (LORBREAL) THEN
               CALL SGEMM('t', 'n', (VBMAX)*(NUNOCC), 1, &
                          (NUNOCC)*VBMAX, KW, VOOV_S_R(1, 1, 1, 1), (NUNOCC)*VBMAX, &
                          T1_R(1, 1, KPTS_MKPTS(KC)), (NUNOCC)*VBMAX, &
                          (1.0_qs), K_KC_R(1, 1, KPTS_MKPTS(KC)), (VBMAX)*(NUNOCC))
            ELSE
               CALL CGEMM('t', 'n', (VBMAX)*(NUNOCC), 1, &
                          (NUNOCC)*VBMAX, KW, VOOV_S(1, 1, 1, 1), (NUNOCC)*VBMAX, &
                          T1(1, 1, KPTS_MKPTS(KC)), (NUNOCC)*VBMAX, &
                          (1.0_qs, 0.0_qs), K_KC(1, 1, KPTS_MKPTS(KC)), (VBMAX)*(NUNOCC))
            END IF

         END DO
      END DO
      END IF

      IF (LORBREAL) THEN
         CALL M_sum_single(WDES%COMM, K_KC_R, SIZE(K_KC_R))
         IF (.not. CANONICAL) THEN
         DO KI = 1, WDES%NKPTS
            K_KC_R(:, :, KI) = K_KC_R(:, :, KI) + F_KC_R(:, :, KI)
         END DO
         END IF
      ELSE
         CALL M_sum_single(WDES%COMM, K_KC, 2*SIZE(K_KC))
         IF (.not. CANONICAL) THEN
         DO KI = 1, WDES%NKPTS
            K_KC(:, :, KI) = K_KC(:, :, KI) + F_KC(:, :, KI)
         END DO
         END IF
      END IF

   END SUBROUTINE FORM_K_KC

!***********************************************************************
!
! The following routine computes an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE FORM_L_KI(WGW, WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      INTEGER :: KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
                 KC, KB_, KD, NC
      COMPLEX(qs) :: KW

      IF (LORBREAL) THEN
         L_KI_R = (0.0_qs)
      ELSE
         L_KI = (0.0_qs, 0.0_qs)
      END IF

      IF (SINGLES) THEN

         KW = (1.0_qs, 0.0_qs)*KPOINTS_FULL%WTKPT(1)
         IF (LORBREAL) THEN
            T1_T_R = (0.0_qs)
         ELSE
            T1_T = (0.0_qs, 0.0_qs)
         END IF

         DO KC = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KC) == 0)) CYCLE
            CALL BCAST2ALL_FTOD_AI(WDES, KC, 2, PW_AI_TMP, OC_AI_TMP)
            DO KI = 1, MY_NKPTS

               KK = KPTS_MKPTS(KI)

               KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KI)) - &
                                        WDES%VKPT(:, KK), KPOINTS_FULL)
               KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KK) - &
                                         WDES%VKPT(:, KC), KPOINTS_FULL)
               KL = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KC) + &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)

               DO NI = 1, VBMAX
               DO NC = 1, (NUNOCC)

                  CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                  FTOD_PW_IJ(1, 1, NI, RKQofKQ(KQ), KI, 1), FTOD_OC_IJ(1, 1, NI, RKQofKQ(KQ), KI, 1), &
                                  (VBMAX), &
                                  PW_AI_TMP(1, 1, NC, RKQofKQ(KQ)), OC_AI_TMP(1, 1, NC, RKQofKQ(KQ)), &
                                  (VBMAX), OOOV(1, 1, NI, NC), OOOV_R(1, 1, NI, NC), zero)

                  IF (LORBREAL) THEN
                     OOOV_R(:, :, NI, NC) = -(2.0_q)*(OOOV_R(:, :, NI, NC))
                  ELSE
                     OOOV(:, :, NI, NC) = -(2.0_q, 0.0_q)*CONJG(OOOV(:, :, NI, NC))
                  END IF

                  CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                  PW_AI_TMP(1, 1, NC, RKQofKQ(KQ_)), OC_AI_TMP(1, 1, NC, RKQofKQ(KQ_)), &
                                  (VBMAX), &
                                  FTOD_PW_IJ(1, 1, NI, RKQofKQ(KQ_), KI, 1), FTOD_OC_IJ(1, 1, NI, RKQofKQ(KQ_), KI, 1), &
                                  (VBMAX), OOOV(1, 1, NI, NC), OOOV_R(1, 1, NI, NC), one)

               END DO
               END DO
               IF (LORBREAL) THEN
                  OOOV_S_R(:, :, :, :) = (-OOOV_R(:, :, :, :)) !klic
               ELSE
                  OOOV_S(:, :, :, :) = (-OOOV(:, :, :, :)) !klic
               END IF

               DO NI = 1, VBMAX
               DO NK = 1, VBMAX
                  DO NL = 1, VBMAX
                     IF (LORBREAL) THEN
                        VO_S_R(:, NL) = OOOV_S_R(NK, NL, NI, :)
                     ELSE
                        VO_S(:, NL) = OOOV_S(NK, NL, NI, :)
                     END IF
                  END DO
                  IF (LORBREAL) THEN
                     CALL SGEMM('t', 'n', 1, 1, (VBMAX)*(NUNOCC), &
                                KW, VO_S_R(1, 1), (VBMAX)*(NUNOCC), &
                                T1_R(1, 1, KC), (VBMAX)*(NUNOCC), &
                                (1._qs), L_KI_R(NK, NI, KPTS_MKPTS(KI)), 1)
                  ELSE
                     CALL CGEMM('t', 'n', 1, 1, (VBMAX)*(NUNOCC), &
                                KW, VO_S(1, 1), (VBMAX)*(NUNOCC), &
                                T1(1, 1, KC), (VBMAX)*(NUNOCC), &
                                (1._qs, 0._qs), L_KI(NK, NI, KPTS_MKPTS(KI)), 1)
                  END IF
               END DO
               END DO

            END DO
         END DO

         IF (LORBREAL) THEN
            CALL M_sum_single(WDES%COMM, L_KI_R, SIZE(L_KI_R))
         ELSE
            CALL M_sum_single(WDES%COMM, L_KI, 2*SIZE(L_KI))
         END IF

      END IF

      IF (LORBREAL) THEN
         L_KI_R = L_KI_R + K_KI_R
         IF (.not. CANONICAL) THEN

            DO KI = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
               OV_S_R = F_KC_R(:, :, RKIofKI(KI))
               CALL SORT_O1V1_V1O1(WDES, OV_S, VO_S, OV_S_R, VO_S_R)

               CALL SGEMM('C', 'n', (VBMAX), (VBMAX), &
                          (NUNOCC), (1._qs), VO_S_R(1, 1), (NUNOCC), &
                          T1_R(1, 1, RKIofKI(KI)), (NUNOCC), &
                          (1._qs), L_KI_R(1, 1, RKIofKI(KI)), (VBMAX))

            END DO
         END IF
      ELSE
         L_KI = L_KI + K_KI
         IF (.not. CANONICAL) THEN

            DO KI = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
               OV_S = F_KC(:, :, KI)
               CALL SORT_O1V1_V1O1(WDES, OV_S, VO_S, OV_S_R, VO_S_R)

               CALL CGEMM('C', 'n', (VBMAX), (VBMAX), &
                          (NUNOCC), (1._qs, 0._qs), VO_S(1, 1), (NUNOCC), &
                          T1(1, 1, RKIofKI(KI)), (NUNOCC), &
                          (1._qs, 0._qs), L_KI(1, 1, RKIofKI(KI)), (VBMAX))

            END DO
         END IF
      END IF

   END SUBROUTINE FORM_L_KI

!***********************************************************************
!
! The following routine computes an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE FORM_L_AC(WGW, WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      INTEGER :: KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
                 KC, KB_, NC, ND, KD
      COMPLEX(qs) :: KW

      IF (LORBREAL) THEN
         L_AC_R = (0.0_qs)
      ELSE
         L_AC = (0.0_qs, 0.0_qs)
      END IF

      IF (SINGLES) THEN

         KW = (1.0_qs, 0.0_qs)*KPOINTS_FULL%WTKPT(1)
         DO KD = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KD) == 0)) CYCLE
            IF (LORBREAL) THEN
               VO_S_R(:, :) = T1_R(:, :, RKIofKI(KD))
            ELSE
               VO_S(:, :) = T1(:, :, RKIofKI(KD))
            END IF
            CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)
            IF (LORBREAL) THEN
               OVVO_S_R(:, :, 1, 1) = OV_S_R(:, :)
            ELSE
               OVVO_S(:, :, 1, 1) = OV_S(:, :)
            END IF
            CALL BCAST2ALL_FTOD_AB(WDES, KD, 2, PW_AB_TMP, OC_AB_TMP)
            CALL BCAST2ALL_FTOD_AI(WDES, KD, 1, PW_AI_TMP, OC_AI_TMP)
            DO KC = 1, MY_NKPTS

               KA = KPTS_MKPTS(KC)

               KQ = KPOINT_IN_FULL_GRID(-WDES%VKPT(:, KPTS_MKPTS(KC)) + &
                                        WDES%VKPT(:, KA), KPOINTS_FULL)
               KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KA) - &
                                         WDES%VKPT(:, KD), KPOINTS_FULL)
               KK = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KD) - &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)

               DO ND = 1, (NUNOCC)
               DO NC = 1, (NUNOCC)

                  CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                  PW_AI_TMP(1, 1, ND, RKQofKQ(KQ)), OC_AI_TMP(1, 1, ND, RKQofKQ(KQ)), &
                                  (VBMAX), &
                                  FTOD_PW_AB(1, 1, NC, RKQofKQ(KQ), KC, 2), FTOD_OC_AB(1, 1, NC, RKQofKQ(KQ), KC, 2), &
                                  (NUNOCC), OVVV(1, 1, NC, ND), OVVV_R(1, 1, NC, ND), zero)
                  IF (LORBREAL) THEN
                     OVVV_R(:, :, NC, ND) = -(2.0_q)*(OVVV_R(:, :, NC, ND))
                  ELSE
                     OVVV(:, :, NC, ND) = -(2.0_q, 0.0_q)*(OVVV(:, :, NC, ND))
                  END IF

                  CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                  FTOD_PW_AI(1, 1, NC, RKQofKQ(KQ_), KC, 1), FTOD_OC_AI(1, 1, NC, RKQofKQ(KQ_), KC, 1), &
                                  (VBMAX), &
                                  PW_AB_TMP(1, 1, ND, RKQofKQ(KQ_)), OC_AB_TMP(1, 1, ND, RKQofKQ(KQ_)), &
                                  (NUNOCC), OVVV(1, 1, NC, ND), OVVV_R(1, 1, NC, ND), one)

               END DO
               END DO
               !ovvv : kacd

               DO NA = 1, (NUNOCC)
               DO NC = 1, (NUNOCC)
                  IF (LORBREAL) THEN
                     OV_S_R(:, :) = -(OVVV_R(:, NA, NC, :))
                  ELSE
                     OV_S(:, :) = -CONJG(OVVV(:, NA, NC, :))
                  END IF
                  IF (LORBREAL) THEN
                     CALL SGEMM('t', 'n', 1, 1, (VBMAX)*(NUNOCC), &
                                REAL(KW, kind=qs), OV_S_R(1, 1), (VBMAX)*(NUNOCC), &
                                OVVO_S_R(1, 1, 1, 1), (VBMAX)*(NUNOCC), &
                                (1._qs), L_AC_R(NA, NC, RKIofKI(KA)), 1)
                  ELSE
                     CALL CGEMM('t', 'n', 1, 1, (VBMAX)*(NUNOCC), &
                                KW, OV_S(1, 1), (VBMAX)*(NUNOCC), &
                                OVVO_S(1, 1, 1, 1), (VBMAX)*(NUNOCC), &
                                (1._qs, 0._qs), L_AC(NA, NC, RKIofKI(KA)), 1)
                  END IF
               END DO
               END DO

            END DO
         END DO

         IF (LORBREAL) THEN
            CALL M_sum_single(WDES%COMM, L_AC_R, SIZE(L_AC_R))
         ELSE
            CALL M_sum_single(WDES%COMM, L_AC, 2*SIZE(L_AC))
         END IF

      END IF

      IF (LORBREAL) THEN
         L_AC_R = L_AC_R + K_AC_R
      ELSE
         L_AC = L_AC + K_AC
      END IF

      IF (.not. CANONICAL) THEN

         IF (LORBREAL) THEN
            DO KI = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
               VO_S_R(:, :) = T1_R(:, :, RKIofKI(KI))
               CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)
               OVVO_S_R(:, :, 1, 1) = OV_S_R(:, :)

               CALL SGEMM('t', 'n', (NUNOCC), (NUNOCC), &
                          (VBMAX), (-1._qs), OVVO_S_R(1, 1, 1, 1), (VBMAX), &
                          F_KC_R(1, 1, RKIofKI(KI)), (VBMAX), &
                          (1._qs), L_AC_R(1, 1, RKIofKI(KI)), (NUNOCC))

            END DO
         ELSE
            DO KI = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
               VO_S(:, :) = T1(:, :, RKIofKI(KI))
               CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)
               OVVO_S(:, :, 1, 1) = OV_S(:, :)

               CALL CGEMM('t', 'n', (NUNOCC), (NUNOCC), &
                          (VBMAX), (-1._qs, 0._qs), OVVO_S(1, 1, 1, 1), (VBMAX), &
                          F_KC(1, 1, RKIofKI(KI)), (VBMAX), &
                          (1._qs, 0._qs), L_AC(1, 1, RKIofKI(KI)), (NUNOCC))

            END DO
         END IF
      END IF

   END SUBROUTINE FORM_L_AC

!***********************************************************************
!
! The following routine contracts an important intermediate and adds the
! contribution to the new guess for the CCSD amplitudes.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE CONTR_L_AC_T2(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB

      IF (LORBREAL) THEN

         DO KA = 1, MY_NKPTS
            DO KI = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
               DO KJ = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE

                  VV2_S_R(:, :) = L_AC_R(:, :, KPTS_MKPTS(KA))
                  CALL SORT_V1V2_V2V1(WDES, VV2_S, VV_S, VV2_S_R, VV_S_R)

                  CALL SGEMM('t', 'n', (NUNOCC), (VBMAX*VBMAX)*(NUNOCC), &
                             (NUNOCC), (1._qs), VV_S_R(1, 1), (NUNOCC), &
                             T2_R(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KJ), KA), (NUNOCC), &
                             (1._qs), T2_N_R(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KJ), KA), (NUNOCC))
               END DO
            END DO
         END DO

      ELSE

         DO KA = 1, MY_NKPTS
            DO KI = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
               DO KJ = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE

                  VV2_S(:, :) = L_AC(:, :, KPTS_MKPTS(KA))
                  CALL SORT_V1V2_V2V1(WDES, VV2_S, VV_S, VV2_S_R, VV_S_R)

                  CALL CGEMM('t', 'n', (NUNOCC), (VBMAX*VBMAX)*(NUNOCC), &
                             (NUNOCC), (1._qs, 0._qs), VV_S(1, 1), (NUNOCC), &
                             T2(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KJ), KA), (NUNOCC), &
                             (1._qs, 0._qs), T2_N(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KJ), KA), (NUNOCC))
               END DO
            END DO
         END DO

      END IF

   END SUBROUTINE CONTR_L_AC_T2

!***********************************************************************
!
! The following routine contracts an important intermediate and adds the
! contribution to the new guess for the CCSD amplitudes.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE CONTR_L_KI_T2(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB

      IF (LORBREAL) THEN

         DO KA = 1, MY_NKPTS
            DO KI = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
               DO KJ = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE

                  VVOO_S_R(:, :, :, :) = T2_R(:, :, :, :, RKIofKI(KI), RKIofKI(KJ), KA)
                  CALL SORT_V1V2O1O2_O1V1V2O2(WDES, VVOO_S, OVVO_S, VVOO_S_R, OVVO_S_R)

                  CALL SGEMM('t', 'n', (NUNOCC)*(NUNOCC)*VBMAX, VBMAX, &
                             (VBMAX), (1._qs), OVVO_S_R(1, 1, 1, 1), (VBMAX), &
                             L_KI_R(1, 1, RKIofKI(KI)), (VBMAX), &
                             (0._qs), VVOO_S_R(1, 1, 1, 1), (NUNOCC)*(NUNOCC)*VBMAX)

                  CALL SORT_V1V2O1O2_V1V2O2O1(WDES, VVOO_S, VVOO2_S, VVOO_S_R, VVOO2_S_R)

         T2_N_R(:, :, :, :, RKIofKI(KI), RKIofKI(KJ), KA) = T2_N_R(:, :, :, :, RKIofKI(KI), RKIofKI(KJ), KA) - VVOO2_S_R(:, :, :, :)

               END DO
            END DO
         END DO

      ELSE

         DO KA = 1, MY_NKPTS
            DO KI = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
               DO KJ = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE

                  VVOO_S(:, :, :, :) = T2(:, :, :, :, RKIofKI(KI), RKIofKI(KJ), KA)
                  CALL SORT_V1V2O1O2_O1V1V2O2(WDES, VVOO_S, OVVO_S, VVOO_S_R, OVVO_S_R)

                  CALL CGEMM('t', 'n', (NUNOCC)*(NUNOCC)*VBMAX, VBMAX, &
                             (VBMAX), (1._qs, 0._qs), OVVO_S(1, 1, 1, 1), (VBMAX), &
                             L_KI(1, 1, RKIofKI(KI)), (VBMAX), &
                             (0._qs, 0._qs), VVOO_S(1, 1, 1, 1), (NUNOCC)*(NUNOCC)*VBMAX)

                  CALL SORT_V1V2O1O2_V1V2O2O1(WDES, VVOO_S, VVOO2_S, VVOO_S_R, VVOO2_S_R)

               T2_N(:, :, :, :, RKIofKI(KI), RKIofKI(KJ), KA) = T2_N(:, :, :, :, RKIofKI(KI), RKIofKI(KJ), KA) - VVOO2_S(:, :, :, :)

               END DO
            END DO
         END DO

      END IF

   END SUBROUTINE CONTR_L_KI_T2

!***********************************************************************
!
! The following routine computes an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE FORM_CHI_KLIJ(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, KC, ND
      COMPLEX(qs) :: KW

      IF (LORBREAL) THEN

         DO KL = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KL) == 0)) CYCLE
            CALL BCAST2ALL_FTOD_IJ(WDES, KL, 1, PW_IJ_TMP, OC_IJ_TMP)
            DO KK = 1, MY_NKPTS
            DO KQ = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0)) CYCLE
               KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KK)) + &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)
               KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) - &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)

               CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                               PW_IJ_TMP(1, 1, 1, RKQofKQ(KQ)), OC_IJ_TMP(1, 1, 1, RKQofKQ(KQ)), &
                               (VBMAX)*VBMAX, &
                               FTOD_PW_IJ(1, 1, 1, RKQofKQ(KQ), KK, 2), FTOD_OC_IJ(1, 1, 1, RKQofKQ(KQ), KK, 2), &
                               (VBMAX)*VBMAX, JLIK(1, 1, 1, 1), JLIK_R(1, 1, 1, 1), zero)
               CALL SORT_JLIK_KLIJ(WDES, JLIK, &
                                   KLIJ, JLIK_R, KLIJ_R)

               CHI_KLIJ_R(:, :, :, :, KI, KJ, KK) = KLIJ_R(:, :, :, :)

            END DO
            END DO
         END DO

         KW = KPOINTS_FULL%WTKPT(1)

         IF (.NOT. LCCD) THEN

            DO KC = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KC) == 0)) CYCLE
               CALL BCAST2ALL_T2(WDES, KC, T2_TMP, T2_TMP_R)
               DO KL = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KL) == 0)) CYCLE
                  CALL BCAST2ALL_FTOD_IA(WDES, KL, 2, PW_IA_TMP, OC_IA_TMP)
                  IF (SINGLES) CALL BCAST2ALL_FTOD_IJ(WDES, KL, 2, PW_IJ_TMP, OC_IJ_TMP)
                  DO KK = 1, MY_NKPTS

                     KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KK)) - &
                                              WDES%VKPT(:, KC), KPOINTS_FULL)

                     DO NL = 1, VBMAX
                     DO NK = 1, VBMAX
                        CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                        FTOD_PW_IA(1, 1, NK, RKQofKQ(KQ), KK, 1), FTOD_OC_IA(1, 1, NK, RKQofKQ(KQ), KK, 1), &
                                        (NUNOCC), &
                                        PW_IA_TMP(1, 1, NL, RKQofKQ(KQ)), OC_IA_TMP(1, 1, NL, RKQofKQ(KQ)), &
                                        (NUNOCC), VVOO(1, 1, NK, NL), VVOO_R(1, 1, NK, NL), (0.0_q, 0.0_q))
                     END DO
                     END DO
                     VVOO_S_R = VVOO_R

                     DO KJ = 1, WDES%NKPTS
                        IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE
                        KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) - &
                                                  WDES%VKPT(:, KJ), KPOINTS_FULL)
                        KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KK)) + &
                                                 WDES%VKPT(:, KQ_), KPOINTS_FULL)

                        IF ((SINGLES) .and. (KC == KI)) THEN
                           VVOO2_S_R(:, :, :, :) = 0.0_qs
                           IF (.not. LDISTING) THEN
                              VVOO2_S_R(:, :, :, :) = T2_TMP_R(:, :, :, :, RKIofKI(KI), RKIofKI(KJ))
                           END IF

                           DO NI = 1, VBMAX
                           DO NJ = 1, VBMAX
                           DO ND = 1, (NUNOCC)
                      VVOO2_S_R(:, ND, NI, NJ) = VVOO2_S_R(:, ND, NI, NJ) + T1_R(:, NI, RKIofKI(KI))*T1_R(ND, NJ, RKIofKI(KJ))/WTKPT
                           END DO
                           END DO
                           END DO

                           CALL SGEMM('t', 'n', (VBMAX*VBMAX), (VBMAX*VBMAX), (NUNOCC)*(NUNOCC), &
                                      REAL(KW, kind=qs), VVOO_S_R(1, 1, 1, 1), (NUNOCC)*(NUNOCC), &
                                      VVOO2_S_R(1, 1, 1, 1), (NUNOCC)*(NUNOCC), &
                                      (1.0_qs), CHI_KLIJ_R(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KJ), KK), (VBMAX*VBMAX))

                           DO NL = 1, VBMAX
                           DO NK = 1, VBMAX

                              KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KK)) - &
                                                       WDES%VKPT(:, KC), KPOINTS_FULL)

                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                              FTOD_PW_IA(1, 1, NK, RKQofKQ(KQ), KK, 1), FTOD_OC_IA(1, 1, NK, RKQofKQ(KQ), KK, 1), &
                                              (NUNOCC), &
                                              PW_IJ_TMP(1, 1, NL, RKQofKQ(KQ)), OC_IJ_TMP(1, 1, NL, RKQofKQ(KQ)), &
                                              (VBMAX), VOVO(1, 1, 1, 1), VOVO_R(1, 1, 1, 1), (0.0_q, 0.0_q))
                              VO_S_R(:, :) = VOVO_R(:, :, 1, 1)

                              DO NJ = 1, VBMAX
                              DO NI = 1, VBMAX
                                 CALL SGEMM('t', 'n', 1, 1, (NUNOCC), &
                                            (1.0_qs), VO_S_R(1, NJ), (NUNOCC), &
                                            T1_R(1, NI, RKIofKI(KI)), (NUNOCC), &
                                            (1.0_qs), CHI_KLIJ_R(NK, NL, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KK), 1)
                              END DO
                              END DO
                           END DO
                           END DO

                        ELSE

                           IF (.not. LDISTING) THEN
                              CALL SGEMM('t', 'n', (VBMAX*VBMAX), (VBMAX*VBMAX), (NUNOCC)*(NUNOCC), &
                                         KW, VVOO_S_R(1, 1, 1, 1), (NUNOCC)*(NUNOCC), &
                                         T2_TMP_R(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KJ)), (NUNOCC)*(NUNOCC), &
                                         (1.0_qs), CHI_KLIJ_R(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KJ), KK), (VBMAX*VBMAX))
                           END IF
                        END IF

                        IF ((SINGLES) .and. (KJ == KC)) THEN

                           KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KK)) - &
                                                    WDES%VKPT(:, KI), KPOINTS_FULL)
                           DO NL = 1, VBMAX
                           DO NK = 1, VBMAX
                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                              FTOD_PW_IJ(1, 1, NK, RKQofKQ(KQ), KK, 1), FTOD_OC_IJ(1, 1, NK, RKQofKQ(KQ), KK, 1), &
                                              (VBMAX), &
                                              PW_IA_TMP(1, 1, NL, RKQofKQ(KQ)), OC_IA_TMP(1, 1, NL, RKQofKQ(KQ)), &
                                              (NUNOCC), OVVO(1, 1, 1, 1), OVVO_R(1, 1, 1, 1), (0.0_q, 0.0_q))

                              OV_S_R(:, :) = OVVO_R(:, :, 1, 1)
                              CALL SORT_O1V1_V1O1(WDES, OV_S, VO_S, OV_S_R, VO_S_R)

                              DO NJ = 1, VBMAX
                              DO NI = 1, VBMAX
                                 CALL SGEMM('t', 'n', 1, 1, (NUNOCC), &
                                            (1.0_qs), VO_S_R(1, NI), (NUNOCC), &
                                            T1_R(1, NJ, RKIofKI(KJ)), (NUNOCC), &
                                            (1.0_qs), CHI_KLIJ_R(NK, NL, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KK), 1)
                              END DO
                              END DO
                           END DO
                           END DO
                        END IF

                     END DO
                  END DO
               END DO
            END DO

         END IF !LCCD

      ELSE ! not LORBREAL

         DO KL = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KL) == 0)) CYCLE
            CALL BCAST2ALL_FTOD_IJ(WDES, KL, 1, PW_IJ_TMP, OC_IJ_TMP)
            DO KK = 1, MY_NKPTS
            DO KQ = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0)) CYCLE
               KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KK)) + &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)
               KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) - &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)

               CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                               PW_IJ_TMP(1, 1, 1, RKQofKQ(KQ)), OC_IJ_TMP(1, 1, 1, RKQofKQ(KQ)), &
                               (VBMAX)*VBMAX, &
                               FTOD_PW_IJ(1, 1, 1, RKQofKQ(KQ), KK, 2), FTOD_OC_IJ(1, 1, 1, RKQofKQ(KQ), KK, 2), &
                               (VBMAX)*VBMAX, JLIK(1, 1, 1, 1), JLIK_R(1, 1, 1, 1), zero)
               CALL SORT_JLIK_KLIJ(WDES, JLIK, &
                                   KLIJ, JLIK_R, KLIJ_R)

               CHI_KLIJ(:, :, :, :, KI, KJ, KK) = KLIJ(:, :, :, :)

            END DO
            END DO
         END DO

         KW = KPOINTS_FULL%WTKPT(1)

         IF (.NOT. LCCD) THEN

            DO KC = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KC) == 0)) CYCLE
               CALL BCAST2ALL_T2(WDES, KC, T2_TMP, T2_TMP_R)
               DO KL = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KL) == 0)) CYCLE
                  CALL BCAST2ALL_FTOD_IA(WDES, KL, 2, PW_IA_TMP, OC_IA_TMP)
                  IF (SINGLES) CALL BCAST2ALL_FTOD_IJ(WDES, KL, 2, PW_IJ_TMP, OC_IJ_TMP)
                  DO KK = 1, MY_NKPTS

                     KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KK)) - &
                                              WDES%VKPT(:, KC), KPOINTS_FULL)

                     !VVOO=zero

                     DO NL = 1, VBMAX
                     DO NK = 1, VBMAX
                        CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                        FTOD_PW_IA(1, 1, NK, RKQofKQ(KQ), KK, 1), FTOD_OC_IA(1, 1, NK, RKQofKQ(KQ), KK, 1), &
                                        (NUNOCC), &
                                        PW_IA_TMP(1, 1, NL, RKQofKQ(KQ)), OC_IA_TMP(1, 1, NL, RKQofKQ(KQ)), &
                                        (NUNOCC), VVOO(1, 1, NK, NL), VVOO_R(1, 1, NK, NL), (0.0_q, 0.0_q))
                     END DO
                     END DO
                     !vvoo : cdkl
                     VVOO_S = VVOO

                     DO KJ = 1, WDES%NKPTS
                        IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE
                        KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) - &
                                                  WDES%VKPT(:, KJ), KPOINTS_FULL)
                        KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KK)) + &
                                                 WDES%VKPT(:, KQ_), KPOINTS_FULL)

                        IF ((SINGLES) .and. (KC == KI)) THEN

                           VVOO2_S(:, :, :, :) = (0.0_qs, 0.0_qs)
                           IF (.not. LDISTING) THEN
                              VVOO2_S(:, :, :, :) = T2_TMP(:, :, :, :, RKIofKI(KI), RKIofKI(KJ))
                           END IF

                           DO NI = 1, VBMAX
                           DO NJ = 1, VBMAX
                           DO ND = 1, (NUNOCC)
                              VVOO2_S(:, ND, NI, NJ) = VVOO2_S(:, ND, NI, NJ) + T1(:, NI, RKIofKI(KI))*T1(ND, NJ, RKIofKI(KJ))/WTKPT
                           END DO
                           END DO
                           END DO

                           CALL CGEMM('t', 'n', (VBMAX*VBMAX), (VBMAX*VBMAX), (NUNOCC)*(NUNOCC), &
                                      KW, VVOO_S(1, 1, 1, 1), (NUNOCC)*(NUNOCC), &
                                      VVOO2_S(1, 1, 1, 1), (NUNOCC)*(NUNOCC), &
                                      (1.0_qs, 0.0_qs), CHI_KLIJ(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KJ), KK), (VBMAX*VBMAX))

                           DO NL = 1, VBMAX
                           DO NK = 1, VBMAX

                              KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KK)) - &
                                                       WDES%VKPT(:, KC), KPOINTS_FULL)

                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                              FTOD_PW_IA(1, 1, NK, RKQofKQ(KQ), KK, 1), FTOD_OC_IA(1, 1, NK, RKQofKQ(KQ), KK, 1), &
                                              (NUNOCC), &
                                              PW_IJ_TMP(1, 1, NL, RKQofKQ(KQ)), OC_IJ_TMP(1, 1, NL, RKQofKQ(KQ)), &
                                              (VBMAX), VOVO(1, 1, 1, 1), VOVO_R(1, 1, 1, 1), (0.0_q, 0.0_q))
                              VO_S(:, :) = VOVO(:, :, 1, 1)

                              DO NJ = 1, VBMAX
                              DO NI = 1, VBMAX
                                 CALL CGEMM('t', 'n', 1, 1, (NUNOCC), &
                                            (1.0_qs, 0.0_qs), VO_S(1, NJ), (NUNOCC), &
                                            T1(1, NI, RKIofKI(KI)), (NUNOCC), &
                                            (1.0_qs, 0.0_qs), CHI_KLIJ(NK, NL, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KK), 1)
                              END DO
                              END DO
                           END DO
                           END DO

                        ELSE

                           IF (.not. LDISTING) THEN
                              CALL CGEMM('t', 'n', (VBMAX*VBMAX), (VBMAX*VBMAX), (NUNOCC)*(NUNOCC), &
                                         KW, VVOO_S(1, 1, 1, 1), (NUNOCC)*(NUNOCC), &
                                         T2_TMP(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KJ)), (NUNOCC)*(NUNOCC), &
                                         (1.0_qs, 0.0_qs), CHI_KLIJ(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KJ), KK), (VBMAX*VBMAX))
                           END IF

                        END IF

                        IF ((SINGLES) .and. (KJ == KC)) THEN

                           KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KK)) - &
                                                    WDES%VKPT(:, KI), KPOINTS_FULL)
                           DO NL = 1, VBMAX
                           DO NK = 1, VBMAX
                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                              FTOD_PW_IJ(1, 1, NK, RKQofKQ(KQ), KK, 1), FTOD_OC_IJ(1, 1, NK, RKQofKQ(KQ), KK, 1), &
                                              (VBMAX), &
                                              PW_IA_TMP(1, 1, NL, RKQofKQ(KQ)), OC_IA_TMP(1, 1, NL, RKQofKQ(KQ)), &
                                              (NUNOCC), OVVO(1, 1, 1, 1), OVVO_R(1, 1, 1, 1), (0.0_q, 0.0_q))

                              OV_S(:, :) = OVVO(:, :, 1, 1)
                              CALL SORT_O1V1_V1O1(WDES, OV_S, VO_S, OV_S_R, VO_S_R)

                              DO NJ = 1, VBMAX
                              DO NI = 1, VBMAX
                                 CALL CGEMM('t', 'n', 1, 1, (NUNOCC), &
                                            (1.0_qs, 0.0_qs), VO_S(1, NI), (NUNOCC), &
                                            T1(1, NJ, RKIofKI(KJ)), (NUNOCC), &
                                            (1.0_qs, 0.0_qs), CHI_KLIJ(NK, NL, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KK), 1)
                              END DO
                              END DO
                           END DO
                           END DO
                        END IF
                     END DO
                  END DO
               END DO
            END DO

         END IF !LCCD

      END IF

   END SUBROUTINE FORM_CHI_KLIJ

!***********************************************************************
!
! The following routine contracts an important intermediate and adds the
! contribution to the new guess for the CCSD amplitudes.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE CONTR_CHI_KLIJ_T(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB

      DO KK = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KK) == 0)) CYCLE
         CALL BCAST2ALL_CHI_KLIJ(WDES, KK, CHI_KLIJ_TMP, CHI_KLIJ_TMP_R)

         DO KB = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
            DO KA = 1, MY_NKPTS
               KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KK) &
                                        - WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
               KL = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)

               IF (LORBREAL) THEN
                  VVOO_R(:, :, :, :) = T2_R(:, :, :, :, RKIofKI(KK), RKIofKI(KL), KA)*WTKPT
                  CALL SORT_V1V2O1O2_O1O2V1V2(WDES, VVOO, OOVV, VVOO_R, OOVV_R)
                  OOVV_S_R = OOVV_R
               ELSE
                  VVOO(:, :, :, :) = T2(:, :, :, :, RKIofKI(KK), RKIofKI(KL), KA)*WTKPT
                  CALL SORT_V1V2O1O2_O1O2V1V2(WDES, VVOO, OOVV, VVOO_R, OOVV_R)
                  OOVV_S = OOVV
               END IF

               IF (SINGLES) THEN
               IF ((KK == KPTS_MKPTS(KA)) .and. (KL == KB)) THEN
               IF (LORBREAL) THEN
                  DO NK = 1, VBMAX
                  DO NL = 1, VBMAX
                  DO NA = 1, (NUNOCC)
                  DO NB = 1, (NUNOCC)
                     OOVV_S_R(NK, NL, NA, NB) = OOVV_S_R(NK, NL, NA, NB) + T1_R(NA, NK, RKIofKI(KK))*T1_R(NB, NL, RKIofKI(KL))
                  END DO
                  END DO
                  END DO
                  END DO
               ELSE
                  DO NK = 1, VBMAX
                  DO NL = 1, VBMAX
                  DO NA = 1, (NUNOCC)
                  DO NB = 1, (NUNOCC)
                     OOVV_S(NK, NL, NA, NB) = OOVV_S(NK, NL, NA, NB) + T1(NA, NK, RKIofKI(KK))*T1(NB, NL, RKIofKI(KL))
                  END DO
                  END DO
                  END DO
                  END DO
               END IF
               END IF
               END IF

               DO KQ = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0)) CYCLE
                  KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) + &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)
                  KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)

!#ifdef gammareal
!                       CALL SGEMM('t','n',(TWOE4ORBITAL_COLS*TWOE4ORBITAL_ROWS),(VBMAX*VBMAX),&
!                          (NSTRIPT*NSTRIPT),1._qs,TWOE4ORBITAL_STRIPT(1,1,1,1),(NSTRIPT*NSTRIPT),&
!                          T2AMPLITUDES_CON_STRIPT(1,1,1,1,mki),(NSTRIPT*NSTRIPT),&
!                          1._qs, T2AMPLITUDES_new(1,1,1,1,mki,mkj,kr),(TWOE4ORBITAL_COLS*TWOE4ORBITAL_ROWS))
!#else
                  IF (LORBREAL) THEN
                     CALL SGEMM('t', 'n', (NUNOCC)*(NUNOCC), (VBMAX*VBMAX), &
                                (VBMAX*VBMAX), (1._qs), OOVV_S_R(1, 1, 1, 1), (VBMAX*VBMAX), &
                                CHI_KLIJ_TMP_R(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KJ)), (VBMAX*VBMAX), &
                                (1._qs), T2_N_R(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KJ), KA), (NUNOCC)*(NUNOCC))
                  ELSE
                     CALL CGEMM('t', 'n', (NUNOCC)*(NUNOCC), (VBMAX*VBMAX), &
                                (VBMAX*VBMAX), (1._qs, 0._qs), OOVV_S(1, 1, 1, 1), (VBMAX*VBMAX), &
                                CHI_KLIJ_TMP(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KJ)), (VBMAX*VBMAX), &
                                (1._qs, 0._qs), T2_N(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KJ), KA), (NUNOCC)*(NUNOCC))
                  END IF
!#endif

               END DO
            END DO
         END DO

      END DO !KK

   END SUBROUTINE CONTR_CHI_KLIJ_T

!***********************************************************************
!
! This is the computationally most intense routine of the CCSD algorithm
! Here we compute and contract the \chi^ab_cd intermediate on-the-fly
! and add its contribution to the new guess for the T amplitudes.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************


   SUBROUTINE FORM_CONTR_ABCD(WDES, KPOINTS)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(kpoints_struct) KPOINTS
      INTEGER :: omp_get_thread_num, omp_get_num_threads, CHUNK, myid
      INTEGER KI, KJ, KA, KB, NI, NJ, NK, NL, KQ, KQ_, KC, KD, &
         NC, ND, NA, NB, RNBLOCKA, RNBLOCKB, RNA, RNB, KC_, NT, NG
      integer :: time_array1(8), time_array2(8), ems

      call BLACS_PINFO(ME, PROCS)

      IF (LORBREAL) THEN
         DO KC_ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KC_) == 0)) CYCLE
            KC = KC_ !KPOINT_IN_FULL_GRID(WDES%VKPT(:,KC_),KPOINTS_FULL)
            CALL BCAST2ALL_T2(WDES, KC, T2_TMP, T2_TMP_R)
            DO KB = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
               call date_and_time(values=time_array1)
               ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                      time_array2(8) - time_array1(8))
               call date_and_time(values=time_array2)
               CALL BCAST2ALL_FTOD_AB(WDES, KB, 2, PW_AB_TMP, OC_AB_TMP) !it wouldn't be necessary to communicate all kq to all processes here
               IF (SINGLES) CALL BCAST2ALL_FTOD_IA(WDES, KB, 2, PW_IA_TMP, OC_IA_TMP) !it wouldn't be necessary to communicate all kq to all processes here
               DO KA = 1, MY_NKPTS

                  DO NA = 1, (NUNOCC), NBLOCKAB
                     RNBLOCKA = MIN((NUNOCC) - NA + 1, NBLOCKAB)
                     DO NB = 1, (NUNOCC), NBLOCKAB
                        RNBLOCKB = MIN((NUNOCC) - NB + 1, NBLOCKAB)
!************************************************************************************************
! OMP Multithreading is experimental and was not very successfull her we no longer use it for the time being.
!************************************************************************************************
!                  CALL OMP_SET_NUM_THREADS(nthreads)

                        KQ = KPOINT_IN_FULL_GRID((WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                                  WDES%VKPT(:, KC)), KPOINTS_FULL)
                        KD = KPOINT_IN_FULL_GRID((WDES%VKPT(:, KB) + &
                                                  WDES%VKPT(:, KQ)), KPOINTS_FULL)

!$omp parallel shared (RNBLOCKB,RNBLOCKA,NGVECTOR,NHVECTOR,FTOD_PW_AB,FTOD_OC_AB,NUNOCC,PW_AB_TMP,OC_AB_TMP,VVVV_R,nthreads), private (RNB,RNA,ND,NC,NG)
!$OMP DO SCHEDULE(STATIC)
                        DO RNB = 1, RNBLOCKB
                           DO RNA = 1, RNBLOCKA

                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                           FTOD_PW_AB(1, 1, NA + RNA - 1, RKQofKQ(KQ), KA, 1), FTOD_OC_AB(1, 1, NA + RNA - 1, RKQofKQ(KQ), KA, 1), &
                                              (NUNOCC), &
                                           PW_AB_TMP(1, 1, NB + RNB - 1, RKQofKQ(KQ)), OC_AB_TMP(1, 1, NB + RNB - 1, RKQofKQ(KQ)), &
                                              (NUNOCC), VVVV(1, 1, RNA, RNB), VVVV_R(1, 1, RNA, RNB), zero)


                           END DO
                        END DO
!$OMP END DO
!$omp end parallel

!                  CALL OMP_SET_NUM_THREADS(1)

                        VVVV_S_R = VVVV_R

                        IF (SINGLES) THEN
                           DO RNA = 1, RNBLOCKA !(NUNOCC)
                           DO NK = 1, VBMAX
                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                           FTOD_PW_AB(1, 1, NA + RNA - 1, RKQofKQ(KQ), KA, 1), FTOD_OC_AB(1, 1, NA + RNA - 1, RKQofKQ(KQ), KA, 1), &
                                              (NUNOCC), &
                                              PW_IA_TMP(1, 1, NK, RKQofKQ(KQ)), OC_IA_TMP(1, 1, NK, RKQofKQ(KQ)), &
                                              (NUNOCC), VV(1, 1), VV_R(1, 1), zero)
                              OVVV_S_R(NK, :, :, NA + RNA - 1) = -VV_R(:, :)
                           END DO
                           END DO
                           VO_S_R(:, :) = T1_R(:, :, RKIofKI(KB))
                           CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)

                           DO RNA = 1, RNBLOCKA !(NUNOCC)
                           DO RNB = 1, RNBLOCKB !(NUNOCC)

                              CALL SGEMM('t', 'n', (NUNOCC)*(NUNOCC), 1, &
                                         (VBMAX), (1._qs), OVVV_S_R(1, 1, 1, NA + RNA - 1), &
                                         (VBMAX), &
                                         OV_S_R(1, NB + RNB - 1), (VBMAX), &
                                         (1._qs), VVVV_S_R(1, 1, RNA, RNB), (NUNOCC)*(NUNOCC))

                           END DO
                           END DO

                           DO RNB = 1, RNBLOCKB !(NUNOCC)
                           DO NK = 1, VBMAX
                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                              FTOD_PW_IA(1, 1, NK, RKQofKQ(KQ), KA, 1), FTOD_OC_IA(1, 1, NK, RKQofKQ(KQ), KA, 1), &
                                              (NUNOCC), &
                                           PW_AB_TMP(1, 1, NB + RNB - 1, RKQofKQ(KQ)), OC_AB_TMP(1, 1, NB + RNB - 1, RKQofKQ(KQ)), &
                                              (NUNOCC), VV(1, 1), VV_R(1, 1), zero)
                              OVVV_S_R(NK, :, :, NB + RNB - 1) = -VV_R(:, :)
                           END DO
                           END DO
                           VO_S_R(:, :) = T1_R(:, :, KPTS_MKPTS(KA))
                           CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)

                           DO RNA = 1, RNBLOCKA !(NUNOCC)
                           DO RNB = 1, RNBLOCKB !(NUNOCC)

                              CALL SGEMM('t', 'n', (NUNOCC)*(NUNOCC), 1, &
                                         (VBMAX), (1._qs), OVVV_S_R(1, 1, 1, NB + RNB - 1), &
                                         (VBMAX), &
                                         OV_S_R(1, NA + RNA - 1), (VBMAX), &
                                         (1._qs), VVVV_S_R(1, 1, RNA, RNB), (NUNOCC)*(NUNOCC))

                           END DO
                           END DO

                        END IF

                        VVVV_S_R = VVVV_S_R*KPOINTS_FULL%WTKPT(1)

!               CALL OMP_SET_NUM_THREADS(nthreads)

!$omp parallel shared (T2_N_R,T2_TMP_R,RNBLOCKA,NUNOCC,VVVV_S_R,NA,NB,KA), private (NI,NJ,myid,NC,ND,KI,KQ,KJ,RNA,RNB,VVOO_S_R)
                        IF (ALLOCATED(VVOO_S_R)) DEALLOCATE (VVOO_S_R)
                        ALLOCATE (VVOO_S_R(NUNOCC, NUNOCC, VBMAX, VBMAX))
!$OMP DO SCHEDULE(STATIC)
                        DO KI = 1, WDES%NKPTS
                           IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE

                           KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                                    WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                           KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                                    WDES%VKPT(:, KQ), KPOINTS_FULL)
                           IF (KJ > KI) CYCLE

                           VVOO_S_R(:, :, :, :) = T2_TMP_R(:, :, :, :, RKIofKI(KI), RKIofKI(KJ))
                           IF (SINGLES) THEN
                           IF ((KI == KC) .and. (KJ == KD)) THEN
                              DO NJ = 1, VBMAX
                              DO NI = 1, VBMAX
                              DO ND = 1, (NUNOCC)
                              DO NC = 1, (NUNOCC)
                        VVOO_S_R(NC, ND, NI, NJ) = VVOO_S_R(NC, ND, NI, NJ) + T1_R(NC, NI, RKIofKI(KI))*T1_R(ND, NJ, RKIofKI(KJ))/ &
                                                            WTKPT
                              END DO
                              END DO
                              END DO
                              END DO
                           END IF
                           END IF

                           DO NJ = 1, VBMAX
                           DO NI = 1, VBMAX

                              DO RNB = 1, RNBLOCKB

                                 CALL SGEMM('t', 'n', RNBLOCKA, 1, &
                                            (NUNOCC)*(NUNOCC), (1._qs), VVVV_S_R(1, 1, 1, RNB), &
                                            (NUNOCC)*(NUNOCC), &
                                            VVOO_S_R(1, 1, NI, NJ), (NUNOCC)*(NUNOCC), &
                                            (1._qs), T2_N_R(NA, NB + RNB - 1, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA), RNBLOCKA)

                              END DO
                           END DO
                           END DO
                        END DO
!$OMP END DO
!$omp end parallel
!               CALL OMP_SET_NUM_THREADS(1)
                     END DO !NA,NBLOCK
                  END DO !NB,NBLOCK
               END DO
            END DO
         END DO

      ELSE
         DO KC_ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KC_) == 0)) CYCLE
            KC = KC_ !KPOINT_IN_FULL_GRID(WDES%VKPT(:,KC_),KPOINTS_FULL)
            CALL BCAST2ALL_T2(WDES, KC, T2_TMP, T2_TMP_R)
            DO KB = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
               call date_and_time(values=time_array1)
               ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                      time_array2(8) - time_array1(8))
               call date_and_time(values=time_array2)
               CALL BCAST2ALL_FTOD_AB(WDES, KB, 2, PW_AB_TMP, OC_AB_TMP) !it wouldn't be necessary to communicate all kq to all processes here
               IF (SINGLES) CALL BCAST2ALL_FTOD_IA(WDES, KB, 2, PW_IA_TMP, OC_IA_TMP) !it wouldn't be necessary to communicate all kq to all processes here
               DO KA = 1, MY_NKPTS

                  DO NA = 1, (NUNOCC), NBLOCKAB
                     RNBLOCKA = MIN((NUNOCC) - NA + 1, NBLOCKAB)
                     DO NB = 1, (NUNOCC), NBLOCKAB
                        RNBLOCKB = MIN((NUNOCC) - NB + 1, NBLOCKAB)
!                  CALL OMP_SET_NUM_THREADS(nthreads)

                        KQ = KPOINT_IN_FULL_GRID((WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                                  WDES%VKPT(:, KC)), KPOINTS_FULL)
                        KD = KPOINT_IN_FULL_GRID((WDES%VKPT(:, KB) + &
                                                  WDES%VKPT(:, KQ)), KPOINTS_FULL)

!$omp parallel shared (RNBLOCKB,RNBLOCKA,NGVECTOR,NHVECTOR,FTOD_PW_AB,FTOD_OC_AB,NUNOCC,PW_AB_TMP,OC_AB_TMP,VVVV,nthreads), private (RNB,RNA,ND,NC,NG)
!$OMP DO SCHEDULE(STATIC)
                        DO RNB = 1, RNBLOCKB
                           DO RNA = 1, RNBLOCKA

                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                           FTOD_PW_AB(1, 1, NA + RNA - 1, RKQofKQ(KQ), KA, 1), FTOD_OC_AB(1, 1, NA + RNA - 1, RKQofKQ(KQ), KA, 1), &
                                              (NUNOCC), &
                                           PW_AB_TMP(1, 1, NB + RNB - 1, RKQofKQ(KQ)), OC_AB_TMP(1, 1, NB + RNB - 1, RKQofKQ(KQ)), &
                                              (NUNOCC), VVVV(1, 1, RNA, RNB), VVVV_R(1, 1, RNA, RNB), zero)

                           END DO
                        END DO
!$OMP END DO
!$omp end parallel

!                  CALL OMP_SET_NUM_THREADS(1)

                        VVVV_S = VVVV

                        IF (SINGLES) THEN
                           DO RNA = 1, RNBLOCKA !(NUNOCC)
                           DO NK = 1, VBMAX
                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                           FTOD_PW_AB(1, 1, NA + RNA - 1, RKQofKQ(KQ), KA, 1), FTOD_OC_AB(1, 1, NA + RNA - 1, RKQofKQ(KQ), KA, 1), &
                                              (NUNOCC), &
                                              PW_IA_TMP(1, 1, NK, RKQofKQ(KQ)), OC_IA_TMP(1, 1, NK, RKQofKQ(KQ)), &
                                              (NUNOCC), VV(1, 1), VV_R(1, 1), zero)
                              OVVV_S(NK, :, :, NA + RNA - 1) = -VV(:, :)
                           END DO
                           END DO
                           VO_S(:, :) = T1(:, :, RKIofKI(KB))
                           CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)

                           DO RNA = 1, RNBLOCKA !(NUNOCC)
                           DO RNB = 1, RNBLOCKB !(NUNOCC)

                              CALL CGEMM('t', 'n', (NUNOCC)*(NUNOCC), 1, &
                                         (VBMAX), (1._qs, 0._qs), OVVV_S(1, 1, 1, NA + RNA - 1), &
                                         (VBMAX), &
                                         OV_S(1, NB + RNB - 1), (VBMAX), &
                                         (1._qs, 0._qs), VVVV_S(1, 1, RNA, RNB), (NUNOCC)*(NUNOCC))

                           END DO
                           END DO

                           DO RNB = 1, RNBLOCKB !(NUNOCC)
                           DO NK = 1, VBMAX
                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                              FTOD_PW_IA(1, 1, NK, RKQofKQ(KQ), KA, 1), FTOD_OC_IA(1, 1, NK, RKQofKQ(KQ), KA, 1), &
                                              (NUNOCC), &
                                           PW_AB_TMP(1, 1, NB + RNB - 1, RKQofKQ(KQ)), OC_AB_TMP(1, 1, NB + RNB - 1, RKQofKQ(KQ)), &
                                              (NUNOCC), VV(1, 1), VV_R(1, 1), zero)
                              OVVV_S(NK, :, :, NB + RNB - 1) = -VV(:, :)
                           END DO
                           END DO
                           VO_S(:, :) = T1(:, :, KPTS_MKPTS(KA))
                           CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)

                           DO RNA = 1, RNBLOCKA !(NUNOCC)
                           DO RNB = 1, RNBLOCKB !(NUNOCC)

                              CALL CGEMM('t', 'n', (NUNOCC)*(NUNOCC), 1, &
                                         (VBMAX), (1._qs, 0._qs), OVVV_S(1, 1, 1, NB + RNB - 1), &
                                         (VBMAX), &
                                         OV_S(1, NA + RNA - 1), (VBMAX), &
                                         (1._qs, 0._qs), VVVV_S(1, 1, RNA, RNB), (NUNOCC)*(NUNOCC))

                           END DO
                           END DO

                        END IF

                        VVVV_S = VVVV_S*WTKPT

!               CALL OMP_SET_NUM_THREADS(nthreads)

!$omp parallel shared (T2_N,T2_TMP,RNBLOCKA,NUNOCC,VVVV_S,NA,NB,KA), private (NI,NJ,myid,NC,ND,KI,KQ,KJ,RNA,RNB,VVOO_S)
                        IF (ALLOCATED(VVOO_S)) DEALLOCATE (VVOO_S)
                        ALLOCATE (VVOO_S(NUNOCC, NUNOCC, VBMAX, VBMAX))
!$OMP DO SCHEDULE(STATIC)
                        DO KI = 1, WDES%NKPTS
                           IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE

                           KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                                    WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                           KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                                    WDES%VKPT(:, KQ), KPOINTS_FULL)
                           IF (KJ > KI) CYCLE

                           VVOO_S(:, :, :, :) = T2_TMP(:, :, :, :, RKIofKI(KI), RKIofKI(KJ))
                           IF (SINGLES) THEN
                           IF ((KI == KC) .and. (KJ == KD)) THEN
                              DO NJ = 1, VBMAX
                              DO NI = 1, VBMAX
                              DO ND = 1, (NUNOCC)
                              DO NC = 1, (NUNOCC)
                                VVOO_S(NC, ND, NI, NJ) = VVOO_S(NC, ND, NI, NJ) + T1(NC, NI, RKIofKI(KI))*T1(ND, NJ, RKIofKI(KJ))/ &
                                                          WTKPT
                              END DO
                              END DO
                              END DO
                              END DO
                           END IF
                           END IF

                           DO NJ = 1, VBMAX
                           DO NI = 1, VBMAX

                              DO RNB = 1, RNBLOCKB

                                 CALL CGEMM('t', 'n', RNBLOCKA, 1, &
                                            (NUNOCC)*(NUNOCC), (1._qs, 0._qs), VVVV_S(1, 1, 1, RNB), &
                                            (NUNOCC)*(NUNOCC), &
                                            VVOO_S(1, 1, NI, NJ), (NUNOCC)*(NUNOCC), &
                                            (1._qs, 0._qs), T2_N(NA, NB + RNB - 1, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA), RNBLOCKA)

                              END DO
                           END DO
                           END DO
                        END DO
!$OMP END DO
!$omp end parallel
!               CALL OMP_SET_NUM_THREADS(1)
                     END DO !NA,NBLOCK
                  END DO !NB,NBLOCK
               END DO
            END DO
         END DO

      END IF

   END SUBROUTINE FORM_CONTR_ABCD

!***********************************************************************
!
! This is an experimental routine needed to compute the particle-particle-ladder
! contribution to the CCSD energy.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE FORM_CONTR_ABCD_PPL(WDES, KPOINTS)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(kpoints_struct) KPOINTS
      INTEGER :: omp_get_thread_num, omp_get_num_threads, CHUNK, myid
      INTEGER KI, KJ, KA, KB, NI, NJ, NK, NL, KQ, KQ_, KC, KD, &
         NC, ND, NA, NB, RNBLOCKA, RNBLOCKB, RNA, RNB, KC_, NT, NG
      integer :: time_array1(8), time_array2(8), ems

      call BLACS_PINFO(ME, PROCS)

      IF (LORBREAL) THEN
         DO KC_ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KC_) == 0)) CYCLE
            KC = KC_ !KPOINT_IN_FULL_GRID(WDES%VKPT(:,KC_),KPOINTS_FULL)
            CALL BCAST2ALL_T2(WDES, KC, T2_TMP, T2_TMP_R)
            DO KB = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
               call date_and_time(values=time_array1)
               ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                      time_array2(8) - time_array1(8))
               call date_and_time(values=time_array2)
               CALL BCAST2ALL_FTOD_AB(WDES, KB, 2, PW_AB_TMP, OC_AB_TMP) !it wouldn't be necessary to communicate all kq to all processes here
               IF (SINGLES) CALL BCAST2ALL_FTOD_IA(WDES, KB, 2, PW_IA_TMP, OC_IA_TMP) !it wouldn't be necessary to communicate all kq to all processes here
               DO KA = 1, MY_NKPTS

                  DO NA = 1, (NUNOCC), NBLOCKAB
                     RNBLOCKA = MIN((NUNOCC) - NA + 1, NBLOCKAB)
                     DO NB = 1, (NUNOCC), NBLOCKAB
                        RNBLOCKB = MIN((NUNOCC) - NB + 1, NBLOCKAB)

                        KQ = KPOINT_IN_FULL_GRID((WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                                  WDES%VKPT(:, KC)), KPOINTS_FULL)
                        KD = KPOINT_IN_FULL_GRID((WDES%VKPT(:, KB) + &
                                                  WDES%VKPT(:, KQ)), KPOINTS_FULL)

!$omp parallel shared (RNBLOCKB,RNBLOCKA,NGVECTOR,NHVECTOR,FTOD_PW_AB,FTOD_OC_AB,NUNOCC,PW_AB_TMP,OC_AB_TMP,VVVV_R,nthreads), private (RNB,RNA,ND,NC,NG)
!$OMP DO SCHEDULE(STATIC)
                        DO RNB = 1, RNBLOCKB
                           DO RNA = 1, RNBLOCKA

                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                           FTOD_PW_AB(1, 1, NA + RNA - 1, RKQofKQ(KQ), KA, 1), FTOD_OC_AB(1, 1, NA + RNA - 1, RKQofKQ(KQ), KA, 1), &
                                              (NUNOCC), &
                                           PW_AB_TMP(1, 1, NB + RNB - 1, RKQofKQ(KQ)), OC_AB_TMP(1, 1, NB + RNB - 1, RKQofKQ(KQ)), &
                                              (NUNOCC), VVVV(1, 1, RNA, RNB), VVVV_R(1, 1, RNA, RNB), zero)

                           END DO
                        END DO
!$OMP END DO
!$omp end parallel

!                  CALL OMP_SET_NUM_THREADS(1)

                        VVVV_S_R = VVVV_R

                        IF (SINGLES) THEN
                           DO RNA = 1, RNBLOCKA !(NUNOCC)
                           DO NK = 1, VBMAX
                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                           FTOD_PW_AB(1, 1, NA + RNA - 1, RKQofKQ(KQ), KA, 1), FTOD_OC_AB(1, 1, NA + RNA - 1, RKQofKQ(KQ), KA, 1), &
                                              (NUNOCC), &
                                              PW_IA_TMP(1, 1, NK, RKQofKQ(KQ)), OC_IA_TMP(1, 1, NK, RKQofKQ(KQ)), &
                                              (NUNOCC), VV(1, 1), VV_R(1, 1), zero)
                              OVVV_S_R(NK, :, :, NA + RNA - 1) = -VV_R(:, :)
                           END DO
                           END DO
                           VO_S_R(:, :) = T1_R(:, :, RKIofKI(KB))
                           CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)

                           DO RNA = 1, RNBLOCKA !(NUNOCC)
                           DO RNB = 1, RNBLOCKB !(NUNOCC)

                              CALL SGEMM('t', 'n', (NUNOCC)*(NUNOCC), 1, &
                                         (VBMAX), (1._qs), OVVV_S_R(1, 1, 1, NA + RNA - 1), &
                                         (VBMAX), &
                                         OV_S_R(1, NB + RNB - 1), (VBMAX), &
                                         (1._qs), VVVV_S_R(1, 1, RNA, RNB), (NUNOCC)*(NUNOCC))

                           END DO
                           END DO

                           DO RNB = 1, RNBLOCKB !(NUNOCC)
                           DO NK = 1, VBMAX
                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                              FTOD_PW_IA(1, 1, NK, RKQofKQ(KQ), KA, 1), FTOD_OC_IA(1, 1, NK, RKQofKQ(KQ), KA, 1), &
                                              (NUNOCC), &
                                           PW_AB_TMP(1, 1, NB + RNB - 1, RKQofKQ(KQ)), OC_AB_TMP(1, 1, NB + RNB - 1, RKQofKQ(KQ)), &
                                              (NUNOCC), VV(1, 1), VV_R(1, 1), zero)
                              OVVV_S_R(NK, :, :, NB + RNB - 1) = -VV_R(:, :)
                           END DO
                           END DO
                           VO_S_R(:, :) = T1_R(:, :, KPTS_MKPTS(KA))
                           CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)

                           DO RNA = 1, RNBLOCKA !(NUNOCC)
                           DO RNB = 1, RNBLOCKB !(NUNOCC)

                              CALL SGEMM('t', 'n', (NUNOCC)*(NUNOCC), 1, &
                                         (VBMAX), (1._qs), OVVV_S_R(1, 1, 1, NB + RNB - 1), &
                                         (VBMAX), &
                                         OV_S_R(1, NA + RNA - 1), (VBMAX), &
                                         (1._qs), VVVV_S_R(1, 1, RNA, RNB), (NUNOCC)*(NUNOCC))

                           END DO
                           END DO

                        END IF

                        VVVV_S_R = VVVV_S_R*KPOINTS_FULL%WTKPT(1)

!               CALL OMP_SET_NUM_THREADS(nthreads)

!$omp parallel shared (T2_N_R,T2_TMP_R,RNBLOCKA,NUNOCC,VVVV_S_R,NA,NB,KA), private (NI,NJ,myid,NC,ND,KI,KQ,KJ,RNA,RNB,VVOO_S_R)
                        IF (ALLOCATED(VVOO_S_R)) DEALLOCATE (VVOO_S_R)
                        ALLOCATE (VVOO_S_R(NUNOCC, NUNOCC, VBMAX, VBMAX))
!$OMP DO SCHEDULE(STATIC)
                        DO KI = 1, WDES%NKPTS
                           IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE

                           KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                                    WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                           KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                                    WDES%VKPT(:, KQ), KPOINTS_FULL)
                           IF (KJ > KI) CYCLE

                           VVOO_S_R(:, :, :, :) = T2_TMP_R(:, :, :, :, RKIofKI(KI), RKIofKI(KJ))
                           IF (SINGLES) THEN
                           IF ((KI == KC) .and. (KJ == KD)) THEN
                              DO NJ = 1, VBMAX
                              DO NI = 1, VBMAX
                              DO ND = 1, (NUNOCC)
                              DO NC = 1, (NUNOCC)
                        VVOO_S_R(NC, ND, NI, NJ) = VVOO_S_R(NC, ND, NI, NJ) + T1_R(NC, NI, RKIofKI(KI))*T1_R(ND, NJ, RKIofKI(KJ))/ &
                                                            WTKPT
                              END DO
                              END DO
                              END DO
                              END DO
                           END IF
                           END IF

                           DO NJ = 1, VBMAX
                           DO NI = 1, VBMAX

                              DO RNB = 1, RNBLOCKB

                                 CALL SGEMM('t', 'n', RNBLOCKA, 1, &
                                            (NUNOCC)*(NUNOCC), (1._qs), VVVV_S_R(1, 1, 1, RNB), &
                                            (NUNOCC)*(NUNOCC), &
                                            VVOO_S_R(1, 1, NI, NJ), (NUNOCC)*(NUNOCC), &
                                            (1._qs), T2_N_R(NA, NB + RNB - 1, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA), RNBLOCKA)

                              END DO
                           END DO
                           END DO
                        END DO
!$OMP END DO
!$omp end parallel
!               CALL OMP_SET_NUM_THREADS(1)
                     END DO !NA,NBLOCK
                  END DO !NB,NBLOCK
               END DO
            END DO
         END DO

      ELSE
         DO KC_ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KC_) == 0)) CYCLE
            KC = KC_ !KPOINT_IN_FULL_GRID(WDES%VKPT(:,KC_),KPOINTS_FULL)
            CALL BCAST2ALL_T2(WDES, KC, T2_TMP, T2_TMP_R)
            DO KB = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
               call date_and_time(values=time_array1)
               ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                      time_array2(8) - time_array1(8))
               call date_and_time(values=time_array2)
               CALL BCAST2ALL_FTOD_AB(WDES, KB, 2, PW_AB_TMP, OC_AB_TMP) !it wouldn't be necessary to communicate all kq to all processes here
               IF (SINGLES) CALL BCAST2ALL_FTOD_IA(WDES, KB, 2, PW_IA_TMP, OC_IA_TMP) !it wouldn't be necessary to communicate all kq to all processes here
               DO KA = 1, MY_NKPTS

                  DO NA = 1, (NUNOCC), NBLOCKAB
                     RNBLOCKA = MIN((NUNOCC) - NA + 1, NBLOCKAB)
                     DO NB = 1, (NUNOCC), NBLOCKAB
                        RNBLOCKB = MIN((NUNOCC) - NB + 1, NBLOCKAB)
!                  CALL OMP_SET_NUM_THREADS(nthreads)

                        KQ = KPOINT_IN_FULL_GRID((WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                                  WDES%VKPT(:, KC)), KPOINTS_FULL)
                        KD = KPOINT_IN_FULL_GRID((WDES%VKPT(:, KB) + &
                                                  WDES%VKPT(:, KQ)), KPOINTS_FULL)

!$omp parallel shared (RNBLOCKB,RNBLOCKA,NGVECTOR,NHVECTOR,FTOD_PW_AB,FTOD_OC_AB,NUNOCC,PW_AB_TMP,OC_AB_TMP,VVVV,nthreads), private (RNB,RNA,ND,NC,NG)
!$OMP DO SCHEDULE(STATIC)
                        DO RNB = 1, RNBLOCKB
                           DO RNA = 1, RNBLOCKA

                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                           FTOD_PW_AB(1, 1, NA + RNA - 1, RKQofKQ(KQ), KA, 1), FTOD_OC_AB(1, 1, NA + RNA - 1, RKQofKQ(KQ), KA, 1), &
                                              (NUNOCC), &
                                           PW_AB_TMP(1, 1, NB + RNB - 1, RKQofKQ(KQ)), OC_AB_TMP(1, 1, NB + RNB - 1, RKQofKQ(KQ)), &
                                              (NUNOCC), VVVV(1, 1, RNA, RNB), VVVV_R(1, 1, RNA, RNB), zero)


                           END DO
                        END DO
!$OMP END DO
!$omp end parallel

!                  CALL OMP_SET_NUM_THREADS(1)

                        VVVV_S = VVVV

                        IF (SINGLES) THEN
                           DO RNA = 1, RNBLOCKA !(NUNOCC)
                           DO NK = 1, VBMAX
                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                           FTOD_PW_AB(1, 1, NA + RNA - 1, RKQofKQ(KQ), KA, 1), FTOD_OC_AB(1, 1, NA + RNA - 1, RKQofKQ(KQ), KA, 1), &
                                              (NUNOCC), &
                                              PW_IA_TMP(1, 1, NK, RKQofKQ(KQ)), OC_IA_TMP(1, 1, NK, RKQofKQ(KQ)), &
                                              (NUNOCC), VV(1, 1), VV_R(1, 1), zero)
                              OVVV_S(NK, :, :, NA + RNA - 1) = -VV(:, :)
                           END DO
                           END DO
                           VO_S(:, :) = T1(:, :, RKIofKI(KB))
                           CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)

                           DO RNA = 1, RNBLOCKA !(NUNOCC)
                           DO RNB = 1, RNBLOCKB !(NUNOCC)

                              CALL CGEMM('t', 'n', (NUNOCC)*(NUNOCC), 1, &
                                         (VBMAX), (1._qs, 0._qs), OVVV_S(1, 1, 1, NA + RNA - 1), &
                                         (VBMAX), &
                                         OV_S(1, NB + RNB - 1), (VBMAX), &
                                         (1._qs, 0._qs), VVVV_S(1, 1, RNA, RNB), (NUNOCC)*(NUNOCC))

                           END DO
                           END DO

                           DO RNB = 1, RNBLOCKB !(NUNOCC)
                           DO NK = 1, VBMAX
                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                              FTOD_PW_IA(1, 1, NK, RKQofKQ(KQ), KA, 1), FTOD_OC_IA(1, 1, NK, RKQofKQ(KQ), KA, 1), &
                                              (NUNOCC), &
                                           PW_AB_TMP(1, 1, NB + RNB - 1, RKQofKQ(KQ)), OC_AB_TMP(1, 1, NB + RNB - 1, RKQofKQ(KQ)), &
                                              (NUNOCC), VV(1, 1), VV_R(1, 1), zero)
                              OVVV_S(NK, :, :, NB + RNB - 1) = -VV(:, :)
                           END DO
                           END DO
                           VO_S(:, :) = T1(:, :, KPTS_MKPTS(KA))
                           CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)

                           DO RNA = 1, RNBLOCKA !(NUNOCC)
                           DO RNB = 1, RNBLOCKB !(NUNOCC)

                              CALL CGEMM('t', 'n', (NUNOCC)*(NUNOCC), 1, &
                                         (VBMAX), (1._qs, 0._qs), OVVV_S(1, 1, 1, NB + RNB - 1), &
                                         (VBMAX), &
                                         OV_S(1, NA + RNA - 1), (VBMAX), &
                                         (1._qs, 0._qs), VVVV_S(1, 1, RNA, RNB), (NUNOCC)*(NUNOCC))

                           END DO
                           END DO

                        END IF

                        VVVV_S = VVVV_S*WTKPT

!               CALL OMP_SET_NUM_THREADS(nthreads)

!$omp parallel shared (T2_N,T2_TMP,RNBLOCKA,NUNOCC,VVVV_S,NA,NB,KA), private (NI,NJ,myid,NC,ND,KI,KQ,KJ,RNA,RNB,VVOO_S)
                        IF (ALLOCATED(VVOO_S)) DEALLOCATE (VVOO_S)
                        ALLOCATE (VVOO_S(NUNOCC, NUNOCC, VBMAX, VBMAX))
!$OMP DO SCHEDULE(STATIC)
                        DO KI = 1, WDES%NKPTS
                           IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE

                           KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                                    WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                           KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                                    WDES%VKPT(:, KQ), KPOINTS_FULL)
                           IF (KJ > KI) CYCLE

                           VVOO_S(:, :, :, :) = T2_TMP(:, :, :, :, RKIofKI(KI), RKIofKI(KJ))
                           IF (SINGLES) THEN
                           IF ((KI == KC) .and. (KJ == KD)) THEN
                              DO NJ = 1, VBMAX
                              DO NI = 1, VBMAX
                              DO ND = 1, (NUNOCC)
                              DO NC = 1, (NUNOCC)
                                VVOO_S(NC, ND, NI, NJ) = VVOO_S(NC, ND, NI, NJ) + T1(NC, NI, RKIofKI(KI))*T1(ND, NJ, RKIofKI(KJ))/ &
                                                          WTKPT
                              END DO
                              END DO
                              END DO
                              END DO
                           END IF
                           END IF

                           DO NJ = 1, VBMAX
                           DO NI = 1, VBMAX

                              DO RNB = 1, RNBLOCKB

                                 CALL CGEMM('t', 'n', RNBLOCKA, 1, &
                                            (NUNOCC)*(NUNOCC), (1._qs, 0._qs), VVVV_S(1, 1, 1, RNB), &
                                            (NUNOCC)*(NUNOCC), &
                                            VVOO_S(1, 1, NI, NJ), (NUNOCC)*(NUNOCC), &
                                            (1._qs, 0._qs), T2_N(NA, NB + RNB - 1, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA), RNBLOCKA)

                              END DO
                           END DO
                           END DO
                        END DO
!$OMP END DO
!$omp end parallel
!               CALL OMP_SET_NUM_THREADS(1)
                     END DO !NA,NBLOCK
                  END DO !NB,NBLOCK
               END DO
            END DO
         END DO

      END IF

   END SUBROUTINE FORM_CONTR_ABCD_PPL

!***********************************************************************
!
! The following routine computes an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE FORM_CHI_AKIC(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER KI, KJ, KA, KK, NI, NJ, NK, NL, KQ, KQ_, KC, KD, KL, NA, &
         ND
      COMPLEX(qs) :: KW

      call BLACS_PINFO(ME, PROCS)
      IF (LORBREAL) THEN
      DO KK = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KK) == 0)) CYCLE
         CALL BCAST2ALL_FTOD_IA(WDES, KK, 2, PW_IA_TMP, OC_IA_TMP)
         DO KA = 1, MY_NKPTS
!               CALL OMP_SET_NUM_THREADS(nthreads)
!$omp parallel shared (KK,KA,PW_IA_TMP,FTOD_PW_AI,OC_IA_TMP,WDES,CHI_CKAI_R,FTOD_PW_IJ,NUNOCC,T1_R), &
!$omp private (KQ,KI,KC,OVVO_R,VOVO_S_R,NK,NL,VO_S_R,OV_S_R,OOVV_S_R,OOOV_S_R,VV_R,VVOO_S_R)
            IF (.not. ALLOCATED(OVVO_R)) ALLOCATE (OVVO_R(VBMAX, NUNOCC, NUNOCC, VBMAX))
            IF (.not. ALLOCATED(VOVO_S_R)) ALLOCATE (VOVO_S_R(NUNOCC, VBMAX, NUNOCC, VBMAX))
            IF (.not. ALLOCATED(VO_S_R)) ALLOCATE (VO_S_R(NUNOCC, VBMAX))
            IF (.not. ALLOCATED(OV_S_R)) ALLOCATE (OV_S_R(VBMAX, NUNOCC))
            IF (.not. ALLOCATED(OOVV_S_R)) ALLOCATE (OOVV_S_R(VBMAX, VBMAX, NUNOCC, NUNOCC))
            IF (.not. ALLOCATED(OOOV_S_R)) ALLOCATE (OOOV_S_R(VBMAX, VBMAX, VBMAX, NUNOCC))
            IF (.not. ALLOCATED(VV_R)) ALLOCATE (VV_R(NUNOCC, NUNOCC))
            IF (.not. ALLOCATED(VVOO_S_R)) ALLOCATE (VVOO_S_R(NUNOCC, NUNOCC, VBMAX, VBMAX))
!$OMP DO SCHEDULE(STATIC)
            DO KQ = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0)) CYCLE
               KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)
               KC = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KK) + &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)

               CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                               FTOD_PW_AI(1, 1, 1, RKQofKQ(KQ), KA, 1), FTOD_OC_AI(1, 1, 1, RKQofKQ(KQ), KA, 1), &
                               (NUNOCC)*VBMAX, &
                               PW_IA_TMP(1, 1, 1, RKQofKQ(KQ)), OC_IA_TMP(1, 1, 1, RKQofKQ(KQ)), &
                               (NUNOCC)*VBMAX, OVVO(1, 1, 1, 1), OVVO_R(1, 1, 1, 1), zero)
               CALL SORT_O1V1V2O2_V2O2V1O1(WDES, OVVO, &
                                           VOVO_S, OVVO_R, VOVO_S_R)

               CHI_CKAI_R(:, :, :, :, RKIofKI(KI), RKIofKI(KK), KA) = VOVO_S_R(:, :, :, :)

               IF (SINGLES) THEN

                  DO NK = 1, VBMAX
                  DO NL = 1, VBMAX

                     CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                     FTOD_PW_IJ(1, 1, NL, RKQofKQ(KQ), KA, 1), FTOD_OC_IJ(1, 1, NL, RKQofKQ(KQ), KA, 1), &
                                     (VBMAX), &
                                     PW_IA_TMP(1, 1, NK, RKQofKQ(KQ)), OC_IA_TMP(1, 1, NK, RKQofKQ(KQ)), &
                                     (NUNOCC), OVVO(1, 1, 1, 1), OVVO_R(1, 1, 1, 1), zero)
                     !ovvo : ic

                     OOOV_S_R(NL, NK, :, :) = -OVVO_R(:, :, 1, 1)
                     !ooov :lkic

                  END DO
                  END DO

                  VO_S_R(:, :) = T1_R(:, :, KPTS_MKPTS(KA))
                  CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)
                  CALL SGEMM('t', 'n', (NUNOCC)*(VBMAX)*(VBMAX), (NUNOCC), &
                             (VBMAX), (1._qs), OOOV_S_R(1, 1, 1, 1), (VBMAX), &
                             OV_S_R(1, 1), (VBMAX), &
                             (0._qs), OOVV_S_R(1, 1, 1, 1), (NUNOCC)*(VBMAX)*(VBMAX))
                  !oovv : kica -> ckai : o1o2v1v2 -> v1o1v2o2

                  CALL SORT_O1O2V1V2_V1O1V2O2(WDES, OOVV_S, VOVO_S, OOVV_S_R, VOVO_S_R)
  CHI_CKAI_R(:, :, :, :, RKIofKI(KI), RKIofKI(KK), KA) = CHI_CKAI_R(:, :, :, :, RKIofKI(KI), RKIofKI(KK), KA) + VOVO_S_R(:, :, :, :)

                  DO NA = 1, (NUNOCC)
                  DO NK = 1, VBMAX

                     CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                     FTOD_PW_AB(1, 1, NA, RKQofKQ(KQ), KA, 1), FTOD_OC_AB(1, 1, NA, RKQofKQ(KQ), KA, 1), &
                                     (NUNOCC), &
                                     PW_IA_TMP(1, 1, NK, RKQofKQ(KQ)), OC_IA_TMP(1, 1, NK, RKQofKQ(KQ)), &
                                     (NUNOCC), VV(1, 1), VV_R(1, 1), zero)
                     !ovvo : dc

                     VVOO_S_R(:, :, NK, 1) = VV_R(:, :)
                     !vvoo :dck

                  END DO

                  DO NI = 1, VBMAX
                     CALL SGEMM('t', 'n', (NUNOCC)*(VBMAX), 1, &
                                (NUNOCC), (1._qs), VVOO_S_R(1, 1, 1, 1), (NUNOCC), &
                                T1_R(1, NI, RKIofKI(KI)), (NUNOCC), &
                                (1._qs), CHI_CKAI_R(1, 1, NA, NI, RKIofKI(KI), RKIofKI(KK), KA), (NUNOCC)*(VBMAX))
                  END DO

                  END DO

               END IF

            END DO
!$OMP END DO
!$omp end parallel

!               CALL OMP_SET_NUM_THREADS(1)
         END DO
      END DO

      KW = WTKPT

      IF (.NOT. LCCD) THEN

         DO KK = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KK) == 0)) CYCLE
            CALL BCAST2ALL_FTOD_IA(WDES, KK, 2, PW_IA_TMP2, OC_IA_TMP2)
            DO KL = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KL) == 0)) CYCLE
               CALL BCAST2ALL_FTOD_IA(WDES, KL, 1, PW_IA_TMP, OC_IA_TMP)

               DO KQ = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0)) CYCLE

                  CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                  PW_IA_TMP(1, 1, 1, RKQofKQ(KQ)), OC_IA_TMP(1, 1, 1, RKQofKQ(KQ)), &
                                  (NUNOCC)*VBMAX, &
                                  PW_IA_TMP2(1, 1, 1, RKQofKQ(KQ)), OC_IA_TMP2(1, 1, 1, RKQofKQ(KQ)), &
                                  (NUNOCC)*VBMAX, VOVO(1, 1, 1, 1), VOVO_R(1, 1, 1, 1), zero)

                  VOVO_S_R = VOVO_R

                  IF (.NOT. RING) THEN

                     KD = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) - &
                                              WDES%VKPT(:, KQ), KPOINTS_FULL)

                     KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KK) - &
                                               WDES%VKPT(:, KD), KPOINTS_FULL)

                     DO KA = 1, MY_NKPTS

                        KD = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) - &
                                                 WDES%VKPT(:, KQ), KPOINTS_FULL)

                        KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) - &
                                                  WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                        KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KD) - &
                                                 WDES%VKPT(:, KQ_), KPOINTS_FULL)

                        VVOO_S_R(:, :, :, :) = T2_R(:, :, :, :, RKIofKI(KL), RKIofKI(KI), KA)
                        IF ((SINGLES) .and. (KL == KPTS_MKPTS(KA))) THEN
                           DO NI = 1, VBMAX
                           DO NL = 1, VBMAX
                           DO ND = 1, (NUNOCC)
                           VVOO_S_R(:, ND, NL, NI) = VVOO_S_R(:, ND, NL, NI) + T1_R(:, NL, RKIofKI(KL))*T1_R(ND, NI, RKIofKI(KI))/ &
                                                        KW*(2.0_qs)
                           END DO
                           END DO
                           END DO
                        END IF

                        CALL SORT_V1V2O1O2_V2O1V1O2(WDES, VVOO_S, VOVO2_S, VVOO_S_R, VOVO2_S_R)
                        CALL SGEMM('t', 'n', (NUNOCC)*VBMAX, &
                                   (NUNOCC)*VBMAX, (NUNOCC)*(VBMAX), &
                                   -(0.5_qs)*KW, VOVO_S_R(1, 1, 1, 1), (NUNOCC)*(VBMAX), &
                                   VOVO2_S_R(1, 1, 1, 1), (NUNOCC)*(VBMAX), &
                                   (1.0_qs), CHI_CKAI_R(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KK), KA), (NUNOCC)*(VBMAX))

                     END DO !KA

                  END IF

                  KD = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) - &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)

                  KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KD) - &
                                            WDES%VKPT(:, KK), KPOINTS_FULL)

                  IF (.NOT. RING) THEN

                     CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                     PW_IA_TMP(1, 1, 1, RKQofKQ(KQ_)), OC_IA_TMP(1, 1, 1, RKQofKQ(KQ_)), &
                                     (NUNOCC)*VBMAX, &
                                     PW_IA_TMP2(1, 1, 1, RKQofKQ(KQ_)), OC_IA_TMP2(1, 1, 1, RKQofKQ(KQ_)), &
                                     (NUNOCC)*VBMAX, VOVO(1, 1, 1, 1), VOVO_R(1, 1, 1, 1), zero)

                     VOVO2_S_R = VOVO_R

                  END IF

                  CALL SORT_V1O1V2O2_V2O1V1O2(WDES, VOVO2_S, VOVO3_S, VOVO2_S_R, VOVO3_S_R)

                  IF ((.NOT. RING) .and. (.not. LDISTING)) VOVO_S_R = (2.0_qs)*VOVO_S_R - VOVO3_S_R
                  IF ((.NOT. RING) .and. LDISTING) VOVO_S_R = (2.0_qs)*VOVO_S_R !-VOVO3_S_R
                  IF (RING) VOVO_S_R = (2.0_qs)*VOVO_S_R

                  DO KA = 1, MY_NKPTS

                     KD = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) - &
                                              WDES%VKPT(:, KQ), KPOINTS_FULL)

                     KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KD) - &
                                               WDES%VKPT(:, KL), KPOINTS_FULL)
                     KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) + &
                                              WDES%VKPT(:, KQ_), KPOINTS_FULL)

                     VVOO_S_R(:, :, :, :) = T2_R(:, :, :, :, RKIofKI(KI), RKIofKI(KL), KA)
                     CALL SORT_V1V2O1O2_V2O2V1O1(WDES, VVOO_S, VOVO2_S, VVOO_S_R, VOVO2_S_R)
                     CALL SGEMM('t', 'n', (NUNOCC)*VBMAX, &
                                (NUNOCC)*VBMAX, (NUNOCC)*(VBMAX), &
                                (0.5_qs)*KW, VOVO_S_R(1, 1, 1, 1), (NUNOCC)*(VBMAX), &
                                VOVO2_S_R(1, 1, 1, 1), (NUNOCC)*(VBMAX), &
                                (1.0_qs), CHI_CKAI_R(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KK), KA), (NUNOCC)*(VBMAX))

                  END DO !KA

               END DO !KQ
            END DO !KL
         END DO

      END IF !LCCD

      ELSE

      DO KK = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KK) == 0)) CYCLE
         CALL BCAST2ALL_FTOD_IA(WDES, KK, 2, PW_IA_TMP, OC_IA_TMP)
         DO KA = 1, MY_NKPTS
!               CALL OMP_SET_NUM_THREADS(nthreads)

!$omp parallel shared (KK,KA,PW_IA_TMP,FTOD_PW_AI,WDES,CHI_CKAI,FTOD_PW_IJ,NUNOCC,T1), &
!$omp private (KQ,KI,KC,OVVO,VOVO_S,NK,NL,VO_S,OV_S,OOVV_S,OOOV_S,VV,VVOO_S)
            IF (ALLOCATED(OVVO)) DEALLOCATE (OVVO)
            ALLOCATE (OVVO(VBMAX, NUNOCC, NUNOCC, VBMAX))
            IF (ALLOCATED(VOVO_S)) DEALLOCATE (VOVO_S)
            ALLOCATE (VOVO_S(NUNOCC, VBMAX, NUNOCC, VBMAX))
            IF (ALLOCATED(VO_S)) DEALLOCATE (VO_S)
            ALLOCATE (VO_S(NUNOCC, VBMAX))
            IF (ALLOCATED(OV_S)) DEALLOCATE (OV_S)
            ALLOCATE (OV_S(VBMAX, NUNOCC))
            IF (ALLOCATED(OOVV_S)) DEALLOCATE (OOVV_S)
            ALLOCATE (OOVV_S(VBMAX, VBMAX, NUNOCC, NUNOCC))
            IF (ALLOCATED(OOOV_S)) DEALLOCATE (OOOV_S)
            ALLOCATE (OOOV_S(VBMAX, VBMAX, VBMAX, NUNOCC))
            IF (ALLOCATED(VV)) DEALLOCATE (VV)
            ALLOCATE (VV(NUNOCC, NUNOCC))
            IF (ALLOCATED(VVOO_S)) DEALLOCATE (VVOO_S)
            ALLOCATE (VVOO_S(NUNOCC, NUNOCC, VBMAX, VBMAX))
!$OMP DO SCHEDULE(STATIC)
            DO KQ = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0)) CYCLE
               KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)
               KC = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KK) + &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)

               CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                               FTOD_PW_AI(1, 1, 1, RKQofKQ(KQ), KA, 1), FTOD_OC_AI(1, 1, 1, RKQofKQ(KQ), KA, 1), &
                               (NUNOCC)*VBMAX, &
                               PW_IA_TMP(1, 1, 1, RKQofKQ(KQ)), OC_IA_TMP(1, 1, 1, RKQofKQ(KQ)), &
                               (NUNOCC)*VBMAX, OVVO(1, 1, 1, 1), OVVO_R(1, 1, 1, 1), zero)
               CALL SORT_O1V1V2O2_V2O2V1O1(WDES, OVVO, &
                                           VOVO_S, OVVO_R, VOVO_S_R)

               CHI_CKAI(:, :, :, :, RKIofKI(KI), RKIofKI(KK), KA) = VOVO_S(:, :, :, :)

               IF (SINGLES) THEN

                  DO NK = 1, VBMAX
                  DO NL = 1, VBMAX

                     CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                     FTOD_PW_IJ(1, 1, NL, RKQofKQ(KQ), KA, 1), FTOD_OC_IJ(1, 1, NL, RKQofKQ(KQ), KA, 1), &
                                     (VBMAX), &
                                     PW_IA_TMP(1, 1, NK, RKQofKQ(KQ)), OC_IA_TMP(1, 1, NK, RKQofKQ(KQ)), &
                                     (NUNOCC), OVVO(1, 1, 1, 1), OVVO_R(1, 1, 1, 1), zero)
                     !ovvo : ic

                     OOOV_S(NL, NK, :, :) = -OVVO(:, :, 1, 1)
                     !ooov :lkic

                  END DO
                  END DO

                  VO_S(:, :) = T1(:, :, KPTS_MKPTS(KA))
                  CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)
                  CALL CGEMM('t', 'n', (NUNOCC)*(VBMAX)*(VBMAX), (NUNOCC), &
                             (VBMAX), (1._qs, 0._qs), OOOV_S(1, 1, 1, 1), (VBMAX), &
                             OV_S(1, 1), (VBMAX), &
                             (0._qs, 0._qs), OOVV_S(1, 1, 1, 1), (NUNOCC)*(VBMAX)*(VBMAX))
                  !oovv : kica -> ckai : o1o2v1v2 -> v1o1v2o2

                  CALL SORT_O1O2V1V2_V1O1V2O2(WDES, OOVV_S, VOVO_S, OOVV_S_R, VOVO_S_R)
        CHI_CKAI(:, :, :, :, RKIofKI(KI), RKIofKI(KK), KA) = CHI_CKAI(:, :, :, :, RKIofKI(KI), RKIofKI(KK), KA) + VOVO_S(:, :, :, :)

                  DO NA = 1, (NUNOCC)
                  DO NK = 1, VBMAX

                     CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                     FTOD_PW_AB(1, 1, NA, RKQofKQ(KQ), KA, 1), FTOD_OC_AB(1, 1, NA, RKQofKQ(KQ), KA, 1), &
                                     (NUNOCC), &
                                     PW_IA_TMP(1, 1, NK, RKQofKQ(KQ)), OC_IA_TMP(1, 1, NK, RKQofKQ(KQ)), &
                                     (NUNOCC), VV(1, 1), VV_R(1, 1), zero)
                     !ovvo : dc

                     VVOO_S(:, :, NK, 1) = VV(:, :)
                     !vvoo :dck

                  END DO

                  DO NI = 1, VBMAX
                     CALL CGEMM('t', 'n', (NUNOCC)*(VBMAX), 1, &
                                (NUNOCC), (1._qs, 0._qs), VVOO_S(1, 1, 1, 1), (NUNOCC), &
                                T1(1, NI, RKIofKI(KI)), (NUNOCC), &
                                (1._qs, 0._qs), CHI_CKAI(1, 1, NA, NI, RKIofKI(KI), RKIofKI(KK), KA), (NUNOCC)*(VBMAX))
                  END DO

                  END DO

               END IF

            END DO
!$OMP END DO
!$omp end parallel

!               CALL OMP_SET_NUM_THREADS(1)
         END DO
      END DO

      KW = WTKPT

      IF (.NOT. LCCD) THEN

         DO KK = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KK) == 0)) CYCLE
            CALL BCAST2ALL_FTOD_IA(WDES, KK, 2, PW_IA_TMP2, OC_IA_TMP2)
            DO KL = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KL) == 0)) CYCLE
               CALL BCAST2ALL_FTOD_IA(WDES, KL, 1, PW_IA_TMP, OC_IA_TMP)

               DO KQ = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0)) CYCLE

                  CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                  PW_IA_TMP(1, 1, 1, RKQofKQ(KQ)), OC_IA_TMP(1, 1, 1, RKQofKQ(KQ)), &
                                  (NUNOCC)*VBMAX, &
                                  PW_IA_TMP2(1, 1, 1, RKQofKQ(KQ)), OC_IA_TMP2(1, 1, 1, RKQofKQ(KQ)), &
                                  (NUNOCC)*VBMAX, VOVO(1, 1, 1, 1), VOVO_R(1, 1, 1, 1), zero)

                  VOVO_S = VOVO

                  IF (.NOT. RING) THEN

                     KD = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) - &
                                              WDES%VKPT(:, KQ), KPOINTS_FULL)

                     KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KK) - &
                                               WDES%VKPT(:, KD), KPOINTS_FULL)

                     DO KA = 1, MY_NKPTS

                        KD = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) - &
                                                 WDES%VKPT(:, KQ), KPOINTS_FULL)

                        KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) - &
                                                  WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                        KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KD) - &
                                                 WDES%VKPT(:, KQ_), KPOINTS_FULL)

                        VVOO_S(:, :, :, :) = T2(:, :, :, :, RKIofKI(KL), RKIofKI(KI), KA)
                        IF ((SINGLES) .and. (KL == KPTS_MKPTS(KA))) THEN
                           DO NI = 1, VBMAX
                           DO NL = 1, VBMAX
                           DO ND = 1, (NUNOCC)
                              VVOO_S(:, ND, NL, NI) = VVOO_S(:, ND, NL, NI) + T1(:, NL, RKIofKI(KL))*T1(ND, NI, RKIofKI(KI))/ &
                                                      KW*(2.0_qs, 0.0_qs)
                           END DO
                           END DO
                           END DO
                        END IF

                        CALL SORT_V1V2O1O2_V2O1V1O2(WDES, VVOO_S, VOVO2_S, VVOO_S_R, VOVO2_S_R)
                        CALL CGEMM('t', 'n', (NUNOCC)*VBMAX, &
                                   (NUNOCC)*VBMAX, (NUNOCC)*(VBMAX), &
                                   -(0.5_qs, 0.0_qs)*KW, VOVO_S(1, 1, 1, 1), (NUNOCC)*(VBMAX), &
                                   VOVO2_S(1, 1, 1, 1), (NUNOCC)*(VBMAX), &
                                   (1.0_qs, 0.0_qs), CHI_CKAI(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KK), KA), (NUNOCC)*(VBMAX))
                     END DO !KA

                  END IF

                  KD = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) - &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)

                  KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KD) - &
                                            WDES%VKPT(:, KK), KPOINTS_FULL)

                  IF (.NOT. RING) THEN

                     CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                     PW_IA_TMP(1, 1, 1, RKQofKQ(KQ_)), OC_IA_TMP(1, 1, 1, RKQofKQ(KQ_)), &
                                     (NUNOCC)*VBMAX, &
                                     PW_IA_TMP2(1, 1, 1, RKQofKQ(KQ_)), OC_IA_TMP2(1, 1, 1, RKQofKQ(KQ_)), &
                                     (NUNOCC)*VBMAX, VOVO(1, 1, 1, 1), VOVO_R(1, 1, 1, 1), zero)

                     VOVO2_S = VOVO

                  END IF

                  CALL SORT_V1O1V2O2_V2O1V1O2(WDES, VOVO2_S, VOVO3_S, VOVO2_S_R, VOVO3_S_R)

                  IF ((.NOT. RING) .and. (.not. LDISTING)) VOVO_S = (2.0_qs, 0.0_qs)*VOVO_S - VOVO3_S
                  IF ((.NOT. RING) .and. LDISTING) VOVO_S = (2.0_qs, 0.0_qs)*VOVO_S !-VOVO3_S
                  IF (RING) VOVO_S = (2.0_qs, 0.0_qs)*VOVO_S

                  DO KA = 1, MY_NKPTS

                     KD = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) - &
                                              WDES%VKPT(:, KQ), KPOINTS_FULL)

                     KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KD) - &
                                               WDES%VKPT(:, KL), KPOINTS_FULL)
                     KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) + &
                                              WDES%VKPT(:, KQ_), KPOINTS_FULL)

                     VVOO_S(:, :, :, :) = T2(:, :, :, :, RKIofKI(KI), RKIofKI(KL), KA)
                     CALL SORT_V1V2O1O2_V2O2V1O1(WDES, VVOO_S, VOVO2_S, VVOO_S_R, VOVO2_S_R)
                     CALL CGEMM('t', 'n', (NUNOCC)*VBMAX, &
                                (NUNOCC)*VBMAX, (NUNOCC)*(VBMAX), &
                                (0.5_qs, 0.0_qs)*KW, VOVO_S(1, 1, 1, 1), (NUNOCC)*(VBMAX), &
                                VOVO2_S(1, 1, 1, 1), (NUNOCC)*(VBMAX), &
                                (1.0_qs, 0.0_qs), CHI_CKAI(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KK), KA), (NUNOCC)*(VBMAX))

                  END DO !KA

               END DO !KQ
            END DO !KL
         END DO

      END IF !LCCD

      END IF

   END SUBROUTINE FORM_CHI_AKIC

!***********************************************************************
!
! The following routine computes an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE FORM_CHI_AKCI(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER KI, KJ, KA, KK, NI, NJ, NK, NL, KQ, KQ_, KC, KD, KL, NA, ND
      COMPLEX(qs) :: KW

      call BLACS_PINFO(ME, PROCS)

      IF (LORBREAL) THEN

         DO KK = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KK) == 0)) CYCLE
            CALL BCAST2ALL_FTOD_IJ(WDES, KK, 2, PW_IJ_TMP, OC_IJ_TMP)
            IF (SINGLES) CALL BCAST2ALL_FTOD_IA(WDES, KK, 2, PW_IA_TMP, OC_IA_TMP)
            DO KA = 1, MY_NKPTS

!               CALL OMP_SET_NUM_THREADS(nthreads)

!$omp parallel shared (KK,KA,PW_IA_TMP,FTOD_PW_AI,WDES,CHI_CKAI_R,PW_IJ_TMP,NUNOCC,T1_R), &
!$omp private (KQ,KI,KC,VVOO_R,NK,NL,VO_S_R,OV_S_R,OVOV_R,OVVO_S_R,OOVV_S_R,VOVO_S_R,VV_R,VVOO_S_R)
               IF (ALLOCATED(VVOO_R)) DEALLOCATE (VVOO_R)
               ALLOCATE (VVOO_R(NUNOCC, NUNOCC, VBMAX, VBMAX))
               IF (ALLOCATED(VOVO_S_R)) DEALLOCATE (VOVO_S_R)
               ALLOCATE (VOVO_S_R(NUNOCC, VBMAX, NUNOCC, VBMAX))
               IF (ALLOCATED(VO_S_R)) DEALLOCATE (VO_S_R)
               ALLOCATE (VO_S_R(NUNOCC, VBMAX))
               IF (ALLOCATED(OV_S_R)) DEALLOCATE (OV_S_R)
               ALLOCATE (OV_S_R(VBMAX, NUNOCC))
               IF (ALLOCATED(OVOV_R)) DEALLOCATE (OVOV_R)
               ALLOCATE (OVOV_R(VBMAX, NUNOCC, VBMAX, NUNOCC))
               IF (ALLOCATED(OVVO_S_R)) DEALLOCATE (OVVO_S_R)
               ALLOCATE (OVVO_S_R(VBMAX, NUNOCC, NUNOCC, VBMAX))
               IF (ALLOCATED(OOVV_S_R)) DEALLOCATE (OOVV_S_R)
               ALLOCATE (OOVV_S_R(VBMAX, VBMAX, NUNOCC, NUNOCC))
               IF (ALLOCATED(VOVO_S_R)) DEALLOCATE (VOVO_S_R)
               ALLOCATE (VOVO_S_R(NUNOCC, VBMAX, NUNOCC, VBMAX))
               IF (ALLOCATED(VV_R)) DEALLOCATE (VV_R)
               ALLOCATE (VV_R(NUNOCC, NUNOCC))
               IF (ALLOCATED(VVOO_S_R)) DEALLOCATE (VVOO_S_R)
               ALLOCATE (VVOO_S_R(NUNOCC, NUNOCC, VBMAX, VBMAX))
!$OMP DO SCHEDULE(STATIC)
               DO KQ = 1, WDES%NKPTS
                  KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KK) + &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)
                  KC = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)

                  IF (EH_SCREENING) THEN
                     CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                     FTOD_PW_S_AB(1, 1, 1, KQ, KA, 1), FTOD_OC_AB(1, 1, 1, KQ, KA, 1), &
                                     (NUNOCC)*(NUNOCC), &
                                     PW_IJ_TMP(1, 1, 1, KQ), OC_IJ_TMP(1, 1, 1, KQ), &
                                     VBMAX*VBMAX, VVOO(1, 1, 1, 1), VVOO_R(1, 1, 1, 1), zero)
                  ELSE
                     CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                     FTOD_PW_AB(1, 1, 1, KQ, KA, 1), FTOD_OC_AB(1, 1, 1, KQ, KA, 1), &
                                     (NUNOCC)*(NUNOCC), &
                                     PW_IJ_TMP(1, 1, 1, KQ), OC_IJ_TMP(1, 1, 1, KQ), &
                                     VBMAX*VBMAX, VVOO(1, 1, 1, 1), VVOO_R(1, 1, 1, 1), zero)
                  END IF

                  CALL SORT_V2V1O1O2_V2O2V1O1(WDES, VVOO, &
                                              VOVO_S, VVOO_R, VOVO_S_R)
                  CHI_CKAI_R(:, :, :, :, KI, KK, KA) = VOVO_S_R(:, :, :, :)

                  IF (SINGLES) THEN

                     VO_S_R(:, :) = T1_R(:, :, KPTS_MKPTS(KA))
                     CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)
                     DO NK = 1, VBMAX
                     DO NL = 1, VBMAX

                        CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                        PW_IJ_TMP(1, 1, NK, KQ), OC_IJ_TMP(1, 1, NK, KQ), &
                                        (VBMAX), &
                                        FTOD_PW_IA(1, 1, NL, KQ, KA, 1), FTOD_OC_IA(1, 1, NL, KQ, KA, 1), &
                                        (NUNOCC), OVOV(1, 1, 1, 1), OVOV_R(1, 1, 1, 1), zero)
                        !vovo : ic

                        OOVV_S_R(NL, :, :, 1) = -(OVOV_R(:, :, 1, 1))
                        !ooov :lic

                     END DO

                     CALL SGEMM('t', 'n', (NUNOCC)*(VBMAX), (NUNOCC), &
                                (VBMAX), (1._qs), OOVV_S_R(1, 1, 1, 1), (VBMAX), &
                                OV_S_R(1, 1), (VBMAX), &
                                (0._qs), OVVO_S_R(1, 1, 1, NK), (NUNOCC)*(VBMAX))
                     !ovvo : icak -> ckai : o2v1v2o1 -> v1o1v2o2

                     END DO

                     CALL SORT_O2V1V2O1_V1O1V2O2(WDES, OVVO_S, VOVO_S, OVVO_S_R, VOVO_S_R)
                     CHI_CKAI_R(:, :, :, :, KI, KK, KA) = CHI_CKAI_R(:, :, :, :, KI, KK, KA) + VOVO_S_R(:, :, :, :)

                     DO NA = 1, (NUNOCC)
                     DO NK = 1, VBMAX

                        CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                        PW_IA_TMP(1, 1, NK, KQ), OC_IA_TMP(1, 1, NK, KQ), &
                                        (NUNOCC), &
                                        FTOD_PW_AB(1, 1, NA, KQ, KA, 1), FTOD_OC_AB(1, 1, NA, KQ, KA, 1), &
                                        (NUNOCC), VV(1, 1), VV_R(1, 1), zero)
                        !vv : dc

                        VVOO_S_R(:, :, NK, 1) = (VV_R(:, :))
                        !vvoo :dck

                     END DO

                     DO NI = 1, VBMAX
                        CALL SGEMM('t', 'n', (NUNOCC)*(VBMAX), 1, &
                                   (NUNOCC), (1._qs), VVOO_S_R(1, 1, 1, 1), (NUNOCC), &
                                   T1_R(1, NI, KI), (NUNOCC), &
                                   (1._qs), CHI_CKAI_R(1, 1, NA, NI, KI, KK, KA), (NUNOCC)*(VBMAX))
                     END DO

                     END DO

                  END IF

               END DO
!$OMP END DO
!$omp end parallel

!               CALL OMP_SET_NUM_THREADS(1)

            END DO
         END DO

         KW = KPOINTS_FULL%WTKPT(1)

         IF (.NOT. LCCD) THEN

            DO KK = 1, WDES%NKPTS
               CALL BCAST2ALL_FTOD_IA(WDES, KK, 2, PW_IA_TMP2, OC_IA_TMP2)
               DO KL = 1, WDES%NKPTS
                  CALL BCAST2ALL_FTOD_IA(WDES, KL, 1, PW_IA_TMP, OC_IA_TMP)

                  DO KQ = 1, WDES%NKPTS

                     CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                     PW_IA_TMP(1, 1, 1, KQ), OC_IA_TMP(1, 1, 1, KQ), &
                                     (NUNOCC)*VBMAX, &
                                     PW_IA_TMP2(1, 1, 1, KQ), OC_IA_TMP2(1, 1, 1, KQ), &
                                     (NUNOCC)*VBMAX, VOVO(1, 1, 1, 1), VOVO_R(1, 1, 1, 1), zero)

                     VOVO2_S_R = VOVO_R
                     CALL SORT_V1O1V2O2_V2O1V1O2(WDES, VOVO2_S, VOVO_S, VOVO2_S_R, VOVO_S_R)

                     KD = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KK) + &
                                              WDES%VKPT(:, KQ), KPOINTS_FULL)

                     DO KA = 1, MY_NKPTS

                        KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) - &
                                                  WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                        KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KD) - &
                                                 WDES%VKPT(:, KQ_), KPOINTS_FULL)

                        VVOO_S_R(:, :, :, :) = T2_R(:, :, :, :, KL, KI, KA)
                        IF ((SINGLES) .and. (KL == KPTS_MKPTS(KA))) THEN
                           DO NI = 1, VBMAX
                           DO NL = 1, VBMAX
                           DO ND = 1, (NUNOCC)
                              VVOO_S_R(:, ND, NL, NI) = VVOO_S_R(:, ND, NL, NI) + T1_R(:, NL, KL)*T1_R(ND, NI, KI)/ &
                                                        KW*(2.0_qs)
                           END DO
                           END DO
                           END DO
                        END IF

                        CALL SORT_V1V2O1O2_V2O1V1O2(WDES, VVOO_S, VOVO2_S, VVOO_S_R, VOVO2_S_R)
                        IF (.not. LDISTING) THEN
                           CALL SGEMM('t', 'n', (NUNOCC)*VBMAX, &
                                      (NUNOCC)*VBMAX, (NUNOCC)*(VBMAX), &
                                      -(0.5_qs)*KW, VOVO_S_R(1, 1, 1, 1), (NUNOCC)*(VBMAX), &
                                      VOVO2_S_R(1, 1, 1, 1), (NUNOCC)*(VBMAX), &
                                      (1.0_qs), CHI_CKAI_R(1, 1, 1, 1, KI, KK, KA), (NUNOCC)*(VBMAX))
                        END IF

                     END DO !KA

                  END DO !KQ
               END DO !KL
            END DO

         END IF !LCCD

      ELSE

         DO KK = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KK) == 0)) CYCLE
            CALL BCAST2ALL_FTOD_IJ(WDES, KK, 2, PW_IJ_TMP, OC_IJ_TMP)
            IF (SINGLES) CALL BCAST2ALL_FTOD_IA(WDES, KK, 2, PW_IA_TMP, OC_IA_TMP)
            DO KA = 1, MY_NKPTS

!               CALL OMP_SET_NUM_THREADS(nthreads)

!$omp parallel shared (KK,KA,PW_IA_TMP,FTOD_PW_AI,WDES,CHI_CKAI,PW_IJ_TMP,NUNOCC,T1), &
!$omp private (KQ,KI,KC,VVOO,NK,NL,VO_S,OV_S,OVOV,OVVO_S,OOVV_S,VOVO_S,VV,VVOO_S)
               IF (ALLOCATED(VVOO)) DEALLOCATE (VVOO)
               ALLOCATE (VVOO(NUNOCC, NUNOCC, VBMAX, VBMAX))
               IF (ALLOCATED(VOVO_S)) DEALLOCATE (VOVO_S)
               ALLOCATE (VOVO_S(NUNOCC, VBMAX, NUNOCC, VBMAX))
               IF (ALLOCATED(VO_S)) DEALLOCATE (VO_S)
               ALLOCATE (VO_S(NUNOCC, VBMAX))
               IF (ALLOCATED(OV_S)) DEALLOCATE (OV_S)
               ALLOCATE (OV_S(VBMAX, NUNOCC))
               IF (ALLOCATED(OVOV)) DEALLOCATE (OVOV)
               ALLOCATE (OVOV(VBMAX, NUNOCC, VBMAX, NUNOCC))
               IF (ALLOCATED(OVVO_S)) DEALLOCATE (OVVO_S)
               ALLOCATE (OVVO_S(VBMAX, NUNOCC, NUNOCC, VBMAX))
               IF (ALLOCATED(OOVV_S)) DEALLOCATE (OOVV_S)
               ALLOCATE (OOVV_S(VBMAX, VBMAX, NUNOCC, NUNOCC))
               IF (ALLOCATED(VOVO_S)) DEALLOCATE (VOVO_S)
               ALLOCATE (VOVO_S(NUNOCC, VBMAX, NUNOCC, VBMAX))
               IF (ALLOCATED(VV)) DEALLOCATE (VV)
               ALLOCATE (VV(NUNOCC, NUNOCC))
               IF (ALLOCATED(VVOO_S)) DEALLOCATE (VVOO_S)
               ALLOCATE (VVOO_S(NUNOCC, NUNOCC, VBMAX, VBMAX))
!$OMP DO SCHEDULE(STATIC)
               DO KQ = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0)) CYCLE
                  KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KK) + &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)
                  KC = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)

                  IF (EH_SCREENING) THEN
                     CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                     FTOD_PW_S_AB(1, 1, 1, RKQofKQ(KQ), KA, 1), FTOD_OC_AB(1, 1, 1, RKQofKQ(KQ), KA, 1), &
                                     (NUNOCC)*(NUNOCC), &
                                     PW_IJ_TMP(1, 1, 1, RKQofKQ(KQ)), OC_IJ_TMP(1, 1, 1, RKQofKQ(KQ)), &
                                     VBMAX*VBMAX, VVOO(1, 1, 1, 1), VVOO_R(1, 1, 1, 1), zero)
                  ELSE
                     CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                     FTOD_PW_AB(1, 1, 1, RKQofKQ(KQ), KA, 1), FTOD_OC_AB(1, 1, 1, RKQofKQ(KQ), KA, 1), &
                                     (NUNOCC)*(NUNOCC), &
                                     PW_IJ_TMP(1, 1, 1, RKQofKQ(KQ)), OC_IJ_TMP(1, 1, 1, RKQofKQ(KQ)), &
                                     VBMAX*VBMAX, VVOO(1, 1, 1, 1), VVOO_R(1, 1, 1, 1), zero)
                  END IF

                  CALL SORT_V2V1O1O2_V2O2V1O1(WDES, VVOO, &
                                              VOVO_S, VVOO_R, VOVO_S_R)
                  CHI_CKAI(:, :, :, :, RKIofKI(KI), RKIofKI(KK), KA) = VOVO_S(:, :, :, :)

                  IF (SINGLES) THEN

                     VO_S(:, :) = T1(:, :, KPTS_MKPTS(KA))
                     CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)
                     DO NK = 1, VBMAX
                     DO NL = 1, VBMAX

                        CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                        PW_IJ_TMP(1, 1, NK, RKQofKQ(KQ)), OC_IJ_TMP(1, 1, NK, RKQofKQ(KQ)), &
                                        (VBMAX), &
                                        FTOD_PW_IA(1, 1, NL, RKQofKQ(KQ), KA, 1), FTOD_OC_IA(1, 1, NL, RKQofKQ(KQ), KA, 1), &
                                        (NUNOCC), OVOV(1, 1, 1, 1), OVOV_R(1, 1, 1, 1), zero)
                        !vovo : ic

                        OOVV_S(NL, :, :, 1) = -CONJG(OVOV(:, :, 1, 1))
                        !ooov :lic

                     END DO

                     CALL CGEMM('t', 'n', (NUNOCC)*(VBMAX), (NUNOCC), &
                                (VBMAX), (1._qs, 0._qs), OOVV_S(1, 1, 1, 1), (VBMAX), &
                                OV_S(1, 1), (VBMAX), &
                                (0._qs, 0._qs), OVVO_S(1, 1, 1, NK), (NUNOCC)*(VBMAX))
                     !ovvo : icak -> ckai : o2v1v2o1 -> v1o1v2o2

                     END DO

                     CALL SORT_O2V1V2O1_V1O1V2O2(WDES, OVVO_S, VOVO_S, OVVO_S_R, VOVO_S_R)
        CHI_CKAI(:, :, :, :, RKIofKI(KI), RKIofKI(KK), KA) = CHI_CKAI(:, :, :, :, RKIofKI(KI), RKIofKI(KK), KA) + VOVO_S(:, :, :, :)

                     DO NA = 1, (NUNOCC)
                     DO NK = 1, VBMAX

                        CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                        PW_IA_TMP(1, 1, NK, RKQofKQ(KQ)), OC_IA_TMP(1, 1, NK, RKQofKQ(KQ)), &
                                        (NUNOCC), &
                                        FTOD_PW_AB(1, 1, NA, RKQofKQ(KQ), KA, 1), FTOD_OC_AB(1, 1, NA, RKQofKQ(KQ), KA, 1), &
                                        (NUNOCC), VV(1, 1), VV_R(1, 1), zero)
                        !vv : dc

                        VVOO_S(:, :, NK, 1) = CONJG(VV(:, :))
                        !vvoo :dck

                     END DO

                     DO NI = 1, VBMAX
                        CALL CGEMM('t', 'n', (NUNOCC)*(VBMAX), 1, &
                                   (NUNOCC), (1._qs, 0._qs), VVOO_S(1, 1, 1, 1), (NUNOCC), &
                                   T1(1, NI, RKIofKI(KI)), (NUNOCC), &
                                   (1._qs, 0._qs), CHI_CKAI(1, 1, NA, NI, RKIofKI(KI), RKIofKI(KK), KA), (NUNOCC)*(VBMAX))
                     END DO

                     END DO

                  END IF

               END DO
!$OMP END DO
!$omp end parallel

!               CALL OMP_SET_NUM_THREADS(1)

            END DO
         END DO

         KW = WTKPT

         IF (.NOT. LCCD) THEN

            DO KK = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KK) == 0)) CYCLE
               CALL BCAST2ALL_FTOD_IA(WDES, KK, 2, PW_IA_TMP2, OC_IA_TMP2)
               DO KL = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KL) == 0)) CYCLE
                  CALL BCAST2ALL_FTOD_IA(WDES, KL, 1, PW_IA_TMP, OC_IA_TMP)

                  DO KQ = 1, WDES%NKPTS
                     IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0)) CYCLE

                     CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                     PW_IA_TMP(1, 1, 1, RKQofKQ(KQ)), OC_IA_TMP(1, 1, 1, RKQofKQ(KQ)), &
                                     (NUNOCC)*VBMAX, &
                                     PW_IA_TMP2(1, 1, 1, RKQofKQ(KQ)), OC_IA_TMP2(1, 1, 1, RKQofKQ(KQ)), &
                                     (NUNOCC)*VBMAX, VOVO(1, 1, 1, 1), VOVO_R(1, 1, 1, 1), zero)

                     VOVO2_S = VOVO
                     CALL SORT_V1O1V2O2_V2O1V1O2(WDES, VOVO2_S, VOVO_S, VOVO2_S_R, VOVO_S_R)

                     KD = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KK) + &
                                              WDES%VKPT(:, KQ), KPOINTS_FULL)

                     DO KA = 1, MY_NKPTS

                        KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KL) - &
                                                  WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                        KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KD) - &
                                                 WDES%VKPT(:, KQ_), KPOINTS_FULL)

                        VVOO_S(:, :, :, :) = T2(:, :, :, :, RKIofKI(KL), RKIofKI(KI), KA)
                        IF ((SINGLES) .and. (KL == KPTS_MKPTS(KA))) THEN
                           DO NI = 1, VBMAX
                           DO NL = 1, VBMAX
                           DO ND = 1, (NUNOCC)
                              VVOO_S(:, ND, NL, NI) = VVOO_S(:, ND, NL, NI) + T1(:, NL, RKIofKI(KL))*T1(ND, NI, RKIofKI(KI))/ &
                                                      KW*(2.0_qs, 0.0_qs)
                           END DO
                           END DO
                           END DO
                        END IF

                        CALL SORT_V1V2O1O2_V2O1V1O2(WDES, VVOO_S, VOVO2_S, VVOO_S_R, VOVO2_S_R)
                        IF (.not. LDISTING) THEN
                           CALL CGEMM('t', 'n', (NUNOCC)*VBMAX, &
                                      (NUNOCC)*VBMAX, (NUNOCC)*(VBMAX), &
                                      -(0.5_qs, 0.0_qs)*KW, VOVO_S(1, 1, 1, 1), (NUNOCC)*(VBMAX), &
                                      VOVO2_S(1, 1, 1, 1), (NUNOCC)*(VBMAX), &
                                      (1.0_qs, 0.0_qs), CHI_CKAI(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KK), KA), (NUNOCC)*(VBMAX))
                        END IF
                     END DO !KA

                  END DO !KQ
               END DO !KL
            END DO

         END IF !LCCD

      END IF  !LORBREAL

   END SUBROUTINE FORM_CHI_AKCI

!***********************************************************************
!
! The following routine contracts an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE CONTR_CHI_AKIC_AS_T2_BCKJ(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
         KC
      COMPLEX(qs) :: KW
      LOGICAL :: TMPLDISTING

      KW = WTKPT*(1.0_qs, 0.0_qs)

      IF (LORBREAL) THEN

         DO KB = 1, WDES%NKPTS
         DO KJ = 1, WDES%NKPTS
            !calculate the antisymmetrized T2 already before it is being sent
            CALL BCAST2ALL_AS_T2_KJ_KA(WDES, KJ, KB, T2_KTMP, T2_KTMP_R)

            CALL SORT_V1V2O1O2_V2O1V1O2_KTMP(WDES, T2_KTMP, T2_KTMP2, T2_KTMP_R, T2_KTMP2_R)

            KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                     WDES%VKPT(:, KJ), KPOINTS_FULL)
            DO KA = 1, MY_NKPTS
               VOVO3_S_R = (0.0_qs)
               KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) + &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)

               DO KK = 1, WDES%NKPTS
                  CALL SGEMM('t', 'n', ((NUNOCC)*(VBMAX)), ((NUNOCC)*VBMAX), &
                             ((NUNOCC)*VBMAX), (1.0_qs), CHI_CKAI_R(1, 1, 1, 1, KI, KK, KA), ((NUNOCC)*VBMAX), &
                             T2_KTMP2_R(1, 1, 1, 1, KK), ((NUNOCC)*VBMAX), &
                             (1.0_qs), VOVO3_S_R(1, 1, 1, 1), ((NUNOCC)*VBMAX))
               END DO

               CALL SORT_V1O1V2O2_V1V2O1O2(WDES, VOVO3_S, VVOO_S, VOVO3_S_R, VVOO_S_R)
               CALL SAXPY((NUNOCC)*(VBMAX)*(NUNOCC)*(VBMAX), REAL(KW, KIND=qs), &
                          VVOO_S_R(1, 1, 1, 1), 1, T2_N_R(1, 1, 1, 1, KI, KJ, KA), 1)

            END DO
         END DO
         END DO !KB

      ELSE
         DO KB = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
            DO KJ = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE
               !calculate the antisymmetrized T2 already before it is being sent
               CALL BCAST2ALL_AS_T2_KJ_KA(WDES, KJ, KB, T2_KTMP, T2_KTMP_R)

               CALL SORT_V1V2O1O2_V2O1V1O2_KTMP(WDES, T2_KTMP, T2_KTMP2, T2_KTMP_R, T2_KTMP2_R)

               KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                        WDES%VKPT(:, KJ), KPOINTS_FULL)
               DO KA = 1, MY_NKPTS
                  VOVO3_S = (0.0_qs, 0.0_qs)
                  KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) + &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)

                  DO KK = 1, WDES%NKPTS
                     IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KK) == 0)) CYCLE
                     CALL CGEMM('t', 'n', ((NUNOCC)*(VBMAX)), ((NUNOCC)*VBMAX), &
                         ((NUNOCC)*VBMAX), (1.0_qs, 0.0_qs), CHI_CKAI(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KK), KA), ((NUNOCC)*VBMAX), &
                                T2_KTMP2(1, 1, 1, 1, RKIofKI(KK)), ((NUNOCC)*VBMAX), &
                                (1.0_qs, 0.0_qs), VOVO3_S(1, 1, 1, 1), ((NUNOCC)*VBMAX))
                  END DO

                  CALL SORT_V1O1V2O2_V1V2O1O2(WDES, VOVO3_S, VVOO_S, VOVO3_S_R, VVOO_S_R)
                  CALL CAXPY((NUNOCC)*(VBMAX)*(NUNOCC)*(VBMAX), KW, &
                             VVOO_S(1, 1, 1, 1), 1, T2_N(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KJ), KA), 1)

               END DO
            END DO
         END DO !KB

      END IF

   END SUBROUTINE CONTR_CHI_AKIC_AS_T2_BCKJ

!***********************************************************************
!
! The following routine contracts an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE CONTR_CHI_AKCI_T2_BCKJ(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
         KC
      COMPLEX(qs) :: KW

      KW = WTKPT*(1.0_qs, 0.0_qs)

      IF (LORBREAL) THEN

         DO KB = 1, WDES%NKPTS
         DO KJ = 1, WDES%NKPTS
            T2_KTMP2_R(:, :, :, :, :) = (0.0_qs)

            DO KK = 1, WDES%NKPTS
               CALL BCAST2ALL_T2_KI_KJ_KA(WDES, KJ, KK, KB, VVOO_S, VVOO_S_R)
               CALL SORT_V1V2O1O2_V2O2V1O1(WDES, VVOO_S, VOVO_S, VVOO_S_R, VOVO_S_R)
               KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                        WDES%VKPT(:, KJ), KPOINTS_FULL)
               DO KA = 1, MY_NKPTS
                  KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) + &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)
                  CALL SGEMM('t', 'n', ((NUNOCC)*(VBMAX)), ((NUNOCC)*VBMAX), &
                             ((NUNOCC)*VBMAX), (1.0_qs), CHI_CKAI_R(1, 1, 1, 1, KI, KK, KA), ((NUNOCC)*VBMAX), &
                             VOVO_S_R(1, 1, 1, 1), ((NUNOCC)*VBMAX), &
                             (1.0_qs), T2_KTMP2_R(1, 1, 1, 1, KA), ((NUNOCC)*VBMAX))
               END DO
            END DO
            KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                     WDES%VKPT(:, KJ), KPOINTS_FULL)
            DO KA = 1, MY_NKPTS
               KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) + &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)
               VOVO3_S_R(:, :, :, :) = T2_KTMP2_R(:, :, :, :, KA)
               CALL SORT_V1O1V2O2_V1V2O1O2(WDES, VOVO3_S, VVOO_S, VOVO3_S_R, VVOO_S_R)
               CALL SAXPY((NUNOCC)*(VBMAX)*(NUNOCC)*(VBMAX), -KW, &
                          VVOO_S_R(1, 1, 1, 1), 1, T2_N_R(1, 1, 1, 1, KI, KJ, KA), 1)

               IF (PRECONDITION) THEN
                  DO NA = 1, NUNOCC
                  DO NB = 1, NUNOCC
                     DO NI = 1, VBMAX
                     DO NJ = 1, VBMAX
              PRET2(NA, NB, NI, NJ, KI, KJ, KA) = PRET2(NA, NB, NI, NJ, KI, KJ, KA) - CONJG(CHI_CKAI(NB, NJ, NA, NI, KI, KJ, KA))*KW
                     END DO
                     END DO
                  END DO
                  END DO
               END IF
            END DO

         END DO
         END DO !KB

      ELSE !LORBREAL

         DO KB = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
            DO KJ = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE
               T2_KTMP2(:, :, :, :, :) = (0.0_qs, 0.0_qs)

               DO KK = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KK) == 0)) CYCLE
                  CALL BCAST2ALL_T2_KI_KJ_KA(WDES, KJ, KK, KB, VVOO_S, VVOO_S_R)
                  CALL SORT_V1V2O1O2_V2O2V1O1(WDES, VVOO_S, VOVO_S, VVOO_S_R, VOVO_S_R)
                  KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                           WDES%VKPT(:, KJ), KPOINTS_FULL)
                  DO KA = 1, MY_NKPTS
                     KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) + &
                                              WDES%VKPT(:, KQ), KPOINTS_FULL)
                     CALL CGEMM('t', 'n', ((NUNOCC)*(VBMAX)), ((NUNOCC)*VBMAX), &
                         ((NUNOCC)*VBMAX), (1.0_qs, 0.0_qs), CHI_CKAI(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KK), KA), ((NUNOCC)*VBMAX), &
                                VOVO_S(1, 1, 1, 1), ((NUNOCC)*VBMAX), &
                                (1.0_qs, 0.0_qs), T2_KTMP2(1, 1, 1, 1, KA), ((NUNOCC)*VBMAX))
                  END DO
               END DO
               KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                        WDES%VKPT(:, KJ), KPOINTS_FULL)
               DO KA = 1, MY_NKPTS
                  KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) + &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)
                  VOVO3_S(:, :, :, :) = T2_KTMP2(:, :, :, :, KA)
                  CALL SORT_V1O1V2O2_V1V2O1O2(WDES, VOVO3_S, VVOO_S, VOVO3_S_R, VVOO_S_R)
                  CALL CAXPY((NUNOCC)*(VBMAX)*(NUNOCC)*(VBMAX), -KW, &
                             VVOO_S(1, 1, 1, 1), 1, T2_N(1, 1, 1, 1, RKIofKI(KI), RKIofKI(KJ), KA), 1)

                  IF (PRECONDITION) THEN
                     DO NA = 1, NUNOCC
                     DO NB = 1, NUNOCC
                        DO NI = 1, VBMAX
                        DO NJ = 1, VBMAX
                           PRET2(NA,NB,NI,NJ,RKIofKI(KI),RKIofKI(KJ),KA)=PRET2(NA,NB,NI,NJ,RKIofKI(KI),RKIofKI(KJ),KA)-CONJG(CHI_CKAI(NB,NJ,NA,NI,RKIofKI(KI),RKIofKI(KJ),KA))*KW
                        END DO
                        END DO
                     END DO
                     END DO
                  END IF
               END DO

            END DO
         END DO !KB

      END IF

   END SUBROUTINE CONTR_CHI_AKCI_T2_BCKJ

!***********************************************************************
!
! The following routine contracts an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE CONTR_CHI_BKCI_T2_ACKJ(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
         KC
      COMPLEX(qs) :: KW

      KW = WTKPT*(1.0_qs)

      IF (LORBREAL) THEN

         DO KB = 1, WDES%NKPTS
         DO KI = 1, WDES%NKPTS
            T2_KTMP2_R(:, :, :, :, :) = (0.0_qs)
            DO KK = 1, WDES%NKPTS
               CALL BCAST2ALL_CHI_CKAI_KJ_KK_KA(WDES, KI, KK, KB, VOVO_S, VOVO_S_R)

               DO KA = 1, MY_NKPTS

                  KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                           WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                  KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)

                  VVOO_S_R(:, :, :, :) = T2_R(:, :, :, :, KK, KJ, KA)

                  CALL SORT_V1V2O1O2_V2O1V1O2(WDES, VVOO_S, VOVO2_S, VVOO_S_R, VOVO2_S_R)

                  CALL SGEMM('t', 'n', ((NUNOCC)*(VBMAX)), ((NUNOCC)*VBMAX), &
                             ((NUNOCC)*VBMAX), REAL(KW, kind=qs), VOVO_S_R(1, 1, 1, 1), ((NUNOCC)*VBMAX), &
                             VOVO2_S_R(1, 1, 1, 1), ((NUNOCC)*VBMAX), &
                             (1.0_qs), T2_KTMP2_R(1, 1, 1, 1, KA), ((NUNOCC)*VBMAX))

                  IF ((PRECONDITION) .and. (KI == KK)) THEN
                     DO NA = 1, NUNOCC
                     DO NB = 1, NUNOCC
                        DO NI = 1, VBMAX
                        DO NJ = 1, VBMAX
                           PRET2(NA, NB, NI, NJ, KI, KJ, KA) = PRET2(NA, NB, NI, NJ, KI, KJ, KA) - CONJG(VOVO_S(NB, NI, NA, NJ))*KW
                        END DO
                        END DO
                     END DO
                     END DO
                  END IF

               END DO !KA
            END DO !KK

            DO KA = 1, MY_NKPTS
               KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                        WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
               KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)

               VOVO3_S_R(:, :, :, :) = T2_KTMP2_R(:, :, :, :, KA)
               CALL SORT_V2O1V1O2_V1V2O1O2(WDES, VOVO3_S, VVOO_S, VOVO3_S_R, VVOO_S_R)
               T2_N_R(:, :, :, :, KI, KJ, KA) = T2_N_R(:, :, :, :, KI, KJ, KA) - VVOO_S_R(:, :, :, :)
            END DO !KA
         END DO
         END DO !KB

      ELSE !LORBREAL

         DO KB = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
            DO KI = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
               T2_KTMP2(:, :, :, :, :) = (0.0_qs, 0.0_qs)
               DO KK = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KK) == 0)) CYCLE
                  CALL BCAST2ALL_CHI_CKAI_KJ_KK_KA(WDES, KI, KK, KB, VOVO_S, VOVO_S_R)

                  DO KA = 1, MY_NKPTS

                     KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                              WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                     KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                              WDES%VKPT(:, KQ), KPOINTS_FULL)

                     VVOO_S(:, :, :, :) = T2(:, :, :, :, RKIofKI(KK), RKIofKI(KJ), KA)

                     CALL SORT_V1V2O1O2_V2O1V1O2(WDES, VVOO_S, VOVO2_S, VVOO_S_R, VOVO2_S_R)

                     CALL CGEMM('t', 'n', ((NUNOCC)*(VBMAX)), ((NUNOCC)*VBMAX), &
                                ((NUNOCC)*VBMAX), KW, VOVO_S(1, 1, 1, 1), ((NUNOCC)*VBMAX), &
                                VOVO2_S(1, 1, 1, 1), ((NUNOCC)*VBMAX), &
                                (1.0_qs, 0.0_qs), T2_KTMP2(1, 1, 1, 1, KA), ((NUNOCC)*VBMAX))

                     IF ((PRECONDITION) .and. (KI == KK)) THEN
                        DO NA = 1, NUNOCC
                        DO NB = 1, NUNOCC
                           DO NI = 1, VBMAX
                           DO NJ = 1, VBMAX
           PRET2(NA,NB,NI,NJ,RKIofKI(KI),RKIofKI(KJ),KA)=PRET2(NA,NB,NI,NJ,RKIofKI(KI),RKIofKI(KJ),KA)-CONJG(VOVO_S(NB,NI,NA,NJ))*KW
                           END DO
                           END DO
                        END DO
                        END DO
                     END IF

                  END DO !KA
               END DO !KK

               DO KA = 1, MY_NKPTS
                  KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                           WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                  KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)

                  VOVO3_S(:, :, :, :) = T2_KTMP2(:, :, :, :, KA)
                  CALL SORT_V2O1V1O2_V1V2O1O2(WDES, VOVO3_S, VVOO_S, VOVO3_S_R, VVOO_S_R)
                T2_N(:, :, :, :, RKIofKI(KI), RKIofKI(KJ), KA) = T2_N(:, :, :, :, RKIofKI(KI), RKIofKI(KJ), KA) - VVOO_S(:, :, :, :)
               END DO !KA
            END DO
         END DO !KB

      END IF !LORBREAL

   END SUBROUTINE CONTR_CHI_BKCI_T2_ACKJ

!***********************************************************************
!
! The following routine swaps the T2 amplitudes between T2_TMP and T2.
! N.B.: The memory footprint of T2 is very large.
!
!***********************************************************************

   SUBROUTINE SWAP_T2(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
                 KC, KB_
      COMPLEX(qs) :: KW, TMP

      KW = WTKPT*(1.0_qs, 0.0_qs)
      KW = WTKPT*(1.0_qs, 0.0_qs)

      IF (LORBREAL) THEN

         DO KB_ = 1, WDES%NKPTS !MY_NKPTS
            CALL BCAST2ALL_T2(WDES, KB_, T2_TMP, T2_TMP_R)
            DO KA = 1, MY_NKPTS
            DO KJ = 1, WDES%NKPTS
               DO KI = 1, WDES%NKPTS
                  KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                           WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                  KB = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KJ) + &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)
                  IF (KB == KB_) THEN

                     DO NI = 1, VBMAX
                     DO NJ = 1, VBMAX
                     DO NA = 1, (NUNOCC)
                     DO NB = 1, (NUNOCC)
                        CHI_CKAI_R(NA, NI, NB, NJ, KI, KJ, KA) = T2_TMP_R(NB, NA, NJ, NI, KJ, KI)
                     END DO
                     END DO
                     END DO
                     END DO

                  END IF

               END DO
            END DO !KJ
            END DO !KA
         END DO !KB

         DO KA = 1, MY_NKPTS
            DO KJ = 1, WDES%NKPTS
               DO KI = 1, WDES%NKPTS
                  KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                           WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                  KB = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KJ) + &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)
                  DO NI = 1, VBMAX
                  DO NJ = 1, VBMAX
                  DO NA = 1, (NUNOCC)
                  DO NB = 1, (NUNOCC)
                     T2_R(NA, NB, NI, NJ, KI, KJ, KA) = CHI_CKAI_R(NA, NI, NB, NJ, KI, KJ, KA)
                  END DO
                  END DO
                  END DO
                  END DO

               END DO
            END DO !KJ
         END DO !KA

      ELSE !LORBREAL

         DO KB_ = 1, WDES%NKPTS !MY_NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB_) == 0)) CYCLE
            CALL BCAST2ALL_T2(WDES, KB_, T2_TMP, T2_TMP_R)
            DO KA = 1, MY_NKPTS
            DO KJ = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE
               DO KI = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
                  KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                           WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                  KB = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KJ) + &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)
                  IF (KB == KB_) THEN

                     DO NI = 1, VBMAX
                     DO NJ = 1, VBMAX
                     DO NA = 1, (NUNOCC)
                     DO NB = 1, (NUNOCC)
                        CHI_CKAI(NA, NI, NB, NJ, RKIofKI(KI), RKIofKI(KJ), KA) = T2_TMP(NB, NA, NJ, NI, RKIofKI(KJ), RKIofKI(KI))
                     END DO
                     END DO
                     END DO
                     END DO

                  END IF

               END DO
            END DO !KJ
            END DO !KA
         END DO !KB

         DO KA = 1, MY_NKPTS
            DO KJ = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE
               DO KI = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
                  KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                           WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                  KB = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KJ) + &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)
                  DO NI = 1, VBMAX
                  DO NJ = 1, VBMAX
                  DO NA = 1, (NUNOCC)
                  DO NB = 1, (NUNOCC)
                     T2(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA) = CHI_CKAI(NA, NI, NB, NJ, RKIofKI(KI), RKIofKI(KJ), KA)
                  END DO
                  END DO
                  END DO
                  END DO

               END DO
            END DO !KJ
         END DO !KA

      END IF

   END SUBROUTINE SWAP_T2

!***********************************************************************
!
! The following routine swaps the T2 amplitudes between T2_TMP and T2_N.
! N.B.: The memory footprint of T2 is very large.
!
!***********************************************************************


   SUBROUTINE SWAP_T2_N(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
                 KC, KB_
      COMPLEX(qs) :: KW, TMP

      KW = WTKPT*(1.0_qs, 0.0_qs)
      KW = WTKPT*(1.0_qs, 0.0_qs)

      IF (LORBREAL) THEN

         DO KB_ = 1, WDES%NKPTS !MY_NKPTS
            CALL BCAST2ALL_T2_N(WDES, KB_, T2_TMP, T2_TMP_R)
            DO KA = 1, MY_NKPTS
            DO KJ = 1, WDES%NKPTS
               DO KI = 1, WDES%NKPTS
                  KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                           WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                  KB = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KJ) + &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)
                  IF (KB == KB_) THEN

                     DO NI = 1, VBMAX
                     DO NJ = 1, VBMAX
                     DO NA = 1, (NUNOCC)
                     DO NB = 1, (NUNOCC)
                        CHI_CKAI_R(NA, NI, NB, NJ, KI, KJ, KA) = T2_TMP_R(NB, NA, NJ, NI, KJ, KI)
                     END DO
                     END DO
                     END DO
                     END DO

                  END IF

               END DO
            END DO !KJ
            END DO !KA
         END DO !KB

         DO KA = 1, MY_NKPTS
            DO KJ = 1, WDES%NKPTS
               DO KI = 1, WDES%NKPTS
                  KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                           WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                  KB = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KJ) + &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)
                  DO NI = 1, VBMAX
                  DO NJ = 1, VBMAX
                  DO NA = 1, (NUNOCC)
                  DO NB = 1, (NUNOCC)
                     T2_N_R(NA, NB, NI, NJ, KI, KJ, KA) = CHI_CKAI_R(NA, NI, NB, NJ, KI, KJ, KA)
                  END DO
                  END DO
                  END DO
                  END DO

               END DO
            END DO !KJ
         END DO !KA

      ELSE !LORBREAL
         DO KB_ = 1, WDES%NKPTS !MY_NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB_) == 0)) CYCLE
            CALL BCAST2ALL_T2_N(WDES, KB_, T2_TMP, T2_TMP_R)
            DO KA = 1, MY_NKPTS
            DO KJ = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE
               DO KI = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
                  KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                           WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                  KB = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KJ) + &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)
                  IF (KB == KB_) THEN

                     DO NI = 1, VBMAX
                     DO NJ = 1, VBMAX
                     DO NA = 1, (NUNOCC)
                     DO NB = 1, (NUNOCC)
                        CHI_CKAI(NA, NI, NB, NJ, RKIofKI(KI), RKIofKI(KJ), KA) = T2_TMP(NB, NA, NJ, NI, RKIofKI(KJ), RKIofKI(KI))
                     END DO
                     END DO
                     END DO
                     END DO

                  END IF

               END DO
            END DO !KJ
            END DO !KA
         END DO !KB

         DO KA = 1, MY_NKPTS
            DO KJ = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE
               DO KI = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
                  KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                           WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                  KB = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KJ) + &
                                           WDES%VKPT(:, KQ), KPOINTS_FULL)
                  DO NI = 1, VBMAX
                  DO NJ = 1, VBMAX
                  DO NA = 1, (NUNOCC)
                  DO NB = 1, (NUNOCC)
                     T2_N(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA) = CHI_CKAI(NA, NI, NB, NJ, RKIofKI(KI), RKIofKI(KJ), KA)
                  END DO
                  END DO
                  END DO
                  END DO

               END DO
            END DO !KJ
         END DO !KA
      END IF

   END SUBROUTINE SWAP_T2_N

!***********************************************************************
!
! The following routine swaps the T2 amplitudes between T2_TMP and T2.
! N.B.: The memory footprint of T2 is very large.
!
!***********************************************************************

   SUBROUTINE P_T2_N(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
                 KC, KB_
      COMPLEX(qs) :: KW

      KW = WTKPT*(1.0_qs, 0.0_qs)

      IF (LORBREAL) THEN

         DO KB = 1, WDES%NKPTS
            CALL BCAST2ALL_T2_N(WDES, KB, T2_TMP, T2_TMP_R)
            DO KA = 1, MY_NKPTS
               DO KJ = 1, WDES%NKPTS
                  DO KI = 1, WDES%NKPTS
                     KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                              WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                     KB_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KJ) + &
                                               WDES%VKPT(:, KQ), KPOINTS_FULL)
                     IF (KB == KB_) THEN
                        DO NI = 1, VBMAX
                        DO NJ = 1, VBMAX
                        DO NA = 1, (NUNOCC)
                        DO NB = 1, (NUNOCC)

                           CHI_CKAI_R(NA, NI, NB, NJ, KI, KJ, KA) = T2_N_R(NA, NB, NI, NJ, KI, KJ, KA) + &
                                                                    T2_TMP_R(NB, NA, NJ, NI, KJ, KI)
                        END DO
                        END DO
                        END DO
                        END DO
                     END IF
                  END DO
               END DO !KJ
            END DO !KA
         END DO !KB

         DO KA = 1, MY_NKPTS
         DO KJ = 1, WDES%NKPTS
            DO KI = 1, WDES%NKPTS
               DO NI = 1, VBMAX
                  DO NJ = 1, VBMAX
                  DO NA = 1, (NUNOCC)
                  DO NB = 1, (NUNOCC)

                     T2_N_R(NA, NB, NI, NJ, KI, KJ, KA) = CHI_CKAI_R(NA, NI, NB, NJ, KI, KJ, KA)

                  END DO
                  END DO
                  END DO
               END DO
            END DO
         END DO !KJ
         END DO !KA

      ELSE !LORBREAL

         DO KB = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
            CALL BCAST2ALL_T2_N(WDES, KB, T2_TMP, T2_TMP_R)
            DO KA = 1, MY_NKPTS
               DO KJ = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE
                  DO KI = 1, WDES%NKPTS
                     IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
                     KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                              WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                     KB_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KJ) + &
                                               WDES%VKPT(:, KQ), KPOINTS_FULL)
                     IF (KB == KB_) THEN
                        DO NI = 1, VBMAX
                        DO NJ = 1, VBMAX
                        DO NA = 1, (NUNOCC)
                        DO NB = 1, (NUNOCC)

                     CHI_CKAI(NA, NI, NB, NJ, RKIofKI(KI), RKIofKI(KJ), KA) = T2_N(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA) + &
                                                                                    T2_TMP(NB, NA, NJ, NI, RKIofKI(KJ), RKIofKI(KI))
                        END DO
                        END DO
                        END DO
                        END DO
                     END IF
                  END DO
               END DO !KJ
            END DO !KA
         END DO !KB

         DO KA = 1, MY_NKPTS
         DO KJ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE
            DO KI = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
               DO NI = 1, VBMAX
               DO NJ = 1, VBMAX
               DO NA = 1, (NUNOCC)
               DO NB = 1, (NUNOCC)

                  T2_N(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA) = CHI_CKAI(NA, NI, NB, NJ, RKIofKI(KI), RKIofKI(KJ), KA)

               END DO
               END DO
               END DO
               END DO
            END DO
         END DO !KJ
         END DO !KA

      END IF

   END SUBROUTINE P_T2_N

   SUBROUTINE UP_T2_N(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
                 KC, KB_
      COMPLEX(qs) :: KW

      KW = KPOINTS_FULL%WTKPT(1)*(1.0_qs, 0.0_qs)

      KW = WTKPT*(1.0_qs, 0.0_qs)

      IF (LORBREAL) THEN

         DO KB = 1, WDES%NKPTS
            CALL BCAST2ALL_T2_N(WDES, KB, T2_TMP, T2_TMP_R)
            DO KA = 1, MY_NKPTS
               DO KJ = 1, WDES%NKPTS
                  DO KI = 1, WDES%NKPTS
                     !T2_KTMP2(:,:,:,:,:)=(0.0_qs,0.0_qs)
                     KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                              WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                     KB_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KJ) + &
                                               WDES%VKPT(:, KQ), KPOINTS_FULL)
                     IF (KB == KB_) THEN
                        DO NI = 1, VBMAX
                        DO NJ = 1, VBMAX
                        DO NA = 1, (NUNOCC)
                        DO NB = 1, (NUNOCC)
                           CHI_CKAI_R(NA, NI, NB, NJ, KI, KJ, KA) = T2_TMP_R(NB, NA, NJ, NI, KJ, KI)
                        END DO
                        END DO
                        END DO
                        END DO
                     END IF
                  END DO
               END DO !KJ
            END DO !KA
         END DO !KB

         DO KA = 1, MY_NKPTS
            DO KI = 1, WDES%NKPTS
               DO KJ = 1, WDES%NKPTS
                  DO NI = 1, VBMAX
                     DO NJ = 1, VBMAX
                     DO NA = 1, (NUNOCC)
                     DO NB = 1, (NUNOCC)
                        IF (KJ > KI) T2_N_R(NA, NB, NI, NJ, KI, KJ, KA) = CHI_CKAI_R(NA, NI, NB, NJ, KI, KJ, KA)
                     END DO
                     END DO
                     END DO
                  END DO
               END DO
            END DO !KJ
         END DO !KA

      ELSE !LORBREAL
         DO KB = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
            CALL BCAST2ALL_T2_N(WDES, KB, T2_TMP, T2_TMP_R)
            DO KA = 1, MY_NKPTS
               DO KJ = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE
                  DO KI = 1, WDES%NKPTS
                     IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
                     KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                              WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
                     KB_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KJ) + &
                                               WDES%VKPT(:, KQ), KPOINTS_FULL)
                     IF (KB == KB_) THEN
                        DO NI = 1, VBMAX
                        DO NJ = 1, VBMAX
                        DO NA = 1, (NUNOCC)
                        DO NB = 1, (NUNOCC)
                           CHI_CKAI(NA, NI, NB, NJ, RKIofKI(KI), RKIofKI(KJ), KA) = T2_TMP(NB, NA, NJ, NI, RKIofKI(KJ), RKIofKI(KI))

                        END DO
                        END DO
                        END DO
                        END DO

                     END IF
                  END DO
               END DO !KJ
            END DO !KA
         END DO !KB

         DO KA = 1, MY_NKPTS
            DO KI = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
               DO KJ = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE
                  DO NI = 1, VBMAX
                     DO NJ = 1, VBMAX
                     DO NA = 1, (NUNOCC)
                     DO NB = 1, (NUNOCC)
            IF (KJ > KI) T2_N(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA) = CHI_CKAI(NA, NI, NB, NJ, RKIofKI(KI), RKIofKI(KJ), KA)
                     END DO
                     END DO
                     END DO
                  END DO
               END DO
            END DO !KJ
         END DO !KA

      END IF

   END SUBROUTINE UP_T2_N

!***********************************************************************
!
! The following routine contracts an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************


   SUBROUTINE CONTR_K_AC_T1(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
         KC
      COMPLEX(qs) :: KW

      IF (LORBREAL) THEN

         DO KI = 1, WDES%NKPTS

            VV2_S_R(:, :) = K_AC_R(:, :, KI) !(1.0_qs,0.0_qs)
            CALL SORT_V1V2_V2V1(WDES, VV2_S, VV_S, VV2_S_R, VV_S_R)

            CALL SGEMM('t', 'n', (NUNOCC), (VBMAX), &
                       (NUNOCC), (1._qs), VV_S_R(1, 1), (NUNOCC), &
                       T1_R(1, 1, KI), (NUNOCC), &
                       (1._qs), T1_N_R(1, 1, KI), (NUNOCC))

         END DO

      ELSE

         DO KI = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE

            VV2_S(:, :) = K_AC(:, :, RKIofKI(KI)) !(1.0_qs,0.0_qs)
            CALL SORT_V1V2_V2V1(WDES, VV2_S, VV_S, VV2_S_R, VV_S_R)

            CALL CGEMM('t', 'n', (NUNOCC), (VBMAX), &
                       (NUNOCC), (1._qs, 0._qs), VV_S(1, 1), (NUNOCC), &
                       T1(1, 1, RKIofKI(KI)), (NUNOCC), &
                       (1._qs, 0._qs), T1_N(1, 1, RKIofKI(KI)), (NUNOCC))

         END DO

      END IF

   END SUBROUTINE CONTR_K_AC_T1

!***********************************************************************
!
! The following routine contracts an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************


   SUBROUTINE CONTR_F_KC_T1T1(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
         KC
      COMPLEX(qs) :: KW

      IF (LORBREAL) THEN

         DO KI = 1, WDES%NKPTS

            VO_S_R(:, :) = T1_R(:, :, KI)
            CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)
            OV2_S_R(:, :) = F_KC_R(:, :, KI)
            CALL SORT_O1V1_V1O1(WDES, OV2_S, VO_S, OV2_S_R, VO_S_R)

            CALL SGEMM('C', 'n', (VBMAX), (VBMAX), &
                       (NUNOCC), (1._qs), VO_S_R(1, 1), (NUNOCC), &
                       T1_R(1, 1, KI), (NUNOCC), &
                       (0._qs), TE4O_R(1, 1), (VBMAX))

            CALL SGEMM('t', 'n', (NUNOCC), (VBMAX), &
                       (VBMAX), (-2._qs), OV_S_R(1, 1), (VBMAX), &
                       TE4O_R(1, 1), (VBMAX), &
                       (1._qs), T1_N_R(1, 1, KI), (NUNOCC))

         END DO

      ELSE

         DO KI = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE

            VO_S(:, :) = T1(:, :, RKIofKI(KI))
            CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)
            OV2_S(:, :) = F_KC(:, :, RKIofKI(KI))
            CALL SORT_O1V1_V1O1(WDES, OV2_S, VO_S, OV2_S_R, VO_S_R)

            CALL CGEMM('C', 'n', (VBMAX), (VBMAX), &
                       (NUNOCC), (1._qs, 0._qs), VO_S(1, 1), (NUNOCC), &
                       T1(1, 1, RKIofKI(KI)), (NUNOCC), &
                       (0._qs, 0._qs), TE4O(1, 1), (VBMAX))

            CALL CGEMM('t', 'n', (NUNOCC), (VBMAX), &
                       (VBMAX), (-2._qs, 0._qs), OV_S(1, 1), (VBMAX), &
                       TE4O(1, 1), (VBMAX), &
                       (1._qs, 0._qs), T1_N(1, 1, RKIofKI(KI)), (NUNOCC))

         END DO

      END IF

   END SUBROUTINE CONTR_F_KC_T1T1

!***********************************************************************
!
! The following routine contracts an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE CONTR_K_KI_T1(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
         KC
      COMPLEX(qs) :: KW

      IF (LORBREAL) THEN

         DO KI = 1, WDES%NKPTS

            VO_S_R = T1_R(:, :, KI)
            CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)

            CALL SGEMM('t', 'n', (NUNOCC), (VBMAX), &
                       (VBMAX), (-1._qs), OV_S_R(1, 1), (VBMAX), &
                       K_KI_R(1, 1, KI), (VBMAX), &
                       (1._qs), T1_N_R(1, 1, KI), (NUNOCC))

         END DO

      ELSE

         DO KI = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE

            VO_S = T1(:, :, RKIofKI(KI))
            CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)

            CALL CGEMM('t', 'n', (NUNOCC), (VBMAX), &
                       (VBMAX), (-1._qs, 0._qs), OV_S(1, 1), (VBMAX), &
                       K_KI(1, 1, RKIofKI(KI)), (VBMAX), &
                       (1._qs, 0._qs), T1_N(1, 1, RKIofKI(KI)), (NUNOCC))

         END DO

      END IF

   END SUBROUTINE CONTR_K_KI_T1

!***********************************************************************
!
! The following routine contracts an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE CONTR_K_KC_T(WGW, WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      INTEGER KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
         KC, NC
      COMPLEX(qs) :: KW

      KW = WTKPT

      IF (LORBREAL) THEN

         T1_T_R = (0.0_qs)
         DO KK = 1, WDES%NKPTS
            DO KA = 1, MY_NKPTS

               DO NA = 1, (NUNOCC)
               DO NI = 1, VBMAX
                  VO_S_R(:, :) = (2.0_qs)*T2_R(NA, :, NI, :, KPTS_MKPTS(KA), KK, KA) - &
                                 T2_R(NA, :, :, NI, KK, KPTS_MKPTS(KA), KA)
                  IF (KK == KPTS_MKPTS(KA)) THEN
                     DO NC = 1, (NUNOCC)
                     DO NK = 1, VBMAX
                        VO_S_R(NC, NK) = VO_S_R(NC, NK) + T1_R(NC, NI, KPTS_MKPTS(KA))*T1_R(NA, NK, KPTS_MKPTS(KA))/ &
                                         KPOINTS_FULL%WTKPT(1)
                     END DO
                     END DO
                  END IF
                  CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)

                  CALL SGEMM('t', 'n', 1, 1, &
                             (NUNOCC)*VBMAX, KW, OV_S_R(1, 1), (NUNOCC)*VBMAX, &
                             K_KC_R(1, 1, KK), (NUNOCC)*VBMAX, &
                             (1.0_qs), T1_T_R(1, 1, KPTS_MKPTS(KA)), 1)

               END DO
               END DO
            END DO
         END DO

         CALL M_sum_single(WDES%COMM, T1_T_R, SIZE(T1_T_R))
         T1_N_R = T1_T_R + T1_N_R

      ELSE

         T1_T = (0.0_qs, 0.0_qs)
         DO KK = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KK) == 0)) CYCLE
            DO KA = 1, MY_NKPTS

               DO NA = 1, (NUNOCC)
               DO NI = 1, VBMAX
                  VO_S(:, :) = (2.0_qs, 0.0_qs)*T2(NA, :, NI, :, RKIofKI(KPTS_MKPTS(KA)), RKIofKI(KK), KA) - &
                               T2(NA, :, :, NI, RKIofKI(KK), RKIofKI(KPTS_MKPTS(KA)), KA)
                  IF (KK == KPTS_MKPTS(KA)) THEN
                     DO NC = 1, (NUNOCC)
                     DO NK = 1, VBMAX
                        VO_S(NC, NK) = VO_S(NC, NK) + T1(NC, NI, RKIofKI(KPTS_MKPTS(KA)))*T1(NA, NK, RKIofKI(KPTS_MKPTS(KA)))/ &
                                       WTKPT
                     END DO
                     END DO
                  END IF
                  CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)

                  CALL CGEMM('t', 'n', 1, 1, &
                             (NUNOCC)*VBMAX, KW, OV_S(1, 1), (NUNOCC)*VBMAX, &
                             K_KC(1, 1, RKIofKI(KK)), (NUNOCC)*VBMAX, &
                             (1.0_qs, 0.0_qs), T1_T(1, 1, RKIofKI(KPTS_MKPTS(KA))), 1)

               END DO
               END DO
            END DO
         END DO

         CALL M_sum_single(WDES%COMM, T1_T, 2*SIZE(T1_T))
         T1_N = T1_T + T1_N

      END IF

   END SUBROUTINE CONTR_K_KC_T

!***********************************************************************
!
! The following routine contracts an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE CONTR_W_AKCD_T(WGW, WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      INTEGER KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
         KC, NC, ND, KD
      COMPLEX(qs) :: KW2

      KW2 = (1.0_qs, 0.0_qs)*WTKPT*WTKPT
      IF (LORBREAL) THEN
         T1_T_R = (0.0_qs)
      ELSE
         T1_T = (0.0_qs, 0.0_qs)
      END IF
      DO KD = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KD) == 0)) CYCLE
         CALL BCAST2ALL_FTOD_AB(WDES, KD, 2, PW_AB_TMP, OC_AB_TMP)
         CALL BCAST2ALL_FTOD_AI(WDES, KD, 1, PW_AI_TMP, OC_AI_TMP)
         DO KA = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA) == 0)) CYCLE
            DO KC = 1, MY_NKPTS

               KQ = KPOINT_IN_FULL_GRID(-WDES%VKPT(:, KPTS_MKPTS(KC)) + &
                                        WDES%VKPT(:, KA), KPOINTS_FULL)
               KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KA) - &
                                         WDES%VKPT(:, KD), KPOINTS_FULL)
               KK = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KD) - &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)

               DO ND = 1, (NUNOCC)
               DO NC = 1, (NUNOCC)

                  CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                  PW_AI_TMP(1, 1, ND, RKQofKQ(KQ)), OC_AI_TMP(1, 1, ND, RKQofKQ(KQ)), &
                                  (VBMAX), &
                                  FTOD_PW_AB(1, 1, NC, RKQofKQ(KQ), KC, 2), FTOD_OC_AB(1, 1, NC, RKQofKQ(KQ), KC, 2), &
                                  (NUNOCC), OVVV(1, 1, NC, ND), OVVV_R(1, 1, NC, ND), zero)

                  IF (LORBREAL) THEN
                     OVVV_R(:, :, NC, ND) = -(2.0_q)*(OVVV_R(:, :, NC, ND))
                  ELSE
                     OVVV(:, :, NC, ND) = -(2.0_q, 0.0_q)*(OVVV(:, :, NC, ND))
                  END IF

                  CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                  FTOD_PW_AI(1, 1, NC, RKQofKQ(KQ_), KC, 1), FTOD_OC_AI(1, 1, NC, RKQofKQ(KQ_), KC, 1), &
                                  (VBMAX), &
                                  PW_AB_TMP(1, 1, ND, RKQofKQ(KQ_)), OC_AB_TMP(1, 1, ND, RKQofKQ(KQ_)), &
                                  (NUNOCC), OVVV(1, 1, NC, ND), OVVV_R(1, 1, NC, ND), one)

               END DO
               END DO

               IF (LORBREAL) THEN
                  DO NA = 1, (NUNOCC)
                  DO NK = 1, VBMAX
                     VV_S_R(:, :) = -(OVVV_R(NK, NA, :, :))
                     DO NI = 1, VBMAX
                        VV2_S_R(:, :) = T2_R(:, :, NI, NK, KA, KK, KC)
                        IF (KPTS_MKPTS(KC) == KA) THEN
                           DO ND = 1, (NUNOCC)
                           DO NC = 1, (NUNOCC)
                              VV2_S_R(NC, ND) = VV2_S_R(NC, ND) + T1_R(NC, NI, KA)*T1_R(ND, NK, KD)/KPOINTS_FULL%WTKPT(1)
                           END DO
                           END DO
                        END IF
                        CALL SGEMM('t', 'n', 1, 1, (NUNOCC)*(NUNOCC), &
                                   REAL(KW2, kind=qs), VV2_S_R, (NUNOCC)*(NUNOCC), &
                                   VV_S_R(1, 1), (NUNOCC)*(NUNOCC), &
                                   (1._qs), T1_T_R(NA, NI, KA), 1)
                     END DO
                  END DO
                  END DO
               ELSE
                  DO NA = 1, (NUNOCC)
                  DO NK = 1, VBMAX
                     VV_S(:, :) = -CONJG(OVVV(NK, NA, :, :))
                     DO NI = 1, VBMAX
                        VV2_S(:, :) = T2(:, :, NI, NK, RKIofKI(KA), RKIofKI(KK), KC)
                        IF (KPTS_MKPTS(KC) == KA) THEN
                           DO ND = 1, (NUNOCC)
                           DO NC = 1, (NUNOCC)
                              VV2_S(NC, ND) = VV2_S(NC, ND) + T1(NC, NI, RKIofKI(KA))*T1(ND, NK, RKIofKI(KD))/WTKPT
                           END DO
                           END DO
                        END IF
                        CALL CGEMM('t', 'n', 1, 1, (NUNOCC)*(NUNOCC), &
                                   KW2, VV2_S, (NUNOCC)*(NUNOCC), &
                                   VV_S(1, 1), (NUNOCC)*(NUNOCC), &
                                   (1._qs, 0._qs), T1_T(NA, NI, RKIofKI(KA)), 1)
                     END DO
                  END DO
                  END DO
               END IF
            END DO
         END DO
      END DO

      IF (LORBREAL) THEN
         CALL M_sum_single(WDES%COMM, T1_T_R, SIZE(T1_T_R))
         T1_N_R = T1_T_R + T1_N_R
      ELSE
         CALL M_sum_single(WDES%COMM, T1_T, 2*SIZE(T1_T))
         T1_N = T1_T + T1_N
      END IF

   END SUBROUTINE CONTR_W_AKCD_T

!***********************************************************************
!
! The following routine contracts an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE CONTR_W_AKIC_T1(WGW, WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      INTEGER KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
         KC, NC, ND, KD
      COMPLEX(qs) :: KW

      KW = (1.0_qs, 0.0_qs)*WTKPT
      IF (LORBREAL) THEN
         T1_T_R = (0.0_qs)
      ELSE
         T1_T = (0.0_qs, 0.0_qs)
      END IF
      DO KC = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KC) == 0)) CYCLE
         CALL BCAST2ALL_FTOD_AB(WDES, KC, 2, PW_AB_TMP, OC_AB_TMP)
         CALL BCAST2ALL_FTOD_AI(WDES, KC, 2, PW_AI_TMP, OC_AI_TMP)
         DO KA = 1, MY_NKPTS

            KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                     WDES%VKPT(:, KPTS_MKPTS(KA)), KPOINTS_FULL)
            KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                      WDES%VKPT(:, KC), KPOINTS_FULL)

            DO NI = 1, VBMAX
            DO NC = 1, (NUNOCC)

               CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                               FTOD_PW_IA(1, 1, NI, RKQofKQ(KQ), KA, 1), FTOD_OC_IA(1, 1, NI, RKQofKQ(KQ), KA, 1), &
                               (NUNOCC), &
                               PW_AI_TMP(1, 1, NC, RKQofKQ(KQ)), OC_AI_TMP(1, 1, NC, RKQofKQ(KQ)), &
                               (VBMAX), VOVO(1, 1, NC, NI), VOVO_R(1, 1, NC, NI), zero)
               IF (LORBREAL) THEN
                  VOVO_R(:, :, NC, NI) = -(2.0_q)*(VOVO_R(:, :, NC, NI))
               ELSE
                  VOVO(:, :, NC, NI) = -(2.0_q, 0.0_q)*CONJG(VOVO(:, :, NC, NI))
               END IF
               CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                               PW_AB_TMP(1, 1, NC, RKQofKQ(KQ_)), OC_AB_TMP(1, 1, NC, RKQofKQ(KQ_)), &
                               (NUNOCC), &
                               FTOD_PW_IJ(1, 1, NI, RKQofKQ(KQ_), KA, 1), FTOD_OC_IJ(1, 1, NI, RKQofKQ(KQ_), KA, 1), &
                               (VBMAX), VOVO(1, 1, NC, NI), VOVO_R(1, 1, NC, NI), one)

            END DO
            END DO

            IF (LORBREAL) THEN
               VOVO_S_R = -VOVO_R
            ELSE
               VOVO_S = -VOVO
            END IF
            CALL SORT_V1O1V2O2_V2O1V1O2(WDES, VOVO_S, VOVO2_S, VOVO_S_R, VOVO2_S_R)

            IF (LORBREAL) THEN
               CALL SGEMM('t', 'n', (NUNOCC)*(VBMAX), 1, (NUNOCC)*(VBMAX), &
                          REAL(KW, kind=qs), VOVO2_S_R(1, 1, 1, 1), (NUNOCC)*(VBMAX), &
                          T1_R(1, 1, KC), (NUNOCC)*(VBMAX), &
                          (1._qs), T1_T_R(1, 1, KPTS_MKPTS(KA)), (NUNOCC)*(VBMAX))
            ELSE
               CALL CGEMM('t', 'n', (NUNOCC)*(VBMAX), 1, (NUNOCC)*(VBMAX), &
                          KW, VOVO2_S(1, 1, 1, 1), (NUNOCC)*(VBMAX), &
                          T1(1, 1, RKIofKI(KC)), (NUNOCC)*(VBMAX), &
                          (1._qs, 0._qs), T1_T(1, 1, RKIofKI(KPTS_MKPTS(KA))), (NUNOCC)*(VBMAX))
            END IF

         END DO
      END DO

      IF (LORBREAL) THEN
         CALL M_sum_single(WDES%COMM, T1_T_R, SIZE(T1_T_R))
         T1_N_R = T1_T_R + T1_N_R
      ELSE
         CALL M_sum_single(WDES%COMM, T1_T, 2*SIZE(T1_T))
         T1_N = T1_T + T1_N
      END IF

   END SUBROUTINE CONTR_W_AKIC_T1

!***********************************************************************
!
! The following routine contracts an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE CONTR_W_KLIC_T(WGW, WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      INTEGER KI, KJ, KK, KL, NI, NJ, NK, NL, KQ, KQ_, NA, NB, KA, KB, &
         KC, NC, ND, KD
      COMPLEX(qs) :: KW2

      KW2 = (1.0_qs, 0.0_qs)*WTKPT*WTKPT
      IF (LORBREAL) THEN
         T1_T_R = (0.0_qs)
      ELSE
         T1_T = (0.0_qs, 0.0_qs)
      END IF

      DO KC = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KC) == 0)) CYCLE
         CALL BCAST2ALL_FTOD_AI(WDES, KC, 2, PW_AI_TMP, OC_AI_TMP)
         DO KK = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KK) == 0)) CYCLE
            DO KI = 1, MY_NKPTS

               KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KI)) - &
                                        WDES%VKPT(:, KK), KPOINTS_FULL)
               KQ_ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KK) - &
                                         WDES%VKPT(:, KC), KPOINTS_FULL)
               KL = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KC) + &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)

               DO NI = 1, VBMAX
               DO NC = 1, (NUNOCC)

                  CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                  FTOD_PW_IJ(1, 1, NI, RKQofKQ(KQ), KI, 1), FTOD_OC_IJ(1, 1, NI, RKQofKQ(KQ), KI, 1), &
                                  (VBMAX), &
                                  PW_AI_TMP(1, 1, NC, RKQofKQ(KQ)), OC_AI_TMP(1, 1, NC, RKQofKQ(KQ)), &
                                  (VBMAX), OOOV(1, 1, NI, NC), OOOV_R(1, 1, NI, NC), zero)
                  IF (LORBREAL) THEN
                     OOOV_R(:, :, NI, NC) = -(2.0_q)*(OOOV_R(:, :, NI, NC))
                  ELSE
                     OOOV(:, :, NI, NC) = -(2.0_q, 0.0_q)*CONJG(OOOV(:, :, NI, NC))
                  END IF

                  CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                  PW_AI_TMP(1, 1, NC, RKQofKQ(KQ_)), OC_AI_TMP(1, 1, NC, RKQofKQ(KQ_)), &
                                  (VBMAX), &
                                  FTOD_PW_IJ(1, 1, NI, RKQofKQ(KQ_), KI, 1), FTOD_OC_IJ(1, 1, NI, RKQofKQ(KQ_), KI, 1), &
                                  (VBMAX), OOOV(1, 1, NI, NC), OOOV_R(1, 1, NI, NC), one)

               END DO
               END DO
               IF (LORBREAL) THEN
                  OOOV_S_R(:, :, :, :) = (OOOV_R(:, :, :, :))
               ELSE
                  OOOV_S(:, :, :, :) = (OOOV(:, :, :, :))
               END IF

               IF (LORBREAL) THEN
                  DO NC = 1, (NUNOCC)
                  DO NA = 1, (NUNOCC)
                     OOVV2_S_R(:, :, 1, 1) = T2_R(NA, NC, :, :, KK, KL, KI)
                     IF (KK == KPTS_MKPTS(KI)) THEN
                        DO NL = 1, VBMAX
                           OOVV2_S_R(:, NL, 1, 1) = OOVV2_S_R(:, NL, 1, 1) + T1_R(NA, :, KK)*T1_R(NC, NL, KL)/ &
                                                    ((1.0_qs)*KPOINTS_FULL%WTKPT(1))
                        END DO
                     END IF
                     DO NI = 1, VBMAX
                        CALL SGEMM('t', 'n', 1, 1, (VBMAX)*(VBMAX), &
                                   REAL(KW2, kind=qs), OOVV2_S_R(1, 1, 1, 1), (VBMAX)*(VBMAX), &
                                   OOOV_S_R(1, 1, NI, NC), (VBMAX)*(VBMAX), &
                                   (1._qs), T1_T_R(NA, NI, KPTS_MKPTS(KI)), 1)
                     END DO
                  END DO
                  END DO

               ELSE !LORBREAL

                  DO NC = 1, (NUNOCC)
                  DO NA = 1, (NUNOCC)
                     OOVV2_S(:, :, 1, 1) = T2(NA, NC, :, :, RKIofKI(KK), RKIofKI(KL), KI)
                     IF (KK == KPTS_MKPTS(KI)) THEN
                        DO NL = 1, VBMAX
                           OOVV2_S(:, NL, 1, 1) = OOVV2_S(:, NL, 1, 1) + T1(NA, :, RKIofKI(KK))*T1(NC, NL, RKIofKI(KL))/ &
                                                  ((1.0_qs, 0.0_qs)*WTKPT)
                        END DO
                     END IF
                     DO NI = 1, VBMAX
                        CALL CGEMM('t', 'n', 1, 1, (VBMAX)*(VBMAX), &
                                   KW2, OOVV2_S(1, 1, 1, 1), (VBMAX)*(VBMAX), &
                                   OOOV_S(1, 1, NI, NC), (VBMAX)*(VBMAX), &
                                   (1._qs, 0._qs), T1_T(NA, NI, RKIofKI(KPTS_MKPTS(KI))), 1)
                     END DO
                  END DO
                  END DO
               END IF

            END DO
         END DO
      END DO

      IF (LORBREAL) THEN
         CALL M_sum_single(WDES%COMM, T1_T_R, SIZE(T1_T_R))
         T1_N_R = T1_T_R + T1_N_R
      ELSE
         CALL M_sum_single(WDES%COMM, T1_T, 2*SIZE(T1_T))
         T1_N = T1_T + T1_N
      END IF

   END SUBROUTINE CONTR_W_KLIC_T

!***********************************************************************
!
! The following routine computes and contracts an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE FORM_CONTR_VT1_VT1T1_C(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER KI, KJ, KA, KB, NI, NJ, NK, NL, KQ, KQ_, KC, KD, &
         NC, ND, NA, NB

      call BLACS_PINFO(ME, PROCS)

      DO KB = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
         CALL BCAST2ALL_FTOD_AB(WDES, KB, 2, PW_AB_TMP, OC_AB_TMP)
         DO KI = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
            DO KA = 1, MY_NKPTS
               KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                        WDES%VKPT(:, KI), KPOINTS_FULL)

               KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) + &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)

               DO NB = 1, (NUNOCC)

                  CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                  FTOD_PW_IJ(1, 1, 1, RKQofKQ(KQ), KA, 1), FTOD_OC_IJ(1, 1, 1, RKQofKQ(KQ), KA, 1), &
                                  (VBMAX*VBMAX), &
                                  PW_AB_TMP(1, 1, NB, RKQofKQ(KQ)), OC_AB_TMP(1, 1, NB, RKQofKQ(KQ)), &
                                  (NUNOCC), OOVV(1, 1, 1, 1), OOVV_R(1, 1, 1, 1), zero)
                  !oovv : ikc_
                  IF (LORBREAL) THEN
                     DO NI = 1, VBMAX
                        VO_S_R = T1_R(:, :, KPTS_MKPTS(KA))
                        CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)
                        OVVO_S_R(:, :, 1, 1) = OV_S_R(:, :)
                        OV_S_R(:, :) = OOVV_R(NI, :, :, 1)

                        CALL SGEMM('t', 'n', (NUNOCC), (NUNOCC), &
                                   (VBMAX), (1._qs), OV_S_R(1, 1), &
                                   (VBMAX), &
                                   OVVO_S_R(1, 1, 1, 1), (VBMAX), &
                                   (0._qs), VVOO_S_R(1, 1, NI, 1), (NUNOCC))
                        !vvoo : cai_
                     END DO
                  ELSE
                     DO NI = 1, VBMAX
                        VO_S = T1(:, :, RKIofKI(KPTS_MKPTS(KA)))
                        CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)
                        OVVO_S(:, :, 1, 1) = OV_S(:, :)
                        OV_S(:, :) = OOVV(NI, :, :, 1)

                        CALL CGEMM('t', 'n', (NUNOCC), (NUNOCC), &
                                   (VBMAX), (1._qs, 0._qs), OV_S(1, 1), &
                                   (VBMAX), &
                                   OVVO_S(1, 1, 1, 1), (VBMAX), &
                                   (0._qs, 0._qs), VVOO_S(1, 1, NI, 1), (NUNOCC))
                        !vvoo : cai_
                     END DO
                  END IF

                  DO NA = 1, (NUNOCC)

                     CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                     FTOD_PW_AI(1, 1, NA, RKQofKQ(KQ), KA, 1), FTOD_OC_AI(1, 1, NA, RKQofKQ(KQ), KA, 1), &
                                     (VBMAX), &
                                     PW_AB_TMP(1, 1, NB, RKQofKQ(KQ)), OC_AB_TMP(1, 1, NB, RKQofKQ(KQ)), &
                                     (NUNOCC), OVVO(1, 1, 1, 1), OVVO_R(1, 1, 1, 1), zero)

                     IF (LORBREAL) THEN
                        !ov : ic -> ci
                        OV_S_R(:, :) = OVVO_R(:, :, 1, 1)
                        CALL SORT_O1V1_V1O1(WDES, OV_S, VO_S, OV_S_R, VO_S_R)

                        DO NI = 1, VBMAX
                           VO_S_R(:, NI) = VO_S_R(:, NI) - VVOO_S_R(:, NA, NI, 1)

                           DO NJ = 1, VBMAX
                              CALL SGEMM('t', 'n', 1, 1, &
                                         (NUNOCC), (1._qs), VO_S_R(1, NI), &
                                         (NUNOCC), &
                                         T1_R(1, NJ, KJ), (NUNOCC), &
                                         (1._qs), T2_N_R(NA, NB, NI, NJ, KI, KJ, KA), 1)
                           END DO
                        END DO
                     ELSE
                        !ov : ic -> ci
                        OV_S(:, :) = OVVO(:, :, 1, 1)
                        CALL SORT_O1V1_V1O1(WDES, OV_S, VO_S, OV_S_R, VO_S_R)

                        DO NI = 1, VBMAX
                           VO_S(:, NI) = VO_S(:, NI) - VVOO_S(:, NA, NI, 1)

                           DO NJ = 1, VBMAX
                              CALL CGEMM('t', 'n', 1, 1, &
                                         (NUNOCC), (1._qs, 0._qs), VO_S(1, NI), &
                                         (NUNOCC), &
                                         T1(1, NJ, RKIofKI(KJ)), (NUNOCC), &
                                         (1._qs, 0._qs), T2_N(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA), 1)
                           END DO
                        END DO
                     END IF
                  END DO !NA
               END DO !NB

            END DO
         END DO
      END DO

   END SUBROUTINE FORM_CONTR_VT1_VT1T1_C

!***********************************************************************
!
! The following routine computes and contracts an important intermediate quantitiy needed in the
! CCSD amplitude equations.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   SUBROUTINE FORM_CONTR_VT1_VT1T1_K(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER KI, KJ, KA, KB, NI, NJ, NK, NL, KQ, KQ_, KC, KD, &
         NC, ND, NA, NB

      call BLACS_PINFO(ME, PROCS)

      DO KB = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
         CALL BCAST2ALL_FTOD_IJ(WDES, KB, 2, PW_IJ_TMP, OC_IJ_TMP)
         CALL BCAST2ALL_FTOD_IA(WDES, KB, 2, PW_IA_TMP, OC_IA_TMP)
         IF (LORBREAL) THEN
            !vo : bk -> kb
            VO_S_R(:, :) = T1_R(:, :, KB)
            CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)
         ELSE
            !vo : bk -> kb
            VO_S(:, :) = T1(:, :, RKIofKI(KB))
            CALL SORT_V1O1_O1V1(WDES, VO_S, OV_S, VO_S_R, OV_S_R)
         END IF
         DO KI = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
            DO KA = 1, MY_NKPTS
               KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                        WDES%VKPT(:, KI), KPOINTS_FULL)

               KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) + &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)

               DO NA = 1, (NUNOCC)

                  CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                  FTOD_PW_AI(1, 1, NA, RKQofKQ(KQ), KA, 1), FTOD_OC_AI(1, 1, NA, RKQofKQ(KQ), KA, 1), &
                                  (VBMAX), &
                                  PW_IA_TMP(1, 1, 1, RKQofKQ(KQ)), OC_IA_TMP(1, 1, 1, RKQofKQ(KQ)), &
                                  (NUNOCC)*(VBMAX), OVOV(1, 1, 1, 1), OVOV_R(1, 1, 1, 1), zero)
                  IF (LORBREAL) THEN
                     !oovv : ick_
                     DO NI = 1, VBMAX
                        VO_S_R(:, :) = OVOV_R(NI, :, :, 1) !ck

                        CALL SGEMM('t', 'n', (VBMAX), (VBMAX), &
                                   (NUNOCC), (1._qs), VO_S_R(1, 1), &
                                   (NUNOCC), &
                                   T1_R(1, 1, KJ), (NUNOCC), &
                                   (0._qs), OOOV_S_R(1, 1, NI, 1), (VBMAX))
                        !ooov_s : kji_
                     END DO
                  ELSE
                     !oovv : ick_
                     DO NI = 1, VBMAX
                        VO_S(:, :) = OVOV(NI, :, :, 1) !ck

                        CALL CGEMM('t', 'n', (VBMAX), (VBMAX), &
                                   (NUNOCC), (1._qs, 0._qs), VO_S(1, 1), &
                                   (NUNOCC), &
                                   T1(1, 1, RKIofKI(KJ)), (NUNOCC), &
                                   (0._qs, 0._qs), OOOV_S(1, 1, NI, 1), (VBMAX))
                        !ooov_s : kji_
                     END DO
                  END IF

                  CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                  FTOD_PW_AI(1, 1, NA, RKQofKQ(KQ), KA, 1), FTOD_OC_AI(1, 1, NA, RKQofKQ(KQ), KA, 1), &
                                  (VBMAX), &
                                  PW_IJ_TMP(1, 1, 1, RKQofKQ(KQ)), OC_IJ_TMP(1, 1, 1, RKQofKQ(KQ)), &
                                  (VBMAX*VBMAX), OOOV(1, 1, 1, 1), OOOV_R(1, 1, 1, 1), zero)
                  !ooov : ijk
                  IF (LORBREAL) THEN
                     DO NI = 1, VBMAX
                        DO NJ = 1, VBMAX
                           OOOV_S_R(:, NJ, NI, 1) = -OOOV_S_R(:, NJ, NI, 1) - OOOV_R(NI, NJ, :, 1)*(1._qs)
                           DO NB = 1, (NUNOCC)
                              CALL SGEMM('t', 'n', 1, 1, &
                                         (VBMAX), (1._qs), OOOV_S_R(1, NJ, NI, 1), &
                                         (VBMAX), &
                                         OV_S_R(1, NB), (VBMAX), &
                                         (1._qs), T2_N_R(NA, NB, NI, NJ, KI, KJ, KA), 1)
                           END DO !NB
                        END DO
                     END DO
                  ELSE
                     DO NI = 1, VBMAX
                        DO NJ = 1, VBMAX
                           OOOV_S(:, NJ, NI, 1) = -OOOV_S(:, NJ, NI, 1) - OOOV(NI, NJ, :, 1)*(1._qs, 0._qs)
                           DO NB = 1, (NUNOCC)
                              CALL CGEMM('t', 'n', 1, 1, &
                                         (VBMAX), (1._qs, 0._qs), OOOV_S(1, NJ, NI, 1), &
                                         (VBMAX), &
                                         OV_S(1, NB), (VBMAX), &
                                         (1._qs, 0._qs), T2_N(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA), 1)
                           END DO !NB
                        END DO
                     END DO
                  END IF
               END DO !NA
            END DO
         END DO
      END DO

   END SUBROUTINE FORM_CONTR_VT1_VT1T1_K

!***********************************************************************
!
! The following routine performs one iteration of the
! T1 and T2 amplitude equations by calling different routines that compute and contract
! intermediate quantities.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
! TODO: Clean output and summarize timings better.
!
!***********************************************************************

   Subroutine CALC_T1_T2_N(WDES, WGW, W, ITERATION, KPOINTS)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      TYPE(wavespin) W
      TYPE(kpoints_struct) KPOINTS
      INTEGER ITERATION
      INTEGER KI, KJ, KA, KB, NI, NJ, NA, NB, KQ, KQ_
      integer :: time_array1(8), time_array2(8), ems

      IF (LORBREAL) THEN
         T2_N_R = (0.0_qs)
         IF (SINGLES) THEN
            T1_N_R = (0.0_qs)
            IF (.not. CANONICAL) T1_N_R = (F_AI_R)
         END IF
         CHI_CKAI_R = (0.0_qs)
      ELSE
         T2_N = (0.0_qs, 0.0_qs)
         IF (SINGLES) THEN
            T1_N = (0.0_qs, 0.0_qs)
            IF (.not. CANONICAL) T1_N = (F_AI)
         END IF
         CHI_CKAI = (0.0_qs, 0.0_qs)
      END IF

      TETSORT = 0


      IF ((.NOT. RING) .and. (.NOT. LCCD)) THEN
         call date_and_time(values=time_array1)
         if (ME == 0) WRITE (*, *) 'forming k_ki'
         IF (LORBREAL) THEN
            K_KI_R = (0.0_qs)
         ELSE
            K_KI = (0.0_qs, 0.0_qs)
         END IF
         CALL FORM_K_KI(WGW, WDES)
         if (ME == 0) WRITE (*, *) 'forming l_ki'
         CALL FORM_L_KI(WGW, WDES)
         if (ME == 0) WRITE (*, *) 'contracting l_ki'
         CALL CONTR_L_KI_T2(WDES)

         if (ME == 0) WRITE (*, *) 'forming k_ac'
         CALL FORM_K_AC(WGW, WDES)
         if (ME == 0) WRITE (*, *) 'forming l_ac'
         CALL FORM_L_AC(WGW, WDES)
         if (ME == 0) WRITE (*, *) 'contracting l_ac'
         CALL CONTR_L_AC_T2(WDES)
      END IF

      IF (CCMP2) GO TO 1001

      IF (SINGLES) THEN
         IF (.not. CANONICAL) CALL CONTR_F_KC_T1T1(WDES)
         CALL FORM_K_KC(WGW, WDES)
         if (ME == 0) write (*, *) 'contracting * T1'

         call date_and_time(values=time_array1)
         CALL CONTR_K_AC_T1(WDES)
         call date_and_time(values=time_array2)
         ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                time_array2(8) - time_array1(8))
         if (ME == 0) write (*, *) 'for CONTR_K_AC_T1 ', ems, 'passed'

         call date_and_time(values=time_array1)
         CALL CONTR_K_KI_T1(WDES)
         call date_and_time(values=time_array2)
         ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                time_array2(8) - time_array1(8))
         if (ME == 0) write (*, *) 'for CONTR_K_KI_T1 ', ems, 'passed'

         call date_and_time(values=time_array1)
         CALL CONTR_K_KC_T(WGW, WDES)
         call date_and_time(values=time_array2)
         ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                time_array2(8) - time_array1(8))
         if (ME == 0) write (*, *) 'for CONTR_K_KC_T ', ems, 'passed'

         call date_and_time(values=time_array1)
         CALL CONTR_W_AKIC_T1(WGW, WDES)
         call date_and_time(values=time_array2)
         ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                time_array2(8) - time_array1(8))
         if (ME == 0) write (*, *) 'for CONTR_W_AKIC_T1 ', ems, 'passed'

         call date_and_time(values=time_array1)
         CALL CONTR_W_AKCD_T(WGW, WDES)
         call date_and_time(values=time_array2)
         ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                time_array2(8) - time_array1(8))
         if (ME == 0) write (*, *) 'for CONTR_W_AKCD_T ', ems, 'passed'

         call date_and_time(values=time_array1)
         CALL CONTR_W_KLIC_T(WGW, WDES)
         call date_and_time(values=time_array2)
         ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                time_array2(8) - time_array1(8))
         if (ME == 0) write (*, *) 'for CONTR_W_KLIC_T ', ems, 'passed'

         if (ME == 0) write (*, *) 'contracting * T1 -> T2'
         call date_and_time(values=time_array1)
         CALL FORM_CONTR_VT1_VT1T1_C(WDES)
         call date_and_time(values=time_array2)
         ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                time_array2(8) - time_array1(8))
         if (ME == 0) write (*, *) 'for FORM_CONTR_VT1_VT1T1_C ', ems, 'passed'
         call date_and_time(values=time_array1)
         CALL FORM_CONTR_VT1_VT1T1_K(WDES)
         call date_and_time(values=time_array2)
         ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                time_array2(8) - time_array1(8))
         if (ME == 0) write (*, *) 'for FORM_CONTR_VT1_VT1T1_K ', ems, 'passed'
      END IF

      if (ME == 0) WRITE (*, *) 'forming chi_akic'

      call date_and_time(values=time_array1)
      CALL FORM_CHI_AKIC(WDES)
      call date_and_time(values=time_array2)
      ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
             time_array2(8) - time_array1(8))
      if (ME == 0) write (*, *) 'for FORM_CHI_AKIC ', ems, 'passed'

      If (ME == 0) WRITE (*, *) 'contracting chi_akic'
      call date_and_time(values=time_array1)
      CALL CONTR_CHI_AKIC_AS_T2_BCKJ(WDES)
      call date_and_time(values=time_array2)
      ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
             time_array2(8) - time_array1(8))
      if (ME == 0) write (*, *) 'for CONTR_CHI_AKIC_AS_T2_BCKJ ', ems, 'passed'

      IF ((.NOT. RING) .or. (EH_SCREENING)) THEN
         if (ME == 0) WRITE (*, *) 'forming chi_akci'
         call date_and_time(values=time_array1)
         CALL FORM_CHI_AKCI(WDES)
         call date_and_time(values=time_array2)
         ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                time_array2(8) - time_array1(8))
         if (ME == 0) write (*, *) 'for FORM_CHI_AKCI ', ems, 'passed'

         if (ME == 0) WRITE (*, *) 'contracting chi_akci'

         call date_and_time(values=time_array1)
         IF ((.not. RING)) CALL CONTR_CHI_AKCI_T2_BCKJ(WDES) !ok
         call date_and_time(values=time_array2)
         ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                time_array2(8) - time_array1(8))
         if (ME == 0) write (*, *) 'for CONTR_CHI_AKCI_T2_BCKJ ', ems, 'passed'

         call date_and_time(values=time_array1)
         CALL CONTR_CHI_BKCI_T2_ACKJ(WDES) !ok
         call date_and_time(values=time_array2)
         ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                time_array2(8) - time_array1(8))
         if (ME == 0) write (*, *) 'for CONTR_CHI_BKCI_T2_ACKJ ', ems, 'passed'

      END IF

1001  CONTINUE


      if (ME == 0) WRITE (*, *) 'calculating P{T2}'
      CALL P_T2_N(WDES)

      IF (LORBREAL) THEN
         DO KB = 1, WDES%NKPTS
            CALL BCAST2ALL_FTOD_AI(WDES, KB, 1, PW_AI_TMP, OC_AI_TMP)
            DO KA = 1, MY_NKPTS
            DO KQ = 1, WDES%NKPTS
               KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) + &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)
               KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)
               CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                               PW_AI_TMP(1, 1, 1, KQ), OC_AI_TMP(1, 1, 1, KQ), &
                               (NUNOCC)*VBMAX, &
                               FTOD_PW_AI(1, 1, 1, KQ, KA, 2), FTOD_OC_AI(1, 1, 1, KQ, KA, 2), &
                               (NUNOCC)*VBMAX, OVOV(1, 1, 1, 1), OVOV_R(1, 1, 1, 1), zero)
               CALL SORT_O1V1O2V2_V2V1O2O1(WDES, OVOV, &
                                           VVOO, OVOV_R, VVOO_R)
               T2_N_R(:, :, :, :, KI, KJ, KA) = T2_N_R(:, :, :, :, KI, KJ, KA) + VVOO_R(:, :, :, :)
            END DO
            END DO
         END DO

      ELSE

         DO KB = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0.0_q)) CYCLE
            CALL BCAST2ALL_FTOD_AI(WDES, KB, 1, PW_AI_TMP, OC_AI_TMP)
            DO KA = 1, MY_NKPTS
            DO KQ = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
               KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) + &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)
               KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) - &
                                        WDES%VKPT(:, KQ), KPOINTS_FULL)
               CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                               PW_AI_TMP(1, 1, 1, RKQofKQ(KQ)), OC_AI_TMP(1, 1, 1, RKQofKQ(KQ)), &
                               (NUNOCC)*VBMAX, &
                               FTOD_PW_AI(1, 1, 1, RKQofKQ(KQ), KA, 2), FTOD_OC_AI(1, 1, 1, RKQofKQ(KQ), KA, 2), &
                               (NUNOCC)*VBMAX, OVOV(1, 1, 1, 1), OVOV_R(1, 1, 1, 1), zero)
               CALL SORT_O1V1O2V2_V2V1O2O1(WDES, OVOV, &
                                           VVOO, OVOV_R, VVOO_R)
               T2_N(:, :, :, :, RKIofKI(KI), RKIofKI(KJ), KA) = T2_N(:, :, :, :, RKIofKI(KI), RKIofKI(KJ), KA) + VVOO(:, :, :, :)
            END DO
            END DO
         END DO

      END IF

      IF (CCMP2) GO TO 1002

      IF (.NOT. (RING)) THEN
         TET2S = 0
         call date_and_time(values=time_array1)
         if (ME == 0) WRITE (*, *) 'forming chi_klij'
         CALL FORM_CHI_KLIJ(WDES)
         if (ME == 0) WRITE (*, *) 'contracting chi_klij'
         CALL CONTR_CHI_KLIJ_T(WDES)
         call date_and_time(values=time_array2)
         ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                time_array2(8) - time_array1(8))
         if (ME == 0) WRITE (*, *) 'ems=', ems

         TETS = 0
         call date_and_time(values=time_array1)
         if (ME == 0) WRITE (*, *) 'forming and contracting chi_abcd'
         CALL FORM_CONTR_ABCD(WDES, KPOINTS)
         IF (ME == 0) WRITE (*, *) '{T2_ijab=T2_jiba}'
         CALL UP_T2_N(WDES)
         call date_and_time(values=time_array2)
         ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                time_array2(8) - time_array1(8))
         if (ME == 0) WRITE (*, *) 'ems=', ems
      END IF

1002  CONTINUE

   END Subroutine CALC_T1_T2_N

!***********************************************************************
!
! The following routine performs one iteration of the
! T1 and T2 amplitude equations by accounting for the ppl term only.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   Subroutine CALC_T1_T2_N_PPL(WDES, WGW, W, ITERATION, KPOINTS)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      TYPE(wavespin) W
      TYPE(kpoints_struct) KPOINTS
      INTEGER ITERATION
      INTEGER KI, KJ, KA, KB, NI, NJ, NA, NB, KQ, KQ_
      integer :: time_array1(8), time_array2(8), ems

      IF (LORBREAL) THEN
         T2_N_R = (0.0_qs)
         IF (SINGLES) THEN
            T1_N_R = (0.0_qs)
            IF (.not. CANONICAL) T1_N_R = (F_AI_R)
         END IF
         CHI_CKAI_R = (0.0_qs)
      ELSE
         T2_N = (0.0_qs, 0.0_qs)
         IF (SINGLES) THEN
            T1_N = (0.0_qs, 0.0_qs)
            IF (.not. CANONICAL) T1_N = (F_AI)
         END IF
         CHI_CKAI = (0.0_qs, 0.0_qs)
      END IF

      IF (.NOT. (RING)) THEN
         TETS = 0
         call date_and_time(values=time_array1)
         if (ME == 0) WRITE (*, *) 'forming and contracting chi_abcd'
         CALL FORM_CONTR_ABCD_PPL(WDES, KPOINTS)
         IF (ME == 0) WRITE (*, *) '{T2_ijab=T2_jiba}'
         CALL UP_T2_N(WDES)
         call date_and_time(values=time_array2)
         ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                time_array2(8) - time_array1(8))
         if (ME == 0) WRITE (*, *) 'ems=', ems
      END IF

   END Subroutine CALC_T1_T2_N_PPL

!***********************************************************************
!
! The following routine performs one iteration of the
! T1 and T2 amplitude equations accounting for the hhl term only.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   Subroutine CALC_T1_T2_N_HHL(WDES, WGW, W, ITERATION, KPOINTS)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      TYPE(wavespin) W
      TYPE(kpoints_struct) KPOINTS
      INTEGER ITERATION
      INTEGER KI, KJ, KA, KB, NI, NJ, NA, NB, KQ, KQ_
      integer :: time_array1(8), time_array2(8), ems

      IF (LORBREAL) THEN
         T2_N_R = (0.0_qs)
         IF (SINGLES) THEN
            T1_N_R = (0.0_qs)
            IF (.not. CANONICAL) T1_N_R = (F_AI_R)
         END IF
         CHI_CKAI_R = (0.0_qs)
      ELSE
         T2_N = (0.0_qs, 0.0_qs)
         IF (SINGLES) THEN
            T1_N = (0.0_qs, 0.0_qs)
            IF (.not. CANONICAL) T1_N = (F_AI)
         END IF
         CHI_CKAI = (0.0_qs, 0.0_qs)
      END IF

      TET2S = 0
      call date_and_time(values=time_array1)
      if (ME == 0) WRITE (*, *) 'forming chi_klij'
      CALL FORM_CHI_KLIJ(WDES)
      if (ME == 0) WRITE (*, *) 'contracting chi_klij'
      CALL CONTR_CHI_KLIJ_T(WDES)
      call date_and_time(values=time_array2)
      ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
             time_array2(8) - time_array1(8))
      if (ME == 0) WRITE (*, *) 'ems=', ems

   END Subroutine CALC_T1_T2_N_HHL

!***********************************************************************
!
! The following routine performs one iteration of the
! T1 and T2 amplitude equations accounting for the phl term only.
! We employ a notation that is similar to the one used in
!       Coupled-cluster singles and doubles for extended systems
!       So Hirata, Rafał Podeszwa, Motoi Tobita, and Rodney J. Bartlett
!       J. Chem. Phys. 120, 2581 (2004)
!
!***********************************************************************

   Subroutine CALC_T1_T2_N_PHL(WDES, WGW, W, ITERATION, KPOINTS)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      TYPE(wavespin) W
      TYPE(kpoints_struct) KPOINTS
      INTEGER ITERATION
      INTEGER KI, KJ, KA, KB, NI, NJ, NA, NB, KQ, KQ_
      integer :: time_array1(8), time_array2(8), ems

      IF (LORBREAL) THEN
         T2_N_R = (0.0_qs)
         IF (SINGLES) THEN
            T1_N_R = (0.0_qs)
            IF (.not. CANONICAL) T1_N_R = (F_AI_R)
         END IF
         CHI_CKAI_R = (0.0_qs)
      ELSE
         T2_N = (0.0_qs, 0.0_qs)
         IF (SINGLES) THEN
            T1_N = (0.0_qs, 0.0_qs)
            IF (.not. CANONICAL) T1_N = (F_AI)
         END IF
         CHI_CKAI = (0.0_qs, 0.0_qs)
      END IF

      IF ((.NOT. RING) .or. (EH_SCREENING)) THEN
         if (ME == 0) WRITE (*, *) 'forming chi_akci'
         call date_and_time(values=time_array1)
         CALL FORM_CHI_AKCI(WDES)
         call date_and_time(values=time_array2)
         ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                time_array2(8) - time_array1(8))
         if (ME == 0) write (*, *) 'for FORM_CHI_AKCI ', ems, 'passed'
!
         if (ME == 0) WRITE (*, *) 'contracting chi_akci'

         call date_and_time(values=time_array1)
         IF ((.not. RING)) CALL CONTR_CHI_AKCI_T2_BCKJ(WDES) !ok
         call date_and_time(values=time_array2)
         ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                time_array2(8) - time_array1(8))
         if (ME == 0) write (*, *) 'for CONTR_CHI_AKCI_T2_BCKJ ', ems, 'passed'

         call date_and_time(values=time_array1)
         CALL CONTR_CHI_BKCI_T2_ACKJ(WDES) !ok
         call date_and_time(values=time_array2)
         ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                time_array2(8) - time_array1(8))
         if (ME == 0) write (*, *) 'for CONTR_CHI_BKCI_T2_ACKJ ', ems, 'passed'

         call date_and_time(values=time_array2)
         ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                time_array2(8) - time_array1(8))
      END IF
      if (ME == 0) WRITE (*, *) 'calculating P{T2}'
      CALL P_T2_N(WDES) !ok

   END Subroutine CALC_T1_T2_N_PHL

!***********************************************************************
!
!
!***********************************************************************

   SUBROUTINE SETUP_KINDEX(WDES)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      INTEGER COUNTER, KI, KJ, KA
      TYPE(wavedes) WDES

      ALLOCATE (KINDEX(WDES%NKPTS, WDES%NKPTS, MY_NKPTS))
      KINDEX = 0
      COUNTER = 0
      DO KA = 1, MY_NKPTS
      DO KI = 1, WDES%NKPTS
      DO KJ = 1, WDES%NKPTS
         COUNTER = COUNTER + 1
         KINDEX(KI, KJ, KA) = COUNTER
      END DO
      END DO
      END DO

   END SUBROUTINE SETUP_KINDEX

!***********************************************************************
!
! This routine reads in the T2CAR file
!
!***********************************************************************

   SUBROUTINE IN_T2(IO, WDES, W, START_IT)
      USE base
      USE wave
      USE lattice
      USE main_mpi
      IMPLICIT NONE
      TYPE(in_struct) IO
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      ! local
      INTEGER NPL_TOT, IU0, IRECLW, MRECL, IREC
      COMPLEX(qs) :: CNUNOCC, CVBMAX, CNKPTS
      REAL(qs) :: RNUNOCC, RVBMAX, RNKPTS
      REAL(q) :: MVKPT
      GDEFS, ALLOCATABLE :: T2_MTMP(:, :, :, :, :)
      REAL(qs), ALLOCATABLE :: T2_MTMP_R(:, :, :, :, :)
      INTEGER :: NI, NJ, NA, NB, KI, KA, KJ, START_IT
      LOGICAL :: file_exists

      INQUIRE (FILE="T2CAR", EXIST=file_exists)
      IF (file_exists) THEN

         IF (ME == 0) WRITE (*, *) 'Found T2CAR file. Reading in T2 amplitudes.'

         START_IT = 2
         CNUNOCC = NUNOCC
         RNUNOCC = NUNOCC
         CVBMAX = VBMAX
         RVBMAX = VBMAX
         CNKPTS = REALNKPTS !WDES%NKPTS
         RNKPTS = REALNKPTS !WDES%NKPTS

         call BLACS_PINFO(ME, PROCS)
         ! RECL : the T2 amplitudes are stored in single precision complex
         MRECL = 8
         IF (ME == 0) THEN

            OPEN (UNIT=12, FILE=DIR_APP(1:DIR_LEN)//'T2CAR', ACCESS='SEQUENTIAL', STATUS='UNKNOWN', FORM='UNFORMATTED')

            !write info about number of orbitals and k-points to header of T2CAR
            READ (12) CNUNOCC
            NUNOCC = CNUNOCC
            READ (12) CVBMAX
            VBMAX = CVBMAX
            READ (12) CNKPTS
            IREC = 4
            IF (REALNKPTS /= CNKPTS) THEN
               WRITE (*, *) 'Number of k-points in T2CAR and KPOINTS files differ.'
               CALL EXIT
            END IF
            DO KA = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA) == 0)) CYCLE
               READ (12) MVKPT
               WDES%VKPT(1, KA) = MVKPT
               IREC = IREC + 1
               READ (12) MVKPT
               WDES%VKPT(2, KA) = MVKPT
               IREC = IREC + 1
               READ (12) MVKPT
               WDES%VKPT(3, KA) = MVKPT
               IREC = IREC + 1

               IF (LMETAL) THEN
               DO NA = 1, NUNOCC
                  READ (12) MVKPT
                  W%CELTOT(NA, KA, 1) = MVKPT*(1.0_q, 0.0_q)
                  IREC = IREC + 1
               END DO
               ELSE
               DO NA = 1, VBMAX + NUNOCC
                  READ (12) MVKPT
                  W%CELTOT(NA, KA, 1) = MVKPT*(1.0_q, 0.0_q)
                  IREC = IREC + 1
               END DO
               END IF

            END DO
         END IF

         ALLOCATE (T2_MTMP(NUNOCC, NUNOCC, VBMAX, VBMAX, REALNKPTS))
         ALLOCATE (T2_MTMP_R(NUNOCC, NUNOCC, VBMAX, VBMAX, WDES%NKPTS))

         !now communicate and write the T2 amplitudes
         DO KA = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA) == 0)) CYCLE
            DO KJ = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE

               IF (ME == 0) THEN
                  DO KI = 1, WDES%NKPTS
                     IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
                     DO NJ = 1, VBMAX
                     DO NI = 1, VBMAX
                     DO NB = 1, NUNOCC
                     DO NA = 1, NUNOCC
                        READ (12) T2_MTMP(NA, NB, NI, NJ, RKIofKI(KI))
                        IREC = IREC + 1
                     END DO
                     END DO
                     END DO
                     END DO
                  END DO
               END IF

               IF (LORBREAL) T2_MTMP_R = REAL(T2_MTMP, kind=qs)

               IF (LORBREAL) THEN
                  IF (ME == 0) THEN
                     CALL SGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC), (VBMAX*VBMAX*WDES%NKPTS), &
                                  T2_MTMP_R(1, 1, 1, 1, 1), (NUNOCC*NUNOCC))
                     IF (PROCS_KPTS(KA) == ME) T2_R(:, :, :, :, :, KJ, MKPTS_KPTS(KA)) = T2_MTMP_R
                  ELSE
                     CALL SGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC), (VBMAX*VBMAX*WDES%NKPTS), &
                                  T2_MTMP_R(1, 1, 1, 1, 1), (NUNOCC*NUNOCC), 0, 0)
                     IF (PROCS_KPTS(KA) == ME) T2_R(:, :, :, :, :, KJ, MKPTS_KPTS(KA)) = T2_MTMP_R
                  END IF
               ELSE
                  IF (ME == 0) THEN
                     CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC), (VBMAX*VBMAX*REALNKPTS), &
                                  T2_MTMP(1, 1, 1, 1, 1), (NUNOCC*NUNOCC))
                     IF (PROCS_KPTS(KA) == ME) T2(:, :, :, :, :, RKIofKI(KJ), MKPTS_KPTS(KA)) = T2_MTMP(:, :, :, :, :)
                  ELSE
                     CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC), (VBMAX*VBMAX*REALNKPTS), &
                                  T2_MTMP(1, 1, 1, 1, 1), (NUNOCC*NUNOCC), 0, 0)
                     IF (PROCS_KPTS(KA) == ME) T2(:, :, :, :, :, RKIofKI(KJ), MKPTS_KPTS(KA)) = T2_MTMP(:, :, :, :, :)
                  END IF
               END IF

            END DO
         END DO

         CLOSE (12)

      END IF

   END SUBROUTINE IN_T2


   SUBROUTINE IN_EMP2_PAIR_CBS(IO, WDES, W, START_IT)
      USE base
      USE wave
      USE lattice
      USE main_mpi
      IMPLICIT NONE
      TYPE(in_struct) IO
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      ! local
      REAL(q) :: E_TEST
      GDEFS, ALLOCATABLE :: T2_MTMP(:, :, :, :, :)
      REAL(qs), ALLOCATABLE :: T2_MTMP_R(:, :, :, :, :)
      INTEGER :: NBI, NBJ, KI, KJ, START_IT, RKI, RKJ
      LOGICAL :: file_exists

      INQUIRE (FILE="Mp2PairEnergies.dat", EXIST=file_exists)
      IF (file_exists) THEN

         IF (ME == 0) WRITE (*, *) 'Found Mp2PairEnergies.dat file. Reading in MP2 pair energies.'

         START_IT = 1
         E_TEST=zero
         call BLACS_PINFO(ME, PROCS)
         ! RECL : the T2 amplitudes are stored in single precision complex
         IF (ME == 0) THEN

            OPEN(unit=12,FILE=DIR_APP(1:DIR_LEN)//'Mp2PairEnergies.dat',FORM='FORMATTED',access='stream',STATUS='UNKNOWN')
!            OPEN (UNIT=12, FILE=DIR_APP(1:DIR_LEN)//'Mp2PairEnergies.dat', ACCESS='SEQUENTIAL', STATUS='UNKNOWN', FORM='UNFORMATTED')

            !write info about number of orbitals and k-points to header of T2CAR
            DO KI=1,WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI)==0)) CYCLE
               
               DO KJ=1,WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ)==0)) CYCLE
               
                  DO NBI=1,VBMAX !loop over valence bands i                     
!                     IF ((LFREEZE) .and. (NBI<=NFREEZE)) CYCLE
                  DO NBJ=1,VBMAX !loop over valence bands i                     
!                     IF ((LFREEZE) .and. (NBJ<=NFREEZE)) CYCLE

                     READ(12,*) EMP2_PAIR_CBS(NBI,NBJ,RKIofKI(KI),RKIofKI(KJ))
                     IF (ME==0) E_TEST=E_TEST+EMP2_PAIR_CBS(NBI,NBJ,RKIofKI(KI),RKIofKI(KJ))


                  ENDDO
                  ENDDO
               ENDDO
            ENDDO

            CLOSE(12)
            WRITE(*,*)'E_MP2 (CBS) =',E_TEST
         ENDIF


         IF (ME == 0) THEN
            CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', VBMAX*VBMAX, (REALNKPTS*REALNKPTS), &
                         EMP2_PAIR_CBS(1, 1, 1, 1), VBMAX*VBMAX)
!            IF (PROCS_KPTS(KA) == ME) T2(:, :, :, :, :, RKIofKI(KJ), MKPTS_KPTS(KA)) = T2_MTMP(:, :, :, :, :)
         ELSE
            CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', VBMAX*VBMAX, (REALNKPTS*REALNKPTS), &
                         EMP2_PAIR_CBS(1, 1, 1, 1), VBMAX*VBMAX, 0, 0)
!            IF (PROCS_KPTS(KA) == ME) T2(:, :, :, :, :, RKIofKI(KJ), MKPTS_KPTS(KA)) = T2_MTMP(:, :, :, :, :)
         END IF


      END IF

   END SUBROUTINE IN_EMP2_PAIR_CBS

!***********************************************************************
!
! This routine reads in the T2CAR file
!
!***********************************************************************

   SUBROUTINE IN_T1(IO, WDES, W, START_IT)
      USE base
      USE wave
      USE lattice
      USE main_mpi
      IMPLICIT NONE
      TYPE(in_struct) IO
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      ! local
      INTEGER NPL_TOT, IU0, IRECLW, MRECL, IREC
      COMPLEX(qs) :: CNUNOCC, CVBMAX, CNKPTS
      REAL(qs) :: RNUNOCC, RVBMAX, RNKPTS
      REAL(q) :: MVKPT
      GDEFS, ALLOCATABLE :: T1_MTMP(:, :, :)
      REAL(qs), ALLOCATABLE :: T1_MTMP_R(:, :, :)
      INTEGER :: NI, NJ, NA, NB, KI, KA, KJ, START_IT
      LOGICAL :: file_exists

      INQUIRE (FILE="T1CAR", EXIST=file_exists)
      IF (file_exists) THEN

         IF (ME == 0) WRITE (*, *) 'Found T1CAR file. Reading in T1 amplitudes.'

         START_IT = 2
         CNUNOCC = NUNOCC
         RNUNOCC = NUNOCC
         CVBMAX = VBMAX
         RVBMAX = VBMAX
         CNKPTS = REALNKPTS
         RNKPTS = WDES%NKPTS

         call BLACS_PINFO(ME, PROCS)
         ! RECL : the T2 amplitudes are stored in single precision complex
         MRECL = 8
         IF (ME == 0) THEN
            OPEN (UNIT=12, FILE=DIR_APP(1:DIR_LEN)//'T1CAR', ACCESS='DIRECT', &
                  FORM='UNFORMATTED', STATUS='UNKNOWN', RECL=MRECL)

            !write info about number of orbitals and k-points to header of T2CAR
            READ (12, REC=1) CNUNOCC
            NUNOCC = CNUNOCC
            READ (12, REC=2) CVBMAX
            VBMAX = CVBMAX
            READ (12, REC=3) CNKPTS
            IREC = 4
            IF (REALNKPTS /= CNKPTS) THEN
               WRITE (*, *) 'Number of k-points in T2CAR and KPOINTS files differ.'
               CALL EXIT
            END IF
            DO KA = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA) == 0)) CYCLE
               READ (12, REC=IREC) MVKPT
               WDES%VKPT(1, KA) = MVKPT
               IREC = IREC + 1
               READ (12, REC=IREC) MVKPT
               WDES%VKPT(2, KA) = MVKPT
               IREC = IREC + 1
               READ (12, REC=IREC) MVKPT
               WDES%VKPT(3, KA) = MVKPT
               IREC = IREC + 1
            END DO
         END IF

         ALLOCATE (T1_MTMP(NUNOCC, VBMAX, REALNKPTS))
         ALLOCATE (T1_MTMP_R(NUNOCC, VBMAX, REALNKPTS))

         !now communicate and write the T2 amplitudes

         IF (ME == 0) THEN
            DO KI = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
               DO NI = 1, VBMAX
               DO NA = 1, NUNOCC
                  READ (12, REC=IREC) T1_MTMP(NA, NI, RKIofKI(KI))
                  IREC = IREC + 1
               END DO
               END DO
            END DO
         END IF

         IF (LORBREAL) T1_MTMP_R = REAL(T1_MTMP, kind=qs)

         IF (LORBREAL) THEN
            IF (ME == 0) THEN
               CALL SGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (VBMAX*REALNKPTS), &
                            T1_MTMP_R(1, 1, 1), (NUNOCC))
               T1_R = T1_MTMP_R
            ELSE
               CALL SGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (VBMAX*REALNKPTS), &
                            T1_MTMP_R(1, 1, 1), (NUNOCC), 0, 0)
               T1_R = T1_MTMP_R
            END IF
         ELSE
            IF (ME == 0) THEN
               CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (VBMAX*REALNKPTS), &
                            T1_MTMP(1, 1, 1), (NUNOCC))
               T1 = T1_MTMP
            ELSE
               CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC), (VBMAX*REALNKPTS), &
                            T1_MTMP(1, 1, 1), (NUNOCC), 0, 0)
               T1 = T1_MTMP
            END IF
         END IF

         CLOSE (12)

      END IF

   END SUBROUTINE IN_T1

!***********************************************************************
!
! This routine writes out the T2 amplitudes to T2CAR
!
!***********************************************************************

   SUBROUTINE OUT_T2(IO, WDES, W)
      USE base
      USE wave
      USE lattice
      USE main_mpi
      IMPLICIT NONE
      TYPE(in_struct) IO
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      ! local
      INTEGER NPL_TOT, IU0, IRECLW, MRECL, IREC
      COMPLEX(qs) :: CNUNOCC, CVBMAX, CNKPTS
      REAL(qs) :: RNUNOCC, RVBMAX, RNKPTS
      GDEFS, ALLOCATABLE :: T2_MTMP(:, :, :, :, :)
      REAL(qs), ALLOCATABLE :: T2_MTMP_R(:, :, :, :, :)
      INTEGER :: NI, NJ, NA, NB, KI, KA, KJ

      CNUNOCC = NUNOCC
      RNUNOCC = NUNOCC
      CVBMAX = VBMAX
      RVBMAX = VBMAX
      CNKPTS = REALNKPTS
      RNKPTS = REALNKPTS

      call BLACS_PINFO(ME, PROCS)
      ! RECL : the T2 amplitudes are stored in single precision complex
      MRECL = 8
      IF (ME == 0) THEN
         OPEN (UNIT=12, FILE=DIR_APP(1:DIR_LEN)//'T2CAR', &
               FORM='UNFORMATTED', STATUS='UNKNOWN')

         !write info about number of orbitals and k-points to header of T2CAR
         WRITE (12) CNUNOCC
         WRITE (12) CVBMAX
         WRITE (12) CNKPTS
         IREC = 4
         DO KA = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA) == 0)) CYCLE
            WRITE (12) WDES%VKPT(1, KA)
            WRITE (12) WDES%VKPT(2, KA)
            WRITE (12) WDES%VKPT(3, KA)
            IREC = IREC + 3

            IF (LMETAL) THEN
            DO NA = 1, NUNOCC
               WRITE (12) REAL(W%CELTOT(NA, KA, 1), kind=q)
               IREC = IREC + 1
            END DO
            ELSE
            DO NA = 1, VBMAX + NUNOCC
               WRITE (12) REAL(W%CELTOT(NA, KA, 1), kind=q)
               IREC = IREC + 1
            END DO
            END IF
         END DO
      END IF

      ALLOCATE (T2_MTMP(NUNOCC, NUNOCC, VBMAX, VBMAX, REALNKPTS))
      ALLOCATE (T2_MTMP_R(NUNOCC, NUNOCC, VBMAX, VBMAX, REALNKPTS))

      !now communicate and write the T2 amplitudes
      DO KA = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA) == 0)) CYCLE
         DO KJ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE

            IF (LORBREAL) THEN
               IF (PROCS_KPTS(KA) == ME) THEN
                  CALL SGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC), (VBMAX*VBMAX*WDES%NKPTS), &
                               T2_R(1, 1, 1, 1, 1, KJ, MKPTS_KPTS(KA)), (NUNOCC*NUNOCC))
                  T2_MTMP_R(:, :, :, :, :) = T2_R(:, :, :, :, :, KJ, MKPTS_KPTS(KA))
               ELSE
                  CALL SGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC), (VBMAX*VBMAX*WDES%NKPTS), &
                               T2_MTMP_R(1, 1, 1, 1, 1), (NUNOCC*NUNOCC), 0, PROCS_KPTS(KA))
               END IF
            ELSE
               IF (PROCS_KPTS(KA) == ME) THEN
                  CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC), (VBMAX*VBMAX*REALNKPTS), &
                               T2(1, 1, 1, 1, 1, RKIofKI(KJ), MKPTS_KPTS(KA)), (NUNOCC*NUNOCC))
                  T2_MTMP(:, :, :, :, :) = T2(:, :, :, :, :, RKIofKI(KJ), MKPTS_KPTS(KA))
               ELSE
                  CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NUNOCC)*(NUNOCC), (VBMAX*VBMAX*REALNKPTS), &
                               T2_MTMP(1, 1, 1, 1, 1), (NUNOCC*NUNOCC), 0, PROCS_KPTS(KA))
               END IF
            END IF

            IF (LORBREAL) T2_MTMP = T2_MTMP_R
            IF (ME == 0) THEN
               DO KI = 1, WDES%NKPTS
                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
                  DO NJ = 1, VBMAX
                  DO NI = 1, VBMAX
                  DO NB = 1, NUNOCC
                  DO NA = 1, NUNOCC
                     WRITE (12) T2_MTMP(NA, NB, NI, NJ, RKIofKI(KI))
                     IREC = IREC + 1
                  END DO
                  END DO
                  END DO
                  END DO
               END DO
            END IF

         END DO
      END DO

      CLOSE (12)

   END SUBROUTINE OUT_T2


!***********************************************************************
!
! This routine writes out the T1 amplitudes to T1CAR
!
!***********************************************************************

   SUBROUTINE OUT_T1(IO, WDES, W)
      USE base
      USE wave
      USE lattice
      USE main_mpi
      IMPLICIT NONE
      TYPE(in_struct) IO
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      ! local
      INTEGER NPL_TOT, IU0, IRECLW, MRECL, IREC
      COMPLEX(qs) :: CNUNOCC, CVBMAX, CNKPTS
      REAL(qs) :: RNUNOCC, RVBMAX, RNKPTS
      GDEFS, ALLOCATABLE :: T1_MTMP(:, :, :)
      INTEGER :: NI, NJ, NA, NB, KI, KA, KJ

      CNUNOCC = NUNOCC
      RNUNOCC = NUNOCC
      CVBMAX = VBMAX
      RVBMAX = VBMAX
      CNKPTS = REALNKPTS
      RNKPTS = WDES%NKPTS

      call BLACS_PINFO(ME, PROCS)
      ! RECL : the T2 amplitudes are stored in single precision complex
      MRECL = 8
      IF (ME == 0) THEN
         OPEN (UNIT=12, FILE=DIR_APP(1:DIR_LEN)//'T1CAR', ACCESS='DIRECT', &
               FORM='UNFORMATTED', STATUS='UNKNOWN', RECL=MRECL)

         !write info about number of orbitals and k-points to header of T2CAR
         WRITE (12, REC=1) CNUNOCC
         WRITE (12, REC=2) CVBMAX
         WRITE (12, REC=3) CNKPTS
         IREC = 4
         DO KA = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA) == 0)) CYCLE
            WRITE (12, REC=IREC) WDES%VKPT(1, KA)
            WRITE (12, REC=IREC + 1) WDES%VKPT(2, KA)
            WRITE (12, REC=IREC + 2) WDES%VKPT(3, KA)
            IREC = IREC + 3
         END DO
      END IF

      ALLOCATE (T1_MTMP(NUNOCC, VBMAX, REALNKPTS))

      !now communicate and write the T2 amplitudes
      IF (LORBREAL) THEN
         T1_MTMP = T1_R
      ELSE
         T1_MTMP = T1
      END IF
      IF (ME == 0) THEN
         DO KI = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE
            DO NI = 1, VBMAX
            DO NA = 1, NUNOCC
               WRITE (12, REC=IREC) T1_MTMP(NA, NI, RKIofKI(KI))
               IREC = IREC + 1
            END DO
            END DO
         END DO
      END IF

      CLOSE (12)

   END SUBROUTINE OUT_T1

!***********************************************************************
!
! This routine writes out the Coulomb-Vertices to FTODCAR that are needed for the calculation of the (T)-energies
!
!***********************************************************************

   SUBROUTINE OUT_FTOD_PW(IO, WDES, W)
      USE base
      USE wave
      USE lattice
      USE main_mpi
      IMPLICIT NONE
      TYPE(in_struct) IO
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      ! local
      INTEGER NPL_TOT, IU0, IRECLW, MRECL, IREC
      COMPLEX(qs) :: CNUNOCC, CVBMAX, CNKPTS, CNGVECTOR
      REAL(qs) :: RNUNOCC, RVBMAX, RNKPTS
      INTEGER :: NI, NJ, NA, NB, KI, KA, KJ, NG, KQ, mncc

      CNUNOCC = NUNOCC
      RNUNOCC = NUNOCC
      CVBMAX = VBMAX
      RVBMAX = VBMAX
      CNKPTS = REALNKPTS
      RNKPTS = REALNKPTS
      CNGVECTOR = NGVECTOR

      call BLACS_PINFO(ME, PROCS)
      ! RECL : the FTODs are stored in double precision complex
      MRECL = 16
      IF (ME == 0) THEN
         OPEN (UNIT=12, FILE=DIR_APP(1:DIR_LEN)//'FTODCAR', ACCESS='DIRECT', &
               FORM='UNFORMATTED', STATUS='UNKNOWN', RECL=MRECL)

         !write info about number of orbitals and k-points to header of T2CAR
         WRITE (12, REC=1) CNUNOCC
         WRITE (12, REC=2) CVBMAX
         WRITE (12, REC=3) CNKPTS

         WRITE (12, REC=4) CNGVECTOR

         IREC = 5
         DO KA = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA) == 0)) CYCLE
            WRITE (12, REC=IREC) WDES%VKPT(1, KA)
            WRITE (12, REC=IREC + 1) WDES%VKPT(2, KA)
            WRITE (12, REC=IREC + 2) WDES%VKPT(3, KA)
            IREC = IREC + 3
         END DO
      END IF

      !now communicate and write the FTOD_PW_AB amplitudes
      DO KA = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA) == 0)) CYCLE

         IF (PROCS_KPTS(KA) == ME) THEN
            CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR)*(NUNOCC), (NUNOCC*REALNKPTS), &
                         FTOD_PW_AB(1, 1, 1, 1, MKPTS_KPTS(KA), 1), (NGVECTOR*NUNOCC))
            PW_AB_TMP(:, :, :, :) = FTOD_PW_AB(:, :, :, :, MKPTS_KPTS(KA), 1)
         ELSE
            CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR)*(NUNOCC), (NUNOCC*REALNKPTS), &
                         PW_AB_TMP(1, 1, 1, 1), (NGVECTOR*NUNOCC), 0, PROCS_KPTS(KA))
         END IF

         IF (ME == 0) THEN
            DO KQ = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
               DO NB = 1, NUNOCC
               DO NA = 1, NUNOCC
               DO NG = 1, NGVECTOR
                  WRITE (12, REC=IREC) PW_AB_TMP(NG, NA, NB, RKQofKQ(KQ))
                  IREC = IREC + 1
               END DO
               END DO
               END DO
            END DO
         END IF

      END DO

      !now communicate and write the FTOD_PW_AI amplitudes
      DO mncc = 1, 2
      DO KA = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA) == 0)) CYCLE

         IF (PROCS_KPTS(KA) == ME) THEN
            CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR)*(VBMAX), (NUNOCC*WDES%NKPTS), &
                         FTOD_PW_AI(1, 1, 1, 1, MKPTS_KPTS(KA), mncc), (NGVECTOR*VBMAX))
            PW_AI_TMP(:, :, :, :) = FTOD_PW_AI(:, :, :, :, MKPTS_KPTS(KA), mncc)
         ELSE
            CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR)*(VBMAX), (NUNOCC*WDES%NKPTS), &
                         PW_AI_TMP(1, 1, 1, 1), (NGVECTOR*VBMAX), 0, PROCS_KPTS(KA))
         END IF

         IF (ME == 0) THEN
            DO KQ = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
               DO NA = 1, NUNOCC
               DO NI = 1, VBMAX
               DO NG = 1, NGVECTOR
                  WRITE (12, REC=IREC) PW_AI_TMP(NG, NI, NA, RKQofKQ(KQ))
                  IREC = IREC + 1
               END DO
               END DO
               END DO
            END DO
         END IF

      END DO
      END DO

      !now communicate and write the FTOD_PW_IJ amplitudes
      DO mncc = 1, 1
      DO KI = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE

         IF (PROCS_KPTS(KI) == ME) THEN
            CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR)*(VBMAX), (VBMAX*REALNKPTS), &
                         FTOD_PW_IJ(1, 1, 1, 1, MKPTS_KPTS(KI), mncc), (NGVECTOR*VBMAX))
            PW_IJ_TMP(:, :, :, :) = FTOD_PW_IJ(:, :, :, :, MKPTS_KPTS(KI), mncc)
         ELSE
            CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NGVECTOR)*(VBMAX), (VBMAX*REALNKPTS), &
                         PW_IJ_TMP(1, 1, 1, 1), (NGVECTOR*VBMAX), 0, PROCS_KPTS(KI))
         END IF

         IF (ME == 0) THEN
            DO KQ = 1, WDES%NKPTS
               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
               DO NI = 1, VBMAX
               DO NJ = 1, VBMAX
               DO NG = 1, NGVECTOR
                  WRITE (12, REC=IREC) PW_IJ_TMP(NG, NJ, NI, RKQofKQ(KQ))
                  IREC = IREC + 1
               END DO
               END DO
               END DO
            END DO
         END IF

      END DO
      END DO

      CLOSE (12)

   END SUBROUTINE OUT_FTOD_PW

!***********************************************************************
!
! This is an experimental routine that computes the (T) contribution to the CCSD(T) energy.
!
!***********************************************************************

   SUBROUTINE CALC_TRIPLES(WDES, WGW, W, ITERATION, IO, KPOINTS)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      TYPE(wavespin) W
      TYPE(in_struct) IO
      TYPE(kpoints_struct) KPOINTS
      INTEGER ITERATION, NBLOCK
      INTEGER KI, KJ, KK, KA, KB, KC, NI, NJ, NK, NA, NB, NC, KQ, KQ_, KMAX, KTMP, MI
      INTEGER RNBLOCKA, RNBLOCKB, RNBLOCKC, RNA, RNB, RNC, MNA, MNB, MNC
      integer :: time_array1(8), time_array2(8), emstriples
      COMPLEX(q) :: KETRIPLES(WDES%NKPTS)
      LOGICAL :: file_exists
      CHARACTER(40) :: filename

      KETRIPLES = zero
      call date_and_time(values=time_array1)

      WRITE (*, *) 'Calculating (T) contribution to correlation energy'

      WRITE (*, *) 'setting up k-point indexing'
      CALL SETUP_KINDEX(WDES)

      NBLOCK = MIN(15, NUNOCC)
      ETRIPLES = (0.0_q, 0.0_q)

      WRITE (*, *) 'Allocating blocked W_ijk^abc intermediates'
      ALLOCATE (W_ijk_abc(NBLOCK, NBLOCK, NBLOCK, VBMAX, VBMAX, VBMAX, WDES%NKPTS*WDES%NKPTS*MY_NKPTS))
      ALLOCATE (Wint_ijk_abc(NBLOCK, NBLOCK, NBLOCK, VBMAX, VBMAX, VBMAX, WDES%NKPTS*WDES%NKPTS*MY_NKPTS))
      ALLOCATE (Wtmp_ijk_abc(NBLOCK, NBLOCK, NBLOCK, VBMAX, VBMAX, VBMAX, WDES%NKPTS*WDES%NKPTS*MY_NKPTS))
      KMAX = WDES%NKPTS*WDES%NKPTS*MY_NKPTS

      WRITE (*, *) 'Calculating blocked W_ijk^abc', KMAX
      DO KA = 1, WDES%NKPTS
         DO KB = 1, WDES%NKPTS
            write (*, *) 'ka,kb', KA, KB
            MI = 1
            DO NA = 1, (NUNOCC), NBLOCK
               RNBLOCKA = MIN((NUNOCC) - NA + 1, NBLOCK)
               DO NB = 1, (NUNOCC), NBLOCK
                  RNBLOCKB = MIN((NUNOCC) - NB + 1, NBLOCK)
                  DO NC = 1, (NUNOCC), NBLOCK
                     RNBLOCKC = MIN((NUNOCC) - NC + 1, NBLOCK)

                     write (filename, *) MI
                     OPEN (UNIT=12, FILE='BLOCK.'//trim(adjustl(filename))//".dat", &
                           FORM='UNFORMATTED', STATUS='UNKNOWN')
                     MI = MI + 1

                     write (12) NA, NB, NC
                     CLOSE (12)
                  END DO
               END DO
            END DO

            DO MNA = 1, (NUNOCC), NBLOCK
            DO MNB = 1, (NUNOCC), NBLOCK
            DO MNC = 1, (NUNOCC), NBLOCK

               NA = MNA
               NB = MNB
               NC = MNC

               INQUIRE (FILE="BLOCK", EXIST=file_exists)
               IF (file_exists) THEN
                  OPEN (UNIT=12, FILE='BLOCK', &
                        FORM='UNFORMATTED', STATUS='UNKNOWN')
                  READ (12) NA, NB, NC
                  CLOSE (12)
               END IF
               RNBLOCKA = MIN((NUNOCC) - NA + 1, NBLOCK)
               RNBLOCKB = MIN((NUNOCC) - NB + 1, NBLOCK)
               RNBLOCKC = MIN((NUNOCC) - NC + 1, NBLOCK)
               write (*, *) 'running NA;NB;NC: ', NA, NB, NC

               DEALLOCATE (W_ijk_abc, Wint_ijk_abc, Wtmp_ijk_abc)
               ALLOCATE (W_ijk_abc(RNBLOCKA, RNBLOCKB, RNBLOCKC, VBMAX, VBMAX, VBMAX, WDES%NKPTS*WDES%NKPTS*MY_NKPTS))
               ALLOCATE (Wint_ijk_abc(RNBLOCKA, RNBLOCKB, RNBLOCKC, VBMAX, VBMAX, VBMAX, WDES%NKPTS*WDES%NKPTS*MY_NKPTS))
               ALLOCATE (Wtmp_ijk_abc(RNBLOCKA, RNBLOCKB, RNBLOCKC, VBMAX, VBMAX, VBMAX, WDES%NKPTS*WDES%NKPTS*MY_NKPTS))

               W_ijk_abc = (0.0_qs, 0.0_qs)
               Wint_ijk_abc = (0.0_qs, 0.0_qs)
               Wtmp_ijk_abc = (0.0_qs, 0.0_qs)
               Wint_ijk_abc = (0.0_qs, 0.0_qs)

               DO KI = 1, WDES%NKPTS
                  EMSCALCWINT = 0
                  DO KJ = 1, WDES%NKPTS
                  DO KK = 1, MY_NKPTS

                     KC = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KI)) + &
                                WDES%VKPT(:, KJ) + WDES%VKPT(:, KPTS_MKPTS(KK)) - WDES%VKPT(:, KA) - WDES%VKPT(:, KB), KPOINTS_FULL)

                     DEALLOCATE (W_ijk_abc)
                     ALLOCATE (W_ijk_abc(RNBLOCKB, RNBLOCKA, RNBLOCKC, VBMAX, VBMAX, VBMAX, WDES%NKPTS*WDES%NKPTS*MY_NKPTS))
                     W_ijk_abc(:, :, :, :, :, :, KMAX) = (0.0_qs, 0.0_qs)
! ijk_abc=ijk_abc+ jik_bac
                     CALL CALC_WINT(KJ, KI, KPTS_MKPTS(KK), KB, KA, KC, NB, NA, NC, RNBLOCKB, RNBLOCKA, RNBLOCKC, WDES, KPOINTS)
                     DO NI = 1, VBMAX
                     DO NJ = 1, VBMAX
                        DO RNA = 1, RNBLOCKA
                        DO RNB = 1, RNBLOCKB

                           Wint_ijk_abc(RNA,RNB,1:RNBLOCKC,NI,NJ,:,KINDEX(KI,KJ,KK))=Wint_ijk_abc(RNA,RNB,1:RNBLOCKC,NI,NJ,:,KINDEX(KI,KJ,KK))+W_ijk_abc(RNB,RNA,1:RNBLOCKC,NJ,NI,:,KMAX)

                        END DO
                        END DO
                     END DO
                     END DO

                     DEALLOCATE (W_ijk_abc)
                     ALLOCATE (W_ijk_abc(RNBLOCKC, RNBLOCKB, RNBLOCKA, VBMAX, VBMAX, VBMAX, WDES%NKPTS*WDES%NKPTS*MY_NKPTS))
                     W_ijk_abc(:, :, :, :, :, :, KMAX) = (0.0_qs, 0.0_qs)
! ijk_abc=ijk_abc+ kji_cba
                     CALL CALC_WINT(KPTS_MKPTS(KK), KJ, KI, KC, KB, KA, NC, NB, NA, RNBLOCKC, RNBLOCKB, RNBLOCKA, WDES, KPOINTS)
                     DO NI = 1, VBMAX
                     DO NJ = 1, VBMAX
                        DO RNA = 1, RNBLOCKA
                        DO RNB = 1, RNBLOCKB

                           Wint_ijk_abc(RNA,RNB,1:RNBLOCKC,NI,NJ,:,KINDEX(KI,KJ,KK))=Wint_ijk_abc(RNA,RNB,1:RNBLOCKC,NI,NJ,:,KINDEX(KI,KJ,KK))+W_ijk_abc(1:RNBLOCKC,RNB,RNA,:,NJ,NI,KMAX)
                        END DO
                        END DO
                     END DO
                     END DO

                     DEALLOCATE (W_ijk_abc)
                     ALLOCATE (W_ijk_abc(RNBLOCKA, RNBLOCKC, RNBLOCKB, VBMAX, VBMAX, VBMAX, WDES%NKPTS*WDES%NKPTS*MY_NKPTS))
                     W_ijk_abc(:, :, :, :, :, :, KMAX) = (0.0_qs, 0.0_qs)
! ijk_abc=ijk_abc+ ikj_acb
                     CALL CALC_WINT(KI, KPTS_MKPTS(KK), KJ, KA, KC, KB, NA, NC, NB, RNBLOCKA, RNBLOCKC, RNBLOCKB, WDES, KPOINTS)
                     DO NI = 1, VBMAX
                     DO NJ = 1, VBMAX
                        DO RNA = 1, RNBLOCKA
                        DO RNB = 1, RNBLOCKB

                           Wint_ijk_abc(RNA,RNB,1:RNBLOCKC,NI,NJ,:,KINDEX(KI,KJ,KK))=Wint_ijk_abc(RNA,RNB,1:RNBLOCKC,NI,NJ,:,KINDEX(KI,KJ,KK))+W_ijk_abc(RNA,1:RNBLOCKC,RNB,NI,:,NJ,KMAX)

                        END DO
                        END DO
                     END DO
                     END DO

                     DEALLOCATE (W_ijk_abc)
                     ALLOCATE (W_ijk_abc(RNBLOCKC, RNBLOCKA, RNBLOCKB, VBMAX, VBMAX, VBMAX, WDES%NKPTS*WDES%NKPTS*MY_NKPTS))
                     W_ijk_abc(:, :, :, :, :, :, KMAX) = (0.0_qs, 0.0_qs)
! ijk_abc=ijk_abc+ kij_cab
                     CALL CALC_WINT(KPTS_MKPTS(KK), KI, KJ, KC, KA, KB, NC, NA, NB, RNBLOCKC, RNBLOCKA, RNBLOCKB, WDES, KPOINTS)
                     DO NI = 1, VBMAX
                     DO NJ = 1, VBMAX
                        DO RNA = 1, RNBLOCKA
                        DO RNB = 1, RNBLOCKB

                           Wint_ijk_abc(RNA,RNB,1:RNBLOCKC,NI,NJ,:,KINDEX(KI,KJ,KK))=Wint_ijk_abc(RNA,RNB,1:RNBLOCKC,NI,NJ,:,KINDEX(KI,KJ,KK))+W_ijk_abc(1:RNBLOCKC,RNA,RNB,:,NI,NJ,KMAX)

                        END DO
                        END DO
                     END DO
                     END DO

                     DEALLOCATE (W_ijk_abc)
                     ALLOCATE (W_ijk_abc(RNBLOCKB, RNBLOCKC, RNBLOCKA, VBMAX, VBMAX, VBMAX, WDES%NKPTS*WDES%NKPTS*MY_NKPTS))
                     W_ijk_abc(:, :, :, :, :, :, KMAX) = (0.0_qs, 0.0_qs)
! ijk_abc=ijk_abc+ jki_bca
                     CALL CALC_WINT(KJ, KPTS_MKPTS(KK), KI, KB, KC, KA, NB, NC, NA, RNBLOCKB, RNBLOCKC, RNBLOCKA, WDES, KPOINTS)
                     DO NI = 1, VBMAX
                     DO NJ = 1, VBMAX
                        DO RNA = 1, RNBLOCKA
                        DO RNB = 1, RNBLOCKB

                           Wint_ijk_abc(RNA,RNB,1:RNBLOCKC,NI,NJ,:,KINDEX(KI,KJ,KK))=Wint_ijk_abc(RNA,RNB,1:RNBLOCKC,NI,NJ,:,KINDEX(KI,KJ,KK))+W_ijk_abc(RNB,1:RNBLOCKC,RNA,NJ,:,NI,KMAX)

                        END DO
                        END DO
                     END DO
                     END DO

                     DEALLOCATE (W_ijk_abc)
                     ALLOCATE (W_ijk_abc(RNBLOCKA, RNBLOCKB, RNBLOCKC, VBMAX, VBMAX, VBMAX, WDES%NKPTS*WDES%NKPTS*MY_NKPTS))
                     W_ijk_abc(:, :, :, :, :, :, KMAX) = (0.0_qs, 0.0_qs)
! ijk_abc=ijk_abc+ijk_abc
                     CALL CALC_WINT(KI, KJ, KPTS_MKPTS(KK), KA, KB, KC, NA, NB, NC, RNBLOCKA, RNBLOCKB, RNBLOCKC, WDES, KPOINTS)
                  Wint_ijk_abc(1:RNBLOCKA,1:RNBLOCKB,1:RNBLOCKC,:,:,:,KINDEX(KI,KJ,KK))=Wint_ijk_abc(1:RNBLOCKA,1:RNBLOCKB,1:RNBLOCKC,:,:,:,KINDEX(KI,KJ,KK))+W_ijk_abc(:,:,:,:,:,:,KMAX)
                  END DO
                  END DO
               END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              add V_ijk_abc to Wint_ijk_abc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               DO KI = 1, WDES%NKPTS

                  CALL BCAST2ALL_FTOD_IA(WDES, KI, 1, PW_IA_TMP, OC_IA_TMP)

                  DO KJ = 1, WDES%NKPTS

                     CALL BCAST2ALL_FTOD_IA(WDES, KJ, 2, PW_IA_TMP2, OC_IA_TMP2)

                     DO KK = 1, MY_NKPTS
                        KC = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KI)) + &
                                WDES%VKPT(:, KJ) + WDES%VKPT(:, KPTS_MKPTS(KK)) - WDES%VKPT(:, KA) - WDES%VKPT(:, KB), KPOINTS_FULL)

                        KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KI)) - &
                                                 WDES%VKPT(:, KA), KPOINTS_FULL)

!!!!!!!!!!!!!!!!!!!!! <ij|ab>t_k^c
                        IF (KC == KPTS_MKPTS(KK)) THEN
                           KTMP = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KI)) + &
                                                      WDES%VKPT(:, KJ) - WDES%VKPT(:, KA), KPOINTS_FULL)
                           IF (KTMP /= KB) WRITE (*, *) 'ahso1'
                           DO NI = 1, VBMAX
                           DO NJ = 1, VBMAX

                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                              PW_IA_TMP(1, 1, NI, KQ), OC_IA_TMP(1, 1, NI, KQ), &
                                              NUNOCC, &
                                              PW_IA_TMP2(1, 1, NJ, KQ), OC_IA_TMP2(1, 1, NJ, KQ), &
                                              NUNOCC, &
                                              VVOO(1, 1, NI, NJ), VVOO_R(1, 1, NI, NJ), zero)

                              IF (LORBREAL) THEN
                              ELSE
                                 VVOO = CONJG(VVOO)
                                 DO NK = 1, VBMAX
                                 DO RNC = 1, RNBLOCKC
                                 DO RNA = 1, RNBLOCKA
                                 DO RNB = 1, RNBLOCKB

                                    Wint_ijk_abc(RNA, RNB, RNC, NI, NJ, NK, KINDEX(KI, KJ, KK)) = Wint_ijk_abc(RNA, RNB, RNC, NI, NJ, NK, KINDEX(KI, KJ, KK)) + &
                                             ((VVOO((NA + RNA - 1), (NB + RNB - 1), NI, NJ))*(T1(NC + RNC - 1, NK, KPTS_MKPTS(KK))))
                                 END DO
                                 END DO
                                 END DO
                                 END DO
                              END IF

                           END DO
                           END DO

                        END IF

!!!!!!!!!!!!!!!!!!!!! <ik|ac>t_j^b
                        IF (KJ == KB) THEN

                           DO NI = 1, VBMAX
                           DO NK = 1, VBMAX

                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                              PW_IA_TMP(1, 1, NI, KQ), OC_IA_TMP(1, 1, NI, KQ), &
                                              NUNOCC, &
                                              FTOD_PW_IA(1, 1, NK, KQ, KK, 2), FTOD_OC_IA(1, 1, NK, KQ, KK, 2), &
                                              NUNOCC, &
                                              VVOO(1, 1, NI, NK), VVOO_R(1, 1, NI, NK), zero)

                              IF (LORBREAL) THEN
                              ELSE
                                 VVOO = CONJG(VVOO)
                                 DO NJ = 1, VBMAX
                                 DO RNC = 1, RNBLOCKC
                                 DO RNB = 1, RNBLOCKB
                                    Wint_ijk_abc(1:RNBLOCKA, RNB, RNC, NI, NJ, NK, KINDEX(KI, KJ, KK)) = Wint_ijk_abc(:, RNB, RNC, NI, NJ, NK, KINDEX(KI, KJ, KK)) + &
                                                       (VVOO(NA:(NA + RNBLOCKA - 1), NC + RNC - 1, NI, NK)*T1(NB + RNB - 1, NJ, KJ))
                                 END DO
                                 END DO
                                 END DO
                              END IF

                           END DO
                           END DO

                        END IF

                        KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KK)) - &
                                                 WDES%VKPT(:, KC), KPOINTS_FULL)

!!!!!!!!!!!!!!!!!!!!! <kj|cb>t_i^a
                        IF (KI == KA) THEN

                           DO NJ = 1, VBMAX
                           DO NK = 1, VBMAX

                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                              FTOD_PW_IA(1, 1, NK, KQ, KK, 1), FTOD_OC_IA(1, 1, NK, KQ, KK, 1), &
                                              NUNOCC, &
                                              PW_IA_TMP2(1, 1, NJ, KQ), OC_IA_TMP2(1, 1, NJ, KQ), &
                                              NUNOCC, &
                                              VVOO(1, 1, NK, NJ), VVOO_R(1, 1, NK, NJ), zero)

                              IF (LORBREAL) THEN
                              ELSE
                                 VVOO = CONJG(VVOO)
                                 DO NI = 1, VBMAX
                                 DO RNC = 1, RNBLOCKC
                                 DO RNB = 1, RNBLOCKB

                                    Wint_ijk_abc(:, RNB, RNC, NI, NJ, NK, KINDEX(KI, KJ, KK)) = Wint_ijk_abc(:, RNB, RNC, NI, NJ, NK, KINDEX(KI, KJ, KK)) + &
                                                       (VVOO(NC + RNC - 1, NB + RNB - 1, NK, NJ)*T1(NA:(NA + RNBLOCKA - 1), NI, KI))
                                 END DO
                                 END DO
                                 END DO
                              END IF

                           END DO
                           END DO

                        END IF

                     END DO
                  END DO
               END DO

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! EVALUATE             4./3.*W_ijk_abc+2./3.*W_kij_abc-2.*W_ikj_abc
!!!!!!!!!!!!!!!!!!!!!!!!!

               W_ijk_abc = (4.0_qs, 0.0_qs)/(3.0_qs, 0.0_qs)*Wint_ijk_abc
               DO KJ = 1, WDES%NKPTS
                  IF (PROCS_KPTS(KJ) == ME) THEN
                      CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (RNBLOCKA)*(RNBLOCKB)*(RNBLOCKC), (VBMAX*VBMAX*VBMAX*WDES%NKPTS*WDES%NKPTS), &
                                  Wint_ijk_abc(1, 1, 1, 1, 1, 1, KINDEX(1, 1, MKPTS_KPTS(KJ))), (RNBLOCKA)*(RNBLOCKB)*(RNBLOCKC))
                     Wtmp_ijk_abc(:,:,:,:,:,:,KINDEX(1,1,1):KINDEX(WDES%NKPTS,WDES%NKPTS,1))=Wint_ijk_abc(:,:,:,:,:,:,KINDEX(1,1,MKPTS_KPTS(KJ)):KINDEX(WDES%NKPTS,WDES%NKPTS,MKPTS_KPTS(KJ)))
                  ELSE
                     CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (RNBLOCKA)*(RNBLOCKB)*(RNBLOCKC), (VBMAX*VBMAX*VBMAX*WDES%NKPTS*WDES%NKPTS), &
                     Wtmp_ijk_abc(1, 1, 1, 1, 1, 1, KINDEX(1, 1, 1)), (RNBLOCKA)*(RNBLOCKB)*(RNBLOCKC), 0, PROCS_KPTS(KJ))
                  END IF

                  DO KI = 1, WDES%NKPTS
                  DO KK = 1, MY_NKPTS

                     DO NI = 1, VBMAX
                     DO NJ = 1, VBMAX
                     DO NK = 1, VBMAX
                        W_ijk_abc(:,:,:,NI,NJ,NK,KINDEX(KI,KJ,KK))=W_ijk_abc(:,:,:,NI,NJ,NK,KINDEX(KI,KJ,KK))+(2.0_qs,0.0_qs)/3.0_qs*Wtmp_ijk_abc(:,:,:,NK,NI,NJ,KINDEX(KPTS_MKPTS(KK),KI,1))
                        W_ijk_abc(:,:,:,NI,NJ,NK,KINDEX(KI,KJ,KK))=W_ijk_abc(:,:,:,NI,NJ,NK,KINDEX(KI,KJ,KK))-2.0_qs*Wtmp_ijk_abc(:,:,:,NI,NK,NJ,KINDEX(KI,KPTS_MKPTS(KK),1))
                     END DO
                     END DO
                     END DO

                  END DO
                  END DO
               END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              remove V_ijk_abc from Wint_ijk_abc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               DO KI = 1, WDES%NKPTS

                  CALL BCAST2ALL_FTOD_IA(WDES, KI, 1, PW_IA_TMP, OC_IA_TMP)

                  DO KJ = 1, WDES%NKPTS

                     CALL BCAST2ALL_FTOD_IA(WDES, KJ, 2, PW_IA_TMP2, OC_IA_TMP2)

                     DO KK = 1, MY_NKPTS
                        KC = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KI)) + &
                                WDES%VKPT(:, KJ) + WDES%VKPT(:, KPTS_MKPTS(KK)) - WDES%VKPT(:, KA) - WDES%VKPT(:, KB), KPOINTS_FULL)

                        KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KI)) - &
                                                 WDES%VKPT(:, KA), KPOINTS_FULL)

!!!!!!!!!!!!!!!!!!!!! <ij|ab>t_k^c
                        IF (KC == KPTS_MKPTS(KK)) THEN

                           DO NI = 1, VBMAX
                           DO NJ = 1, VBMAX

                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                              PW_IA_TMP(1, 1, NI, KQ), OC_IA_TMP(1, 1, NI, KQ), &
                                              NUNOCC, &
                                              PW_IA_TMP2(1, 1, NJ, KQ), OC_IA_TMP2(1, 1, NJ, KQ), &
                                              NUNOCC, &
                                              VVOO(1, 1, NI, NJ), VVOO_R(1, 1, NI, NJ), zero)

                              IF (LORBREAL) THEN
                              ELSE
                                 VVOO = CONJG(VVOO)
                                 DO NK = 1, VBMAX
                                 DO RNC = 1, RNBLOCKC
                                 DO RNA = 1, RNBLOCKA
                                 DO RNB = 1, RNBLOCKB

                                    Wint_ijk_abc(RNA, RNB, RNC, NI, NJ, NK, KINDEX(KI, KJ, KK)) = Wint_ijk_abc(RNA, RNB, RNC, NI, NJ, NK, KINDEX(KI, KJ, KK)) - &
                                             ((VVOO((NA + RNA - 1), (NB + RNB - 1), NI, NJ))*(T1(NC + RNC - 1, NK, KPTS_MKPTS(KK))))
                                 END DO
                                 END DO
                                 END DO
                                 END DO
                              END IF

                           END DO
                           END DO

                        END IF

!!!!!!!!!!!!!!!!!!!!! <ik|ac>t_j^b
                        IF (KJ == KB) THEN

                           DO NI = 1, VBMAX
                           DO NK = 1, VBMAX

                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                              PW_IA_TMP(1, 1, NI, KQ), OC_IA_TMP(1, 1, NI, KQ), &
                                              NUNOCC, &
                                              FTOD_PW_IA(1, 1, NK, KQ, KK, 2), FTOD_OC_IA(1, 1, NK, KQ, KK, 2), &
                                              NUNOCC, &
                                              VVOO(1, 1, NI, NK), VVOO_R(1, 1, NI, NK), zero)

                              IF (LORBREAL) THEN
                              ELSE
                                 VVOO = CONJG(VVOO)
                                 DO NJ = 1, VBMAX
                                 DO RNC = 1, RNBLOCKC
                                 DO RNB = 1, RNBLOCKB

                                    Wint_ijk_abc(1:RNBLOCKA, RNB, RNC, NI, NJ, NK, KINDEX(KI, KJ, KK)) = Wint_ijk_abc(:, RNB, RNC, NI, NJ, NK, KINDEX(KI, KJ, KK)) - &
                                                       (VVOO(NA:(NA + RNBLOCKA - 1), NC + RNC - 1, NI, NK)*T1(NB + RNB - 1, NJ, KJ))
                                 END DO
                                 END DO
                                 END DO
                              END IF

                           END DO
                           END DO

                        END IF

                        KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KK)) - &
                                                 WDES%VKPT(:, KC), KPOINTS_FULL)

!!!!!!!!!!!!!!!!!!!!! <kj|cb>t_i^a
                        IF (KI == KA) THEN

                           DO NJ = 1, VBMAX
                           DO NK = 1, VBMAX

                              CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                                              FTOD_PW_IA(1, 1, NK, KQ, KK, 1), FTOD_OC_IA(1, 1, NK, KQ, KK, 1), &
                                              NUNOCC, &
                                              PW_IA_TMP2(1, 1, NJ, KQ), OC_IA_TMP2(1, 1, NJ, KQ), &
                                              NUNOCC, &
                                              VVOO(1, 1, NK, NJ), VVOO_R(1, 1, NK, NJ), zero)

                              IF (LORBREAL) THEN
                              ELSE
                                 VVOO = CONJG(VVOO)
                                 DO NI = 1, VBMAX
                                 DO RNC = 1, RNBLOCKC
                                 DO RNB = 1, RNBLOCKB

                                    Wint_ijk_abc(:, RNB, RNC, NI, NJ, NK, KINDEX(KI, KJ, KK)) = Wint_ijk_abc(:, RNB, RNC, NI, NJ, NK, KINDEX(KI, KJ, KK)) - &
                                                       (VVOO(NC + RNC - 1, NB + RNB - 1, NK, NJ)*T1(NA:(NA + RNBLOCKA - 1), NI, KI))
                                 END DO
                                 END DO
                                 END DO
                              END IF

                           END DO
                           END DO

                        END IF

                     END DO
                  END DO
               END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!     CALCULATE E_T
!!!!!!!!!!!!!!!!!!!!!!!!!!!

               DO KI = 1, WDES%NKPTS
               DO KJ = 1, WDES%NKPTS
               DO KK = 1, MY_NKPTS

                  KC = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KI)) + &
                                WDES%VKPT(:, KJ) + WDES%VKPT(:, KPTS_MKPTS(KK)) - WDES%VKPT(:, KA) - WDES%VKPT(:, KB), KPOINTS_FULL)

! ijk_abc=ijk_abc+ jik_bac
                  DO NK = 1, VBMAX
                  DO NJ = 1, VBMAX
                  DO NI = 1, VBMAX
                     DO RNC = 1, RNBLOCKC
                     DO RNB = 1, RNBLOCKB
                     DO RNA = 1, RNBLOCKA

                        ETRIPLES = ETRIPLES + &
   (Wint_ijk_abc(RNA, RNB, RNC, NI, NJ, NK, KINDEX(KI, KJ, KK)))*GCONJG(W_ijk_abc(RNA, RNB, RNC, NI, NJ, NK, KINDEX(KI, KJ, KK)))/ &
                           REAL(W%CELTOT(NI,KI,1)+W%CELTOT(NJ,KJ,1)+W%CELTOT(NK,KPTS_MKPTS(KK),1)-W%CELTOT(RNA+NA-1+VBMAX,KA,1)-W%CELTOT(RNB+NB-1+VBMAX,KB,1)-W%CELTOT(RNC+NC-1+VBMAX,KC,1),kind=q)

                        KETRIPLES(KA) = KETRIPLES(KA) + &
   (Wint_ijk_abc(RNA, RNB, RNC, NI, NJ, NK, KINDEX(KI, KJ, KK)))*GCONJG(W_ijk_abc(RNA, RNB, RNC, NI, NJ, NK, KINDEX(KI, KJ, KK)))/ &
                           REAL(W%CELTOT(NI,KI,1)+W%CELTOT(NJ,KJ,1)+W%CELTOT(NK,KPTS_MKPTS(KK),1)-W%CELTOT(RNA+NA-1+VBMAX,KA,1)-W%CELTOT(RNB+NB-1+VBMAX,KB,1)-W%CELTOT(RNC+NC-1+VBMAX,KC,1),kind=q)

                     END DO
                     END DO
                     END DO
                  END DO
                  END DO
                  END DO

               END DO
               END DO
               END DO

            END DO
            END DO
            END DO

         END DO
      END DO

      call date_and_time(values=time_array2)
      EMSTRIPLES = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                    time_array2(8) - time_array1(8))

      CALLMPI(M_sum_d(WGW%COMM_INTER, ETRIPLES, 2))
         WRITE(*,*)'(T) contribution to the correlation energy is',ETRIPLES*REAL(1.0_q/WDES%NKPTS/WDES%NKPTS/WDES%NKPTS/WDES%NKPTS/WDES%NKPTS,kind=q) !,T2
         WRITE(IO%IU6,*)'(T) contribution to the correlation energy is',ETRIPLES*REAL(1.0_q/WDES%NKPTS/WDES%NKPTS/WDES%NKPTS/WDES%NKPTS/WDES%NKPTS,kind=q) !,T2
      WRITE (*, *) 'calculating triples costs', emstriples
      WRITE (IO%IU6, *) 'calculating triples costs', emstriples
WRITE(*,*)'k-point dependent triples energy is ',KETRIPLES*REAL(1.0_q/WDES%NKPTS/WDES%NKPTS/WDES%NKPTS/WDES%NKPTS/WDES%NKPTS,kind=q)
         WRITE(IO%IU6,*)'k-point dependent triples energy is ',KETRIPLES*REAL(1.0_q/WDES%NKPTS/WDES%NKPTS/WDES%NKPTS/WDES%NKPTS/WDES%NKPTS,kind=q)

   END Subroutine CALC_TRIPLES

!***********************************************************************
!
! This routine computes the W_ijk^abc intermediate needed for the (T) contribution to the CCSD(T) energy.
!
!***********************************************************************

   SUBROUTINE CALC_WINT(KI, KJ, KK, KA, KB, KC, NA, NB, NC, RNBLOCKA, RNBLOCKB, RNBLOCKC, WDES, KPOINTS)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(kpoints_struct) KPOINTS
      INTEGER ITERATION, NBLOCK, MKIMAX, MKIMIN
      INTEGER KI, KJ, KK, KA, KB, KC, NI, NJ, NK, NM, NA, NB, NC, KQ, KQ_, KF, KM, KMAX, KE, KE2, KM2
      INTEGER RNBLOCKA, RNBLOCKB, RNBLOCKC, RNA, RNB, RNC, ierr
      COMPLEX(q) :: VROR(NUNOCC, RNBLOCKB, VBMAX, RNBLOCKC)
      COMPLEX(q) :: MPW_AB(NGVECTOR, NUNOCC, NUNOCC), MOC_AB(NHVECTOR, NUNOCC, NUNOCC)
      REAL(q) :: VROR_R(NUNOCC, RNBLOCKB, VBMAX, RNBLOCKC)
      COMPLEX(qs) :: VRRO_S(NUNOCC, RNBLOCKB, RNBLOCKC, VBMAX)
      COMPLEX(qs) :: VROR_S(NUNOCC, RNBLOCKB, VBMAX, RNBLOCKC)
      REAL(qs) :: VRRO_S_R(NUNOCC, RNBLOCKA, RNBLOCKB, VBMAX)

      COMPLEX(q) :: OROO(VBMAX, RNBLOCKC, VBMAX, VBMAX)
      REAL(q) :: OROO_R(VBMAX, RNBLOCKC, VBMAX, VBMAX)
      COMPLEX(qs) :: OROO_S(VBMAX, RNBLOCKC, VBMAX, VBMAX)
      REAL(qs) :: OROO_S_R(VBMAX, RNBLOCKC, VBMAX, VBMAX)
      COMPLEX(q) KW
      INTEGER :: time_array1(8), time_array2(8)

      call date_and_time(values=time_array1)
      KW = (1.0_qs, 0.0_qs)*KPOINTS_FULL%WTKPT(1)

      KMAX = WDES%NKPTS*WDES%NKPTS*MY_NKPTS

      KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KC)) - &
                               WDES%VKPT(:, KK), KPOINTS_FULL)
      KE = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KB)) + &
                               WDES%VKPT(:, KQ), KPOINTS_FULL)
      KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KI)) - &
                               WDES%VKPT(:, KA), KPOINTS_FULL)

      KM = KPOINT_IN_FULL_GRID(-WDES%VKPT(:, (KQ)) + &
                               WDES%VKPT(:, KB), KPOINTS_FULL)
      KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KC)) - &
                               WDES%VKPT(:, KK), KPOINTS_FULL)

      CALL MPI_ALLREDUCE(KK, MKIMAX, 1, MPI_INTEGER, MPI_MAX, WDES%COMM, ierr)
      CALL MPI_ALLREDUCE(KK, MKIMIN, 1, MPI_INTEGER, MPI_MIN, WDES%COMM, ierr)
      IF (MKIMAX /= MKIMIN) THEN
         PW_IA_TMP(:, :, :, KQ) = FTOD_PW_IA(:, :, :, KQ, MKPTS_KPTS(KK), 2)
         IF (ASSOCIATED(H)) THEN
            OC_IA_TMP(:, :, :, KQ) = FTOD_OC_IA(:, :, :, KQ, MKPTS_KPTS(KK), 2)
         END IF
      ELSE
         CALL BCAST2ALL_FTOD_IA(WDES, KK, 2, PW_IA_TMP, OC_IA_TMP) !problem for KI=KK ! solve KI/KK locally
      END IF

      KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KC)) - &
                               WDES%VKPT(:, KK), KPOINTS_FULL)
      DO nm = 1, VBMAX
         PW_AI_TMP(:, nm, :, KQ) = (PW_IA_TMP(:, :, nm, KQ))
      END DO
      IF (ASSOCIATED(H)) THEN
      DO nm = 1, VBMAX
         OC_AI_TMP(:, nm, :, KQ) = (OC_IA_TMP(:, :, nm, KQ))
      END DO
      END IF

      CALL MPI_ALLREDUCE(KE, MKIMAX, 1, MPI_INTEGER, MPI_MAX, WDES%COMM, ierr)
      CALL MPI_ALLREDUCE(KE, MKIMIN, 1, MPI_INTEGER, MPI_MIN, WDES%COMM, ierr)
      IF (MKIMAX /= MKIMIN) THEN
         CALL TRANSFER_FTOD_AB(WDES, KE, KQ, 1, PW_AB_TMP, OC_AB_TMP)
      ELSE
         CALL BCAST2ALL_FTOD_AB(WDES, KE, 1, PW_AB_TMP, OC_AB_TMP) !problem for KB=KC
      END IF

      KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KC)) - &
                               WDES%VKPT(:, KK), KPOINTS_FULL)
      MPW_AB(:, :, :) = PW_AB_TMP(:, :, :, KQ)
      IF (ASSOCIATED(H)) THEN
         MOC_AB(:, :, :) = OC_AB_TMP(:, :, :, KQ)
      END IF
      DO nm = 1, NUNOCC
         PW_AB_TMP(:, :, nm, KQ) = MPW_AB(:, nm, :)
      END DO
      IF (ASSOCIATED(H)) THEN
      DO nm = 1, NUNOCC
         OC_AB_TMP(:, :, nm, KQ) = MOC_AB(:, nm, :)
      END DO
      END IF

      CALL TRANSFER_T2_KF_KJ_KK(WDES, KE, KJ, KI, VVOO_S, VVOO_S_R)

      CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                      PW_AB_TMP(1, 1, NB, KQ), OC_AB_TMP(1, 1, NB, KQ), &
                      (NUNOCC)*(RNBLOCKB), &
                      PW_AI_TMP(1, 1, NC, KQ), OC_AI_TMP(1, 1, NC, KQ), &
                      VBMAX*RNBLOCKC, &
                      VROR(1, 1, 1, 1), VROR_R(1, 1, 1, 1), zero)
! ebkc
      IF (LORBREAL) THEN
      ELSE
         DO RNC = 1, RNBLOCKC
         DO RNB = 1, RNBLOCKB
            VRRO_S(:, RNB, RNC, :) = CONJG(VROR(:, RNB, :, RNC))
         END DO
         END DO
         DO NI = 1, VBMAX
         DO NJ = 1, VBMAX
         DO NK = 1, VBMAX
            CALL CGEMM('t', 'n', (RNBLOCKA), (RNBLOCKB*RNBLOCKC), (NUNOCC), &
                       (1.0_qs, 0.0_qs), VVOO_S(1, NA, NJ, NI), (NUNOCC), &
                       VRRO_S(1, 1, 1, NK), (NUNOCC), &
                       (1.0_qs, 0.0_qs), W_ijk_abc(1, 1, 1, NI, NJ, NK, KMAX), (RNBLOCKA))
         END DO
         END DO
         END DO
      END IF

      KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                               WDES%VKPT(:, KA), KPOINTS_FULL)

!        -\sum_m
      CALL MPI_ALLREDUCE(KJ, MKIMAX, 1, MPI_INTEGER, MPI_MAX, WDES%COMM, ierr)
      CALL MPI_ALLREDUCE(KJ, MKIMIN, 1, MPI_INTEGER, MPI_MIN, WDES%COMM, ierr)
      IF (MKIMAX /= MKIMIN) THEN
         PW_IJ_TMP(:, :, :, :) = FTOD_PW_IJ(:, :, :, :, MKPTS_KPTS(KJ), 1)
         IF (ASSOCIATED(H)) THEN
            OC_IJ_TMP(:, :, :, :) = FTOD_OC_IJ(:, :, :, :, MKPTS_KPTS(KJ), 1)
         END IF
      ELSE
         CALL BCAST2ALL_FTOD_IJ(WDES, KJ, 1, PW_IJ_TMP, OC_IJ_TMP)  !problem KJ=KK !solve -> KI/KK stored locally
      END IF

      CALL MPI_ALLREDUCE(KK, MKIMAX, 1, MPI_INTEGER, MPI_MAX, WDES%COMM, ierr)
      CALL MPI_ALLREDUCE(KK, MKIMIN, 1, MPI_INTEGER, MPI_MIN, WDES%COMM, ierr)
      IF (MKIMAX /= MKIMIN) THEN
         PW_IA_TMP(:, :, :, :) = FTOD_PW_IA(:, :, :, :, MKPTS_KPTS(KK), 2)
         IF (ASSOCIATED(H)) THEN
            OC_IA_TMP(:, :, :, :) = FTOD_OC_IA(:, :, :, :, MKPTS_KPTS(KK), 2)
         END IF
      ELSE
         CALL BCAST2ALL_FTOD_IA(WDES, KK, 2, PW_IA_TMP, OC_IA_TMP) !problem for KI=KK ! solve KI/KK locally
      END IF

      KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KJ)) - &
                               WDES%VKPT(:, (KM)), KPOINTS_FULL)

      CALL TRANSFER_T2_KF_KJ_KK(WDES, KA, KI, KM, VVOO_S, VVOO_S_R)

      DO NJ = 1, VBMAX
      DO NK = 1, VBMAX
         CALL CONTR_FTOD(NGVECTOR, NHVECTOR, &
                         PW_IJ_TMP(1, 1, NJ, KQ), OC_IJ_TMP(1, 1, NJ, KQ), &
                         VBMAX, &
                         PW_IA_TMP(1, NC, NK, KQ), OC_IA_TMP(1, NC, NK, KQ), &
                         RNBLOCKC, &
                         OROO(1, 1, NJ, NK), OROO_R(1, 1, NJ, NK), zero)
      END DO
      END DO

      IF (LORBREAL) THEN
      ELSE
         OROO_S = CONJG(OROO)
         DO NI = 1, VBMAX
         DO NJ = 1, VBMAX
         DO NK = 1, VBMAX
         DO NM = 1, VBMAX
            DO RNC = 1, RNBLOCKC
                     W_ijk_abc(:,:,RNC,NI,NJ,NK,KMAX)=W_ijk_abc(:,:,RNC,NI,NJ,NK,KMAX)-OROO_S(NM,RNC,NJ,NK)*VVOO_S((NA):(NA+RNBLOCKA-1),(NB):(NB+RNBLOCKB-1),NI,NM)
            END DO
         END DO
         END DO
         END DO
         END DO
      END IF

      call date_and_time(values=time_array2)
      EMSCALCWINT = EMSCALCWINT + ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
                                   time_array2(8) - time_array1(8))

   END SUBROUTINE CALC_WINT

!***********************************************************************
!
! Routine to broadcast the W_ijk^abc intermediate needed for the (T) contribution to the CCSD(T) energy.
!
!***********************************************************************

   Subroutine BCAST2ALL_WINT(WDES, NBLOCK, KK, W_MTMP, W_MTMP_R)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KK, CHIL, CHIM, NBLOCK
      COMPLEX(qs) :: W_MTMP(:, :, :, :, :, :, :)
      REAL(qs) :: W_MTMP_R(:, :, :, :, :, :, :)
      integer :: time_array1(8), time_array2(8), ems

      call BLACS_PINFO(ME, PROCS)
      call date_and_time(values=time_array1)

      IF (LORBREAL) THEN
      ELSE
         IF (PROCS_KPTS(KK) == ME) THEN
            CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NBLOCK)*(NBLOCK)*(NBLOCK), (VBMAX*VBMAX*VBMAX*WDES%NKPTS*WDES%NKPTS), &
                         Wint_ijk_abc(1, 1, 1, 1, 1, 1, KINDEX(1, 1, MKPTS_KPTS(KK))), (NBLOCK)*(NBLOCK)*(NBLOCK))
            W_MTMP(:, :, :, :, :, :, 1) = Wint_ijk_abc(:, :, :, :, :, :, KINDEX(1, 1, MKPTS_KPTS(KK)))
         ELSE
            CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NBLOCK)*(NBLOCK)*(NBLOCK), (VBMAX*VBMAX*VBMAX*WDES%NKPTS*WDES%NKPTS), &
                        Wint_ijk_abc(1, 1, 1, 1, 1, 1, KINDEX(1, 1, MKPTS_KPTS(KK))), (NBLOCK)*(NBLOCK)*(NBLOCK), 0, PROCS_KPTS(KK))
         END IF
      END IF

      call date_and_time(values=time_array2)
      ems = ((time_array2(6)*60 + time_array2(7))*1000 - (time_array1(6)*60 + time_array1(7))*1000 + &
             time_array2(8) - time_array1(8))
      TET2S = TET2S + ems

   END Subroutine BCAST2ALL_WINT

!***********************************************************************
!
! Communication routine for a particular k-vector tuple of the T2 amplitude.
!
!***********************************************************************

   SUBROUTINE TRANSFER_T2_KF_KJ_KK(WDES, KF, KJ, KK, XY, XY_R)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KI, MNCC, ATPRO, KF, KJ, KK
      GDEFS :: XY(:, :, :, :)
      REAL(qs) :: XY_R(:, :, :, :)
      INTEGER :: ORDERLIST(3, PROCS), tag, ierror, myrank
      INTEGER :: req(PROCS), recreq
      INTEGER :: status(MPI_STATUS_SIZE, PROCS)

      tag = 1
      call BLACS_PINFO(ME, PROCS)
      CALL MPI_COMM_RANK(WDES%comm, myrank, ierror)
      me = myrank
      ORDERLIST = -1
      IF (ME == 0) THEN
         ORDERLIST(1, 1) = KF
         ORDERLIST(2, 1) = KJ
         ORDERLIST(3, 1) = KK
         DO ATPRO = 1, PROCS - 1
            CALL MPI_RECV(ORDERLIST(1, ATPRO + 1), 3, MPI_INTEGER, ATPRO, tag, &
                          WDES%COMM%MPI_COMM, status(1, ATPRO), ierror)
         END DO
      ELSE
         call MPI_ISEND((/KF, KJ, KK/), 3, MPI_INTEGER, 0, tag, &
                        WDES%COMM%MPI_COMM, req(ME + 1), ierror)
         CALL MPI_WAIT(req(ME + 1), status(1, ME + 1), ierror)
      END IF

      IF (ME == 0) THEN
         CALL IGEBS2D(CONTXT_COLS, 'All', 'i-ring', PROCS*3, 1, ORDERLIST(1, 1), PROCS*3)
      ELSE
         CALL IGEBR2D(CONTXT_COLS, 'All', 'i-ring', PROCS*3, 1, ORDERLIST(1, 1), PROCS*3, 0, 0)
      END IF

      IF (LORBREAL) THEN
         DO ATPRO = 1, PROCS
            IF ((PROCS_KPTS(ORDERLIST(1, ATPRO)) == ME) .and. ((ATPRO - 1) /= ME)) THEN
               call MPI_ISEND(T2_R(1, 1, 1, 1, ORDERLIST(2, ATPRO), ORDERLIST(3, ATPRO), MKPTS_KPTS(ORDERLIST(1, ATPRO))), &
                              (PROCS*WDES%NBANDS - VBMAX)*(PROCS*WDES%NBANDS - VBMAX)*VBMAX*VBMAX, MPI_REAL8, ATPRO - 1, tag, &
                              WDES%COMM%MPI_COMM, req(ME + 1), ierror)
               CALL MPI_WAIT(req(ME + 1), status(1, ME + 1), ierror)
            END IF
            IF ((PROCS_KPTS(KF) /= ME) .and. ((ATPRO - 1) == ME)) THEN
               call MPI_RECV(XY_R(1, 1, 1, 1), &
                             (PROCS*WDES%NBANDS - VBMAX)*(PROCS*WDES%NBANDS - VBMAX)*VBMAX*VBMAX, MPI_REAL8, PROCS_KPTS(KF), tag, &
                             WDES%COMM%MPI_COMM, status(1, ME + 1), ierror)
            END IF
         END DO
         IF ((PROCS_KPTS(KF) == ME)) THEN
            XY_R(:, :, :, :) = T2_R(:, :, :, :, KJ, KK, MKPTS_KPTS(KF)) 
         END IF
      ELSE
         DO ATPRO = 1, PROCS
            IF ((PROCS_KPTS(ORDERLIST(1, ATPRO)) == ME) .and. ((ATPRO - 1) /= ME)) THEN
               call MPI_ISEND(T2(1, 1, 1, 1, ORDERLIST(2, ATPRO), ORDERLIST(3, ATPRO), MKPTS_KPTS(ORDERLIST(1, ATPRO))), &
                              (PROCS*WDES%NBANDS - VBMAX)*(PROCS*WDES%NBANDS - VBMAX)*VBMAX*VBMAX, MPI_COMPLEX8, ATPRO - 1, tag, &
                              WDES%COMM%MPI_COMM, req(ME + 1), ierror)
               CALL MPI_WAIT(req(ME + 1), status(1, ME + 1), ierror)
            END IF
            IF ((PROCS_KPTS(KF) /= ME) .and. ((ATPRO - 1) == ME)) THEN
               call MPI_RECV(XY(1, 1, 1, 1), &
                           (PROCS*WDES%NBANDS - VBMAX)*(PROCS*WDES%NBANDS - VBMAX)*VBMAX*VBMAX, MPI_COMPLEX8, PROCS_KPTS(KF), tag, &
                             WDES%COMM%MPI_COMM, status(1, ME + 1), ierror) !, recreq, ierror)
            END IF
         END DO
         IF ((PROCS_KPTS(KF) == ME)) THEN
            XY(:, :, :, :) = T2(:, :, :, :, KJ, KK, MKPTS_KPTS(KF)) 
         END IF

      END IF

   END SUBROUTINE TRANSFER_T2_KF_KJ_KK

!***********************************************************************
!
! Communication routine for a particular k-vector pair of the Coulomb-Vertex.
!
!***********************************************************************

   SUBROUTINE TRANSFER_FTOD_AB(WDES, KB, KQ, mmncc, MPW_AB_TMP, MOC_AB_TMP)
      use mkpoints
      use base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: KB, KQ, MNCC, ATPRO
      GDEF :: MPW_AB_TMP(:, :, :, :)
      GDEF :: MOC_AB_TMP(:, :, :, :)
      INTEGER :: ORDERLIST(2, PROCS), tag, ierror, myrank
      INTEGER :: req(PROCS), recreq, mmncc, ierr
      INTEGER :: status(MPI_STATUS_SIZE, PROCS)

      tag = 1
      call BLACS_PINFO(ME, PROCS)
      CALL MPI_COMM_RANK(WDES%comm, myrank, ierror)
      me = myrank
      ORDERLIST = -1
      IF (ME == 0) THEN
         ORDERLIST(1, 1) = KB
         ORDERLIST(2, 1) = KQ
         DO ATPRO = 1, PROCS - 1
            CALL MPI_RECV(ORDERLIST(1, ATPRO + 1), 2, MPI_INTEGER, ATPRO, tag, &
                          WDES%COMM%MPI_COMM, status(1, ATPRO), ierror)
         END DO
      ELSE
         call MPI_ISEND((/KB, KQ/), 2, MPI_INTEGER, 0, tag, &
                        WDES%COMM%MPI_COMM, req(ME + 1), ierror)
         CALL MPI_WAIT(req(ME + 1), status(1, ME + 1), ierror)
      END IF

      IF (ME == 0) THEN
         CALL IGEBS2D(CONTXT_COLS, 'All', 'i-ring', PROCS*2, 1, ORDERLIST(1, 1), PROCS*2)
      ELSE
         CALL IGEBR2D(CONTXT_COLS, 'All', 'i-ring', PROCS*2, 1, ORDERLIST(1, 1), PROCS*2, 0, 0)
      END IF

      DO ATPRO = 1, PROCS
         IF ((PROCS_KPTS(ORDERLIST(1, ATPRO)) == ME) .and. ((ATPRO - 1) /= ME)) THEN
            call MPI_ISEND(FTOD_PW_AB(1, 1, 1, ORDERLIST(2, ATPRO), MKPTS_KPTS(ORDERLIST(1, ATPRO)), mmncc), &
                           (PROCS*WDES%NBANDS - VBMAX)*(PROCS*WDES%NBANDS - VBMAX)*NGVECTOR, MPI_COMPLEX16, ATPRO - 1, tag, &
                           WDES%COMM%MPI_COMM, req(ME + 1), ierror)
            CALL MPI_WAIT(req(ME + 1), status(1, ME + 1), ierror)
         END IF
         IF ((PROCS_KPTS(KB) /= ME) .and. ((ATPRO - 1) == ME)) THEN
            call MPI_RECV(MPW_AB_TMP(1, 1, 1, KQ), &
                          (PROCS*WDES%NBANDS - VBMAX)*(PROCS*WDES%NBANDS - VBMAX)*NGVECTOR, MPI_COMPLEX16, PROCS_KPTS(KB), tag, &
                          WDES%COMM%MPI_COMM, status(1, ME + 1), ierror) 
         END IF
      END DO
      IF ((PROCS_KPTS(KB) == ME)) THEN
         MPW_AB_TMP(:, :, :, KQ) = FTOD_PW_AB(:, :, :, KQ, MKPTS_KPTS(KB), mmncc)
      END IF

      IF (ASSOCIATED(H)) THEN
         DO ATPRO = 1, PROCS
            IF ((PROCS_KPTS(ORDERLIST(1, ATPRO)) == ME) .and. ((ATPRO - 1) /= ME)) THEN
               call MPI_ISEND(FTOD_OC_AB(1, 1, 1, ORDERLIST(2, ATPRO), MKPTS_KPTS(ORDERLIST(1, ATPRO)), mmncc), &
                              (PROCS*WDES%NBANDS - VBMAX)*(PROCS*WDES%NBANDS - VBMAX)*NHVECTOR, MPI_COMPLEX16, ATPRO - 1, tag, &
                              WDES%COMM%MPI_COMM, req(ME + 1), ierror)
               CALL MPI_WAIT(req(ME + 1), status(1, ME + 1), ierror)
            END IF
            IF ((PROCS_KPTS(KB) /= ME) .and. ((ATPRO - 1) == ME)) THEN
               call MPI_RECV(MOC_AB_TMP(1, 1, 1, KQ), &
                             (PROCS*WDES%NBANDS - VBMAX)*(PROCS*WDES%NBANDS - VBMAX)*NHVECTOR, MPI_COMPLEX16, PROCS_KPTS(KB), tag, &
                             WDES%COMM%MPI_COMM, status(1, ME + 1), ierror)
            END IF
         END DO
         IF ((PROCS_KPTS(KB) == ME)) THEN
            MOC_AB_TMP(:, :, :, KQ) = FTOD_OC_AB(:, :, :, KQ, MKPTS_KPTS(KB), mmncc)
         END IF
      END IF

   END SUBROUTINE TRANSFER_FTOD_AB

!***********************************************************************
!
! Fock Martix setup routine
!
!***********************************************************************

   SUBROUTINE SETUP_FOCKMATRIX(W, WGW, WDES, IO)
      USE prec
      USE poscar
      USE pseudo
      USE wave_high
      USE full_kpoints
      USE mkpoints
      USE lattice
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      TYPE(in_struct) IO
      TYPE(wavespin) W
      INTEGER :: NI, NJ, NA, NB, KI, KJ, KQ, KQ_

      F_KI(:, :, :) = (0.0_qs, 0.0_qs)
      F_BA(:, :, :) = (0.0_qs, 0.0_qs)
      F_AI(:, :, :) = (0.0_qs, 0.0_qs)
      F_KC(:, :, :) = (0.0_qs, 0.0_qs)

      CALL FOCKM_IN(W,WGW,IO,WDES)
   END SUBROUTINE SETUP_FOCKMATRIX

!***********************************************************************
!
! Fock Martix reader routine
!
!***********************************************************************

   SUBROUTINE FOCKM_IN(W, WGW, IO, WDES)
      USE base
      IMPLICIT NONE
      TYPE(in_struct) IO
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      TYPE(wavespin) W
      INTEGER :: NB_TOT, NBI, IREC, NISP, NKPTS, KI, ISP, NBA
      COMPLEX(q), ALLOCATABLE :: FOCKM_LINE(:)
      LOGICAL :: CONSISTENT, ex

      CONSISTENT = .TRUE.
      F_KI = (0.0_qs, 0.0_qs)
      F_BA = (0.0_qs, 0.0_qs)

      call BLACS_GRIDINFO(CONTXT_GRID, NPROW, NPCOL, MYROW, MYCOL)

      IF ((MYROW == 0) .and. (MYCOL == 0)) THEN
         inquire (file="FOCKM", exist=ex)
         IF (.not. ex) THEN
            DO NBI = 1, VBMAX
               F_KI(NBI, NBI, :) = W%CELTOT(NBI, :, 1)
            END DO
            DO NBI = VBMAX + 1, WDES%NB_TOT
               F_BA(NBI - VBMAX, NBI - VBMAX, :) = W%CELTOT(NBI, :, 1)
            END DO
            RETURN
         END IF

         OPEN (UNIT=19, FILE='FOCKM', ACCESS='DIRECT', &
               FORM='UNFORMATTED', STATUS='UNKNOWN', RECL=((3)*IO%ICMPLX*2))
         READ (19, REC=1) NISP, NKPTS, NB_TOT
         IF (WDES%ISPIN /= NISP) THEN
            WRITE (IO%IU0, *) 'ERROR: number of SPINORS in FOCKM file differs from&
            &settings in the INCAR file'
            RETURN
         END IF
         IF (WDES%NKPTS /= NKPTS) THEN
            WRITE (IO%IU0, *) 'ERROR: number of kpoints in FOCKM file differs from&
            &settings in the KPOINTS file'
            RETURN
         END IF
         IF (WDES%NB_TOT /= NB_TOT) THEN
            WRITE (IO%IU0, *) 'WARNING: number of bands in FOCKM :', NB_TOT, '.  &
            & Now using ', WDES%NB_TOT, ' bands.'
         END IF

         CLOSE (19)
         ALLOCATE (FOCKM_LINE(NB_TOT))
         OPEN (UNIT=19, FILE='FOCKM', ACCESS='DIRECT', &
               FORM='UNFORMATTED', STATUS='UNKNOWN', RECL=((NB_TOT + 3)*IO%ICMPLX*2))
         IREC = 2
         DO ISP = 1, WDES%ISPIN
            DO KI = 1, WDES%NKPTS
               DO NBI = 1, NB_TOT
                  READ (19, REC=IREC) FOCKM_LINE

                  IREC = IREC + 1

                  IF (NBI > WDES%NB_TOT) CYCLE

                  IF (REAL(W%CELTOT(NBI, KI, ISP), KIND=q) /= REAL(FOCKM_LINE(NBI), KIND=q)) CONSISTENT = .FALSE.

                  W%CELTOT(NBI, KI, ISP) = CONJG(FOCKM_LINE(NBI))

                  DO NBA = VBMAX + 1, WDES%NB_TOT
                     IF ((NBI > VBMAX) .and. (NBI /= NBA)) THEN
                        F_BA(NBA - VBMAX, NBI - VBMAX, KI) = CONJG(FOCKM_LINE(NBA))
                     END IF
                  END DO
               END DO
            END DO
         END DO
         CLOSE (19)
      ELSE
         W%CELTOT = zero
      END IF

      CALLMPI(M_sum_single(WGW%COMM_INTER, F_KI, 2*VBMAX*VBMAX*WDES%NKPTS))
      CALLMPI(M_sum_single(WGW%COMM_INTER, F_BA, 2*(WDES%NB_TOT - VBMAX)*(WDES%NB_TOT - VBMAX)*WDES%NKPTS))

      CALLMPI(M_sum_z(WGW%COMM_INTER, W%CELTOT(1, 1, 1), WDES%NB_TOT*WDES%NKPTS*WDES%ISPIN))
      IF (CONSISTENT) WRITE (*, *) "Eigenvalues in WAVECAR and FOCKM agree. You &
           & are probably working with NATURAL ORBITALS. Everything seems to be OK!"

      IF (.not. CONSISTENT) THEN
         WRITE (*, *) "Eigenvalues in WAVECAR and FOCKM DISAGREE!! You &
          & should remove FOCKM if you are working with canonical orbitals."
      END IF

   END SUBROUTINE FOCKM_IN

   SUBROUTINE SETUP_DIIS(NKPTS, WDES)
      INTEGER :: NKPTS
      TYPE(wavedes) WDES

      NVECL = NUNOCC*NUNOCC*VBMAX*VBMAX*REALNKPTS*REALNKPTS*MY_NKPTS
      NVECLS = VBMAX*NUNOCC*REALNKPTS
      IF (LORBREAL) THEN
         ALLOCATE (MF_R(NVECL, DIISMAXSD), MMU_R(NVECL, DIISMAXSD))
         ALLOCATE (T2TMP_R(NVECL))
         MF_R = 0.0_q
         MMU_R = 0.0_q
      ELSE
         ALLOCATE (MF(NVECL, DIISMAXSD), MMU(NVECL, DIISMAXSD))
         ALLOCATE (T2TMP(NVECL))
         MF = zero
         MMU = zero
         IF (SINGLES) THEN
            ALLOCATE (MFS(NVECLS, DIISMAXSD), MMUS(NVECLS, DIISMAXSD))
            ALLOCATE (T1TMP(NVECLS))
            MFS = zero
            MMUS = zero
         END IF
      END IF
      IF (LBRUECKNER) THEN
         NVECLW = WDES%NRPLWV*WDES%NBANDS*WDES%NKPTS*WDES%ISPIN
         ALLOCATE (WF(NVECLW, DIISMAXSDW), WMUS(NVECLW, DIISMAXSDW))
         ALLOCATE (WTMP(NVECLW))
         WF = zero
         WMUS = zero
      END IF

   END SUBROUTINE

!***********************************************************************
!
! Setup plane wave vector grid used to represent electronic transition structure factor
!
!***********************************************************************

   SUBROUTINE SET_GVEC(WDES1, LATT_CUR, FSG, POTFAK, GVX, GVY, GVZ, ENCUT, ENCUTSOFT)
      USE constant
      USE wave
      USE lattice
      USE full_kpoints
      IMPLICIT NONE

      TYPE(wavedes1) WDES1
      TYPE(latt) :: LATT_CUR
      REAL(q)     :: FSG
      REAL(q) :: POTFAK(WDES1%NGVECTOR)
      REAL(q) :: GVX(WDES1%NGVECTOR)
      REAL(q) :: GVY(WDES1%NGVECTOR)
      REAL(q) :: GVZ(WDES1%NGVECTOR)
      REAL(q), OPTIONAL :: ENCUT, ENCUTSOFT
      ! local
      INTEGER NI, NP
      REAL(q) :: DKX, DKY, DKZ, GX, GY, GZ, GSQU, GSQUP, GQUAD, SCALE, SCALEFSG, RHOEFF, FOMEGAP2
      REAL(q) :: E, FSG_TMP
      REAL(q) :: Q1, Q2

      ! input needed for attenuated Coulomb kernel
      IF (PRESENT(ENCUT) .AND. PRESENT(ENCUTSOFT)) THEN
         Q1 = SQRT(ENCUTSOFT/HSQDTM)
         Q2 = SQRT(ENCUT/HSQDTM)
      END IF

      ! effective electron density in a.u.^-3
      RHOEFF = (HFSCREEN*HFSCREEN*AUTOA*AUTOA/(4._q*EXP(LOG(3._q/PI)/3._q)))**3
      ! plasma frequency
      FOMEGAP2 = 16._q*PI*RHOEFF/(AUTOA*AUTOA*AUTOA*AUTOA)

      SCALE = EDEPS/LATT_CUR%OMEGA/TPI**2

      IF (KPOINTS_FULL%WTKPT(1) == 0) THEN
         WRITE (*, *) 'internal error in  SET_GFAC_WAVEFUN: division by zero in SCALEFSG, and FSG anyway most likely wrong'
         WRITE (*, *) ' needs to be fixed in  CALCULATE_LOCAL_FIELD_PREPARE as well'
         STOP
      END IF
      SCALEFSG = 1/(KPOINTS_FULL%WTKPT(1)*NKREDX*NKREDY*NKREDZ)
      IF (ODDONLY .OR. EVENONLY) SCALEFSG = SCALEFSG/2

      DKX = WDES1%VKPT(1)*LATT_CUR%B(1, 1) + &
            WDES1%VKPT(2)*LATT_CUR%B(1, 2) + &
            WDES1%VKPT(3)*LATT_CUR%B(1, 3)
      DKY = WDES1%VKPT(1)*LATT_CUR%B(2, 1) + &
            WDES1%VKPT(2)*LATT_CUR%B(2, 2) + &
            WDES1%VKPT(3)*LATT_CUR%B(2, 3)
      DKZ = WDES1%VKPT(1)*LATT_CUR%B(3, 1) + &
            WDES1%VKPT(2)*LATT_CUR%B(3, 2) + &
            WDES1%VKPT(3)*LATT_CUR%B(3, 3)

      NP = WDES1%NGVECTOR

      DO NI = 1, NP

         GX = (WDES1%IGX(NI)*LATT_CUR%B(1, 1) + WDES1%IGY(NI)* &
               LATT_CUR%B(1, 2) + WDES1%IGZ(NI)*LATT_CUR%B(1, 3))
         GY = (WDES1%IGX(NI)*LATT_CUR%B(2, 1) + WDES1%IGY(NI)* &
               LATT_CUR%B(2, 2) + WDES1%IGZ(NI)*LATT_CUR%B(2, 3))
         GZ = (WDES1%IGX(NI)*LATT_CUR%B(3, 1) + WDES1%IGY(NI)* &
               LATT_CUR%B(3, 2) + WDES1%IGZ(NI)*LATT_CUR%B(3, 3))

         POTFAK(NI) = ((DKX + GX)**2 + (DKY + GY)**2 + (DKZ + GZ)**2)**0.5
         GVX(NI) = DKX + GX
         GVY(NI) = DKY + GY
         GVZ(NI) = DKZ + GZ

      END DO

   END SUBROUTINE SET_GVEC

   SUBROUTINE INTERPOLATE_REGULAR_SFAC(LATT_CUR,WDES,WGW,IOIU,ENCUTGW,ENCUTGWSOFT,ENERGY)
      USE constant
      USE wave
      USE lattice
      USE full_kpoints
!      USE bspline_module
      USE mathtools, ONLY : INVERT_REAL_MATRIX
      IMPLICIT NONE

      TYPE(latt) :: LATT_CUR
      TYPE(wavedes) WDES
      TYPE(wavedes) WGW
      INTEGER :: IOIU
      REAL(q) :: ENCUTGW, ENCUTGWSOFT
      ! local
      TYPE(wavedes1) WGWQ
      REAL(q) :: TINY = 1.E-7
      REAL(q) :: LARGE = 1.E+7
      INTEGER :: NG, KQ, NX, NY, NZ, I, IX, IY, IZ, KX, KY, KZ, IFLAG, IKNOT
      REAL(q), ALLOCATABLE :: TX(:), TY(:), TZ(:), BCOEF(:,:,:), W2(:,:),W1(:),W0(:)
      REAL(q) :: PROJA, PROJB, PROJC, AMIN, BMIN, CMIN, AMAX, BMAX,CMAX,LA,LB,LC
      REAL(q) :: val_int, XVAL, YVAL, ZVAL, GLEN, GX,GY,GZ, DX,DY,DZ
      REAL(q) :: REC_VEC(3,3), SCALE, POT, E, ENERGY, ENERGY_NOINT, GSQU, ESH
      INTEGER :: inbvx,inbvy,inbvz,iloy,iloz
      INTEGER :: NXFINE, NYFINE, NZFINE, IXF, IYF, IZF
      LOGICAL :: INIT, LEXTRAP


      SCALE = EDEPS/LATT_CUR%OMEGA/TPI**2

      REC_VEC(:,:)=LATT_CUR%B(:, :)
      CALL INVERT_REAL_MATRIX( REC_VEC, IOIU ) 

! determine lattice vectors and number of grid points
      LA=SQRT(LATT_CUR%B(1, 1)**2 + LATT_CUR%B(2, 1)**2 + LATT_CUR%B(3, 1)**2)
      LB=SQRT(LATT_CUR%B(1, 2)**2 + LATT_CUR%B(2, 2)**2 + LATT_CUR%B(3, 2)**2)
      LC=SQRT(LATT_CUR%B(1, 3)**2 + LATT_CUR%B(2, 3)**2 + LATT_CUR%B(3, 3)**2)

!      WRITE(*,*)'LATT VEC A', LATT_CUR%B(1, 1), LATT_CUR%B(2, 1), LATT_CUR%B(3, 1)
!      WRITE(*,*)'LATT VEC B', LATT_CUR%B(1, 2), LATT_CUR%B(2, 2), LATT_CUR%B(3, 2)
!      WRITE(*,*)'LATT VEC C', LATT_CUR%B(1, 3), LATT_CUR%B(2, 3), LATT_CUR%B(3, 3)

      AMIN=LARGE
      BMIN=LARGE
      CMIN=LARGE
      AMAX=TINY
      BMAX=TINY
      CMAX=TINY
      DO KQ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
         DO NG = 1, NGVECTOR
            IF ((GVECLEN(NG, RKQofKQ(KQ))) < TINY) CYCLE
!            WRITE(7,'(7G20.12)') GVECX(NG,RKQofKQ(KQ)),GVECY(NG,RKQofKQ(KQ)),GVECZ(NG,RKQofKQ(KQ)),GVECLEN(NG,RKQofKQ(KQ)),POTFAK_FULL(NG,RKQofKQ(KQ)),REAL(SDFACTOR(NG,RKQofKQ(KQ)),kind=q),REAL(SXFACTOR(NG,RKQofKQ(KQ)),kind=q)
             PROJA=(GVECX(NG,RKQofKQ(KQ))*REC_VEC(1, 1) + &
                             GVECY(NG,RKQofKQ(KQ))*REC_VEC(1, 2) + &
                             GVECZ(NG,RKQofKQ(KQ))*REC_VEC(1, 3))

             PROJB=(GVECX(NG,RKQofKQ(KQ))*REC_VEC(2, 1) + &
                             GVECY(NG,RKQofKQ(KQ))*REC_VEC(2, 2) + &
                             GVECZ(NG,RKQofKQ(KQ))*REC_VEC(2, 3))

             PROJC=(GVECX(NG,RKQofKQ(KQ))*REC_VEC(3, 1) + &
                             GVECY(NG,RKQofKQ(KQ))*REC_VEC(3, 2) + &
                             GVECZ(NG,RKQofKQ(KQ))*REC_VEC(3, 3))




             IF ((AMIN>PROJA) .AND. (PROJA>TINY)) AMIN=PROJA
             IF ((BMIN>PROJB) .AND. (PROJB>TINY)) BMIN=PROJB
             IF ((CMIN>PROJC) .AND. (PROJC>TINY)) CMIN=PROJC

             IF ((AMAX<PROJA)) AMAX=PROJA
             IF ((BMAX<PROJB)) BMAX=PROJB
             IF ((CMAX<PROJC)) CMAX=PROJC
         END DO
      END DO

!      WRITE(*,*)'MAX', AMAX, BMAX, CMAX
!      WRITE(*,*)'MIN', AMIN, BMIN, CMIN
!      WRITE(*,*)'IMAX', NINT(AMAX/AMIN), NINT(BMAX/BMIN), NINT(CMAX/CMIN)

      NX=NINT(AMAX/AMIN)*2+1
      NY=NINT(BMAX/BMIN)*2+1
      NZ=NINT(CMAX/CMIN)*2+1

      ALLOCATE(REG_SFACTOR(NX,NY,NZ))
!      ALLOCATE(REG_SXFACTOR(NX,NY,NZ))
      ALLOCATE(REG_X(NX))
      ALLOCATE(REG_Y(NY))
      ALLOCATE(REG_Z(NZ))

! determine regular grid points
      DO I=1,NX
         REG_X(I)=-AMAX+AMIN*(I-1)
      ENDDO
      DO I=1,NY
         REG_Y(I)=-BMAX+BMIN*(I-1)
      ENDDO
      DO I=1,NZ
         REG_Z(I)=-CMAX+CMIN*(I-1)
      ENDDO

!      WRITE(*,*)'REG_X', REG_X
!      WRITE(*,*)'REG_Y', REG_Y
!      WRITE(*,*)'REG_Z', REG_Z


! map structure factor to regular grid
      REG_SFACTOR=0.0_q
!      REG_SXFACTOR=0.0_q

      DO KQ = 1, WDES%NKPTS
         CALL SETWDES(WGW, WGWQ, KQ)
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
!         WRITE(*,*)'WGWQ%GRID%NPLWV=',WGWQ%GRID%NPLWV
         DO NG = 1, NGVECTOR
             PROJA=(GVECX(NG,RKQofKQ(KQ))*REC_VEC(1, 1) + &
                             GVECY(NG,RKQofKQ(KQ))*REC_VEC(1, 2) + &
                             GVECZ(NG,RKQofKQ(KQ))*REC_VEC(1, 3))

             PROJB=(GVECX(NG,RKQofKQ(KQ))*REC_VEC(2, 1) + &
                             GVECY(NG,RKQofKQ(KQ))*REC_VEC(2, 2) + &
                             GVECZ(NG,RKQofKQ(KQ))*REC_VEC(2, 3))

             PROJC=(GVECX(NG,RKQofKQ(KQ))*REC_VEC(3, 1) + &
                             GVECY(NG,RKQofKQ(KQ))*REC_VEC(3, 2) + &
                             GVECZ(NG,RKQofKQ(KQ))*REC_VEC(3, 3))

             IX=NINT((PROJA+AMAX)/AMIN)+1
             IY=NINT((PROJB+BMAX)/BMIN)+1
             IZ=NINT((PROJC+CMAX)/CMIN)+1

             REG_SFACTOR(IX,IY,IZ)=(2.0_q*REAL(SDFACTOR(NG,RKQofKQ(KQ)),kind=q)-REAL(SXFACTOR(NG,RKQofKQ(KQ)),kind=q))/WGWQ%GRID%NPLWV

!             WRITE(*,'(A,4F10.4)') 'ix,iy,iz,v: ',REG_X(ix),REG_Y(iy),REG_Z(iz), REG_SDFACTOR(IX,IY,IZ)
!             WRITE(*,*) 'ix,iy,iz,v: ',(ix),(iy),(iz), REG_SDFACTOR(IX,IY,IZ)
!             WRITE(*,'(A,4F10.4)') 'ix,iy,iz,v: ',(PROJA+AMAX)/AMIN+1,(PROJB+BMAX)/BMIN+1, (PROJC+CMAX)/CMIN+1,REAL(SDFACTOR(NG,RKQofKQ(KQ)),kind=q)

         END DO
      END DO


      KX=2
      KY=2
      KZ=2
      IKNOT=0

      ALLOCATE(TX(NX+KX))
      ALLOCATE(TY(NY+KY))
      ALLOCATE(TZ(NZ+KZ))
      ALLOCATE(BCOEF(NX,NY,NZ))

      inbvx=1
      inbvy=1
      inbvz=1
      iloy=1
      iloz=1

      ALLOCATE(W2(KY,KZ))
      ALLOCATE(W1(KZ))
      ALLOCATE(W0(3*MAX(KX,KY,KZ)))

!!      call db2ink(x,nx,y,ny,fcn_2d,kx,ky,iknot,tx,ty,fcn_2d,iflag)
!      CALL db3ink(REG_X,NX,REG_Y,NY,REG_Z,NZ,REG_SFACTOR,kx,ky,kz,iknot,tx,ty,tz,bcoef,iflag)
      if (iflag/=0) error stop 'error calling db3ink'

!     set paramters of finer grid for interpolation
!

      NXFINE=CCNINT
      NYFINE=CCNINT
      NZFINE=CCNINT
      LEXTRAP=.TRUE.

      DX=REG_X(2)-REG_X(1)
      DY=REG_Y(2)-REG_Y(1)
      DZ=REG_Z(2)-REG_Z(1)

      ENERGY=0.0_q
      ESH=ENCUTGW

      DO IX=1,NX
      DO IY=1,NY
      DO IZ=1,NZ

         DO IXF=1,NXFINE
         DO IYF=1,NYFINE
         DO IZF=1,NZFINE

            XVAL=REG_X(IX)-DX/2.0_q+DX/2.0_q/NXFINE+(IXF-1)*DX/NXFINE !+ABS(REG_X(1)-REG_X(2))*REAL(IXF,kind=q)/NXFINE/2.0_q
            YVAL=REG_Y(IY)-DY/2.0_q+DY/2.0_q/NYFINE+(IYF-1)*DY/NYFINE !+ABS(REG_X(1)-REG_X(2))*REAL(IXF,kind=q)/NXFINE/2.0_q
            ZVAL=REG_Z(IZ)-DZ/2.0_q+DZ/2.0_q/NZFINE+(IZF-1)*DZ/NZFINE !+ABS(REG_X(1)-REG_X(2))*REAL(IXF,kind=q)/NXFINE/2.0_q
!            YVAL=REG_Y(IY) !+ABS(REG_Y(1)-REG_Y(2))*REAL(IYF,kind=q)/NYFINE/2.0_q
!            ZVAL=REG_Z(IZ) !+ABS(REG_Z(1)-REG_Z(2))*REAL(IZF,kind=q)/NZFINE/2.0_q

!            CALL db3val(xval,yval,zval,0,0,0,tx,ty,tz,nx,ny,nz,kx,ky,kz,bcoef,val_int,iflag,inbvx,inbvy,inbvz,iloy,iloz,w2,w1,w0,lextrap)

            GX=XVAL*LATT_CUR%B(1, 1)+YVAL*LATT_CUR%B(1, 2)+ZVAL*LATT_CUR%B(1, 3)
            GY=XVAL*LATT_CUR%B(2, 1)+YVAL*LATT_CUR%B(2, 2)+ZVAL*LATT_CUR%B(2, 3)
            GZ=XVAL*LATT_CUR%B(3, 1)+YVAL*LATT_CUR%B(3, 2)+ZVAL*LATT_CUR%B(3, 3)
            !GY=XVAL*LATT_CUR%B(2, 1)+YVAL*LATT_CUR%B(2, 2)+ZVAL*LATT_CUR%B(2, 3)
            !GZ=XVAL*LATT_CUR%B(3, 1)+YVAL*LATT_CUR%B(3, 2)+ZVAL*LATT_CUR%B(3, 3)

!            WRITE(*,*)'gx,gy,gz ',XVAL, YVAL,ZVAL

            GLEN=SQRT(GX**2+GY**2+GZ**2)
            GSQU=GLEN**2
            E=HSQDTM*(GSQU*TPI**2)
            POT=SCALE/(GSQU)
            IF (E>ENCUTGWSOFT) THEN
               POT=POT*(1+COS((E-ENCUTGWSOFT)/(ENCUTGW-ENCUTGWSOFT)*PI))/2
            ENDIF
            IF (E>ENCUTGW) POT=0.0_q
!            IF (ISNAN(POT)) THEN
!             IF (E<=ENCUTGW) THEN
!                IF (IOIU >= 0) WRITE (*, *) 'nan', E,POT,GSQU,val_int !,ENCUTGWSOFT,ENCUTGW
!             ENDIF
!             IF (E<=ESH) THEN
!                ESH=E
!                IF (IOIU >= 0) WRITE (*, *) 'esh', ESH !,ENCUTGWSOFT,ENCUTGW
!             ENDIF
             IF (E<1.E-08) CYCLE
!            END IF
!            IF (ISNAN(val_int)) THEN
!               WRITE (*, *) 'isnan val_int', E,GSQU,ENCUTGWSOFT,ENCUTGW
!            END IF



            ENERGY=ENERGY+POT*val_int/(NXFINE*NYFINE*NZFINE)

!            WRITE(*,*)'glen,val ',GLEN,val_int !, REG_SDFACTOR(IX,IY,IZ)
            !WRITE(*,*)'glen,val ',GLEN,REG_SDFACTOR(IX,IY,IZ)


            if (iflag/=0) THEN
!               WRITE(*,*) 'iflag ',iflag,xval,yval,zval,GX,GY,GZ,GLEN
               !error stop 'error calling db3val'
            ENDIF
!            WRITE(*,*)'ix,iy,iz,v: ',ix,iy,iz,val_int
!         call db3val(xval,yval,zval,idx,idy,idz,tx,ty,tz,nx,ny,nz,kx,ky,kz,bcoef,f,iflag,inbvx,inbvy,inbvz,iloy,iloz,w2,w1,w0,extrap)


         ENDDO
         ENDDO
         ENDDO

      ENDDO
      ENDDO
      ENDDO

!      WRITE(*,*) 'energy from interpolated SFAC ',ENERGY

   END SUBROUTINE INTERPOLATE_REGULAR_SFAC



!***********************************************************************
!
! Add Coulomb kernel contribution to co-densities to get Coulomb-Vertex
!
!***********************************************************************

   SUBROUTINE APPLY_POTFAK_ALL(WDES)
      USE constant
      USE wave
      USE lattice
      USE full_kpoints
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER :: NG, KQ

      DO KQ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
         DO NG = 1, NGVECTOR
            FTOD_PW_IA(NG, :, :, RKQofKQ(KQ), :, 1) = FTOD_PW_IA(NG, :, :, RKQofKQ(KQ), :, 1)*POTFAK_FULL(NG, RKQofKQ(KQ))
            FTOD_PW_AI(NG, :, :, RKQofKQ(KQ), :, 1) = FTOD_PW_AI(NG, :, :, RKQofKQ(KQ), :, 1)*POTFAK_FULL(NG, RKQofKQ(KQ))
            FTOD_PW_AB(NG, :, :, RKQofKQ(KQ), :, 1) = FTOD_PW_AB(NG, :, :, RKQofKQ(KQ), :, 1)*POTFAK_FULL(NG, RKQofKQ(KQ))
            FTOD_PW_IJ(NG, :, :, RKQofKQ(KQ), :, 1) = FTOD_PW_IJ(NG, :, :, RKQofKQ(KQ), :, 1)*POTFAK_FULL(NG, RKQofKQ(KQ))
         END DO
      END DO

   END SUBROUTINE

!***********************************************************************
!
! This routine computes the electronic transition structure factor from
! a set of amplitudes and the Coulomb-Vertices.
!
!***********************************************************************

   SUBROUTINE CALC_SFACTOR(LATT_CUR, WGW, WDES, IO, ENCUTGW, ENCUTGWSOFT)
      USE constant
      USE wave
      USE lattice
      USE full_kpoints
      IMPLICIT NONE
      TYPE(latt) LATT_CUR
      TYPE(wavedes) WGW
      TYPE(wavedes) WDES
      TYPE(in_struct) IO
      REAL(q) :: ENCUTGW, ENCUTGWSOFT
! local variables
      INTEGER :: NG, KQ, NI, NJ, KI, KJ, NA, NB, KA, KB
      COMPLEX(q) :: ETEST
      REAL(q) ::  ENERGYINT

      IF (IO%IU0 >= 0) WRITE (*, *) 'Calculating T(G)'

      SDFACTOR = zero
      SXFACTOR = zero

      DO KB = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
         CALL BCAST2ALL_FTOD_AI(WDES, KB, 2, PW_AI_TMP, OC_AI_TMP)
         DO KA = 1, MY_NKPTS
         DO KQ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
            KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
            KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KB)) + &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)

            IF (.not. LORBREAL) THEN
               DO NG = 1, NGVECTOR
               DO NI = 1, VBMAX
               DO NJ = 1, VBMAX
               DO NA = 1, NUNOCC
               DO NB = 1, NUNOCC

                  IF (SINGLES .and. (KI == KPTS_MKPTS(KA)) .and. (KJ == KB)) THEN
                     SDFACTOR(NG, RKQofKQ(KQ)) = SDFACTOR(NG, RKQofKQ(KQ)) + (T1(NA, NI, RKIofKI(KI))*T1(NB, NJ, RKIofKI(KJ))) &
                                                 *KPOINTS_FULL%WTKPT(KI)*KPOINTS_FULL%WTKPT(KI)* &
                                        (FTOD_PW_AI_NOPOT(NG, NI, NA, RKQofKQ(KQ), KA, 1))*CONJG(PW_AI_TMP(NG, NJ, NB, RKQofKQ(KQ)))
                  END IF
  SDFACTOR(NG, RKQofKQ(KQ)) = SDFACTOR(NG, RKQofKQ(KQ)) + (T2(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA))*KPOINTS_FULL%WTKPT(1) &
                                              *KPOINTS_FULL%WTKPT(KI)*KPOINTS_FULL%WTKPT(KI)* &
                                        (FTOD_PW_AI_NOPOT(NG, NI, NA, RKQofKQ(KQ), KA, 1))*CONJG(PW_AI_TMP(NG, NJ, NB, RKQofKQ(KQ)))
               END DO
               END DO
               END DO
               END DO
               END DO
            ELSE
               DO NG = 1, NGVECTOR
               DO NI = 1, VBMAX
               DO NJ = 1, VBMAX
               DO NA = 1, NUNOCC
               DO NB = 1, NUNOCC

                  IF (SINGLES .and. (KI == KPTS_MKPTS(KA)) .and. (KJ == KB)) THEN
                     SDFACTOR(NG, KQ) = SDFACTOR(NG, KQ) + (T1_R(NA, NI, KI)*T1_R(NB, NJ, KJ)) &
                                        *KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1)* &
                                        (FTOD_PW_AI_NOPOT(NG, NI, NA, KQ, KA, 1))*CONJG(PW_AI_TMP(NG, NJ, NB, KQ))
                  END IF
                  SDFACTOR(NG, KQ) = SDFACTOR(NG, KQ) + (T2_R(NA, NB, NI, NJ, KI, KJ, KA))*KPOINTS_FULL%WTKPT(1) &
                                     *KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1)* &
                                     (FTOD_PW_AI_NOPOT(NG, NI, NA, KQ, KA, 1))*CONJG(PW_AI_TMP(NG, NJ, NB, KQ))
               END DO
               END DO
               END DO
               END DO
               END DO
            END IF

         END DO
         END DO
      END DO

      CALLMPI(M_sum_z(WDES%COMM, SDFACTOR, SIZE(SDFACTOR, 1)*SIZE(SDFACTOR, 2)))

      IF (IO%IU0 >= 0) WRITE (*, *) 'Calculating direct correlation energy ... '
      ETEST = zero
      DO KQ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
         DO NG = 1, NGVECTOR
            ETEST = ETEST + SDFACTOR(NG, RKQofKQ(KQ))*POTFAK_FULL(NG, RKQofKQ(KQ))
         END DO
      END DO
      IF (IO%IU0 >= 0) WRITE (*, *) 'Direct correlation energy: ', ETEST*2

      IF (IO%IU0 >= 0) WRITE (*, *) 'Calculating X(G)'

      SXFACTOR = zero

      DO KB = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
         CALL BCAST2ALL_FTOD_AI(WDES, KB, 2, PW_AI_TMP, OC_AI_TMP)
         DO KA = 1, MY_NKPTS
         DO KQ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0)) CYCLE
            KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
            KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) + &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)

            IF (.not. LORBREAL) THEN
               DO NG = 1, NGVECTOR
               DO NI = 1, VBMAX
               DO NJ = 1, VBMAX
               DO NA = 1, NUNOCC
               DO NB = 1, NUNOCC

                  IF (SINGLES .and. (KJ == KPTS_MKPTS(KA)) .and. (KI == KB)) THEN
                     SXFACTOR(NG, RKQofKQ(KQ)) = SXFACTOR(NG, RKQofKQ(KQ)) + (T1(NA, NJ, RKIofKI(KJ))*T1(NB, NI, RKIofKI(KI))) &
                                                 *KPOINTS_FULL%WTKPT(KI)*KPOINTS_FULL%WTKPT(KI)* &
                                        (FTOD_PW_AI_NOPOT(NG, NI, NA, RKQofKQ(KQ), KA, 1))*CONJG(PW_AI_TMP(NG, NJ, NB, RKQofKQ(KQ)))
                  END IF
 SXFACTOR(NG, RKQofKQ(KQ)) = SXFACTOR(NG, RKQofKQ(KQ)) + (T2(NA, NB, NJ, NI, RKIofKI(KJ), RKIofKI(KI), KA))*KPOINTS_FULL%WTKPT(KI) &
                                              *KPOINTS_FULL%WTKPT(KI)*KPOINTS_FULL%WTKPT(KI)* &
                                        (FTOD_PW_AI_NOPOT(NG, NI, NA, RKQofKQ(KQ), KA, 1))*CONJG(PW_AI_TMP(NG, NJ, NB, RKQofKQ(KQ)))
               END DO
               END DO
               END DO
               END DO
               END DO
            ELSE
               DO NG = 1, NGVECTOR
               DO NI = 1, VBMAX
               DO NJ = 1, VBMAX
               DO NA = 1, NUNOCC
               DO NB = 1, NUNOCC
                  IF (SINGLES .and. (KJ == KPTS_MKPTS(KA)) .and. (KI == KB)) THEN
                     SXFACTOR(NG, KQ) = SXFACTOR(NG, KQ) + (T1_R(NA, NJ, KJ)*T1_R(NB, NI, KI)) &
                                        *KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1)* &
                                        (FTOD_PW_AI_NOPOT(NG, NI, NA, KQ, KA, 1))*CONJG(PW_AI_TMP(NG, NJ, NB, KQ))
                  END IF

                  SXFACTOR(NG, KQ) = SXFACTOR(NG, KQ) + (T2_R(NA, NB, NJ, NI, KJ, KI, KA))*KPOINTS_FULL%WTKPT(1) &
                                     *KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1)* &
                                     (FTOD_PW_AI_NOPOT(NG, NI, NA, KQ, KA, 1))*CONJG(PW_AI_TMP(NG, NJ, NB, KQ))
               END DO
               END DO
               END DO
               END DO
               END DO
            END IF

         END DO
         END DO
      END DO

      CALLMPI(M_sum_z(WDES%COMM, SXFACTOR, SIZE(SXFACTOR, 1)*SIZE(SXFACTOR, 2)))

      IF (IO%IU0 >= 0) WRITE (*, *) 'Calculatung exchange-like correlation energy ... '
      ETEST = zero
      DO KQ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
         DO NG = 1, NGVECTOR
            ETEST = ETEST + SXFACTOR(NG, RKQofKQ(KQ))*POTFAK_FULL(NG, RKQofKQ(KQ))
         END DO
      END DO
      IF (IO%IU0 >= 0) WRITE (*, *) 'Exchange-like correlation energy: ', ETEST

      IF (IO%IU0 >= 0) THEN
         WRITE (*, *) 'Writing to CORRofG ...'
         OPEN (unit=7, file="CORRofG")
         WRITE (7, '(7G20.12)') 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
         DO KQ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
            DO NG = 1, NGVECTOR
               IF (GVECLEN(NG, RKQofKQ(KQ)) == 0.0_q) CYCLE
             WRITE(7,'(7G20.12)') GVECX(NG,RKQofKQ(KQ)),GVECY(NG,RKQofKQ(KQ)),GVECZ(NG,RKQofKQ(KQ)),GVECLEN(NG,RKQofKQ(KQ)),POTFAK_FULL(NG,RKQofKQ(KQ)),REAL(SDFACTOR(NG,RKQofKQ(KQ)),kind=q),REAL(SXFACTOR(NG,RKQofKQ(KQ)),kind=q)
            END DO
         END DO
         CLOSE (7)
      END IF

      CALL INTERPOLATE_REGULAR_SFAC(LATT_CUR,WDES,WGW,IO%IU0,ENCUTGW,ENCUTGWSOFT,ENERGYINT)
      E_CCSD_FS=ENERGYINT
      IF (IO%IU0 >= 0) WRITE (*, *) 'Interpolated correlation energy: ', ENERGYINT
      IF (IO%IU0 >= 0) WRITE (IO%IU6, *) 'Interpolated correlation energy: ', ENERGYINT
!      IF (IO%IU0 >= 0) WRITE (IO%IU6, *) 'CCSD correlation energy:', E_CCSD


   END SUBROUTINE

!***********************************************************************
!
! This routine computes the pair specific ppl correction to the approximate the
! CBS limit of the CCSD correlation energy
!
!***********************************************************************

   SUBROUTINE CALC_PSPPL(LATT_CUR, WGW, WDES, IO, ENCUTGW, ENCUTGWSOFT)
      USE constant
      USE wave
      USE lattice
      USE full_kpoints
      IMPLICIT NONE
      TYPE(latt) LATT_CUR
      TYPE(wavedes) WGW
      TYPE(wavedes) WDES
      TYPE(in_struct) IO
      REAL(q) :: ENCUTGW, ENCUTGWSOFT
! local variables
      INTEGER :: NG, KQ, NI, NJ, KI, KJ, NA, NB, KA, KB
      COMPLEX(q) :: ETEST, DPSPPL, DMP2
      REAL(q) ::  ENERGYINT

!D2_VVOO((NUNOC
!D2_OOOO(VBMAX,
!EMP2_PAIR_CBS(
!EMP2_PAIR_FBS(
!GMP2_PAIR(VBMA
!GCC_PAIR(VBMAX


      IF (IO%IU0 >= 0) WRITE (*, *) 'Delta_ijij'
      D2_OOOO= zero

      DO KJ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE
         CALL BCAST2ALL_FTOD_IJ(WDES, KJ, 2, PW_IJ_TMP, OC_IJ_TMP)
         DO KI = 1, MY_NKPTS
            KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KI)) - &
                                     WDES%VKPT(:, KPTS_MKPTS(KI)), KPOINTS_FULL)

!            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
!            KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
!                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
!            KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KB)) + &
!                                     WDES%VKPT(:, KQ), KPOINTS_FULL)

            DO NG = 1, NGVECTOR
            DO NI = 1, VBMAX
            DO NJ = 1, VBMAX

  D2_OOOO(NI,NJ,RKIofKI(KPTS_MKPTS(KI)), RKIofKI(KJ)) = D2_OOOO(NI,NJ,RKIofKI(KPTS_MKPTS(KI)), RKIofKI(KJ)) + &
                                             (FTOD_PW_IJ_NOPOT(NG, NI, NI, RKQofKQ(KQ), KI, 1))*CONJG(PW_IJ_TMP(NG, NJ, NJ, RKQofKQ(KQ)))
            END DO
            END DO
            END DO

         END DO
      END DO

      CALLMPI(M_sum_z(WDES%COMM, D2_OOOO, VBMAX**2*REALNKPTS**2 ))

      IF (IO%IU0 >= 0) WRITE (*, *) 'Computing gCC_ij'
      GCC_PAIR= zero

      DO KB = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
         CALL BCAST2ALL_FTOD_AI(WDES, KB, 2, PW_AI_TMP, OC_AI_TMP)
         DO KA = 1, MY_NKPTS
         DO KQ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
            KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
            KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KB)) + &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)

            DO NI = 1, VBMAX
            DO NJ = 1, VBMAX
            DO NA = 1, NUNOCC
            DO NB = 1, NUNOCC
!            DO NG = 1, NGVECTOR

!              IF (SINGLES .and. (KI == KPTS_MKPTS(KA)) .and. (KJ == KB)) THEN
!                 GCC_PAIR(NI,NJ,RKIofKI(KI),RKIofKI(KJ)) = &
!                         GCC_PAIR(NI,NJ,RKIofKI(KI),RKIofKI(KJ)) + &
!                         (T1(NA, NI, RKIofKI(KI))*T1(NB, NJ, RKIofKI(KJ))) &
!                         *(FTOD_PW_AI_NOPOT(NG, NI, NA,RKQofKQ(KQ), KA, 1))*CONJG(PW_AI_TMP(NG, NJ, NB, RKQofKQ(KQ))) &
!                         /D2_OOOO(NI,NJ,RKIofKI(KI), RKIofKI(KJ))
!
!              END IF
!                 GCC_PAIR(NI,NJ,RKIofKI(KI),RKIofKI(KJ)) = &
!                         GCC_PAIR(NI,NJ,RKIofKI(KI),RKIofKI(KJ)) + &
!                         (T2(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA))*KPOINTS_FULL%WTKPT(1) &
!                         *(FTOD_PW_AI_NOPOT(NG, NI, NA, RKQofKQ(KQ), KA, 1))*CONJG(PW_AI_TMP(NG, NJ, NB, RKQofKQ(KQ)))&
!                         /D2_OOOO(NI,NJ,RKIofKI(KI), RKIofKI(KJ))
!            END DO !NG


              IF (SINGLES .and. (KI == KPTS_MKPTS(KA)) .and. (KJ == KB)) THEN
                 GCC_PAIR(NI,NJ,RKIofKI(KI),RKIofKI(KJ)) = &
                         GCC_PAIR(NI,NJ,RKIofKI(KI),RKIofKI(KJ)) + &
                         (T1(NA, NI, RKIofKI(KI))*T1(NB, NJ, RKIofKI(KJ))) &
                         *(D2PAW_VVOO(NA,NB,NI,NJ,RKIofKI(KI), RKIofKI(KJ),KA))&
                         /D2PAW_OOOO(NI,NJ,RKIofKI(KI), RKIofKI(KJ))


              END IF
                 GCC_PAIR(NI,NJ,RKIofKI(KI),RKIofKI(KJ)) = &
                         GCC_PAIR(NI,NJ,RKIofKI(KI),RKIofKI(KJ)) + &
                         (T2(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA))*KPOINTS_FULL%WTKPT(1) &
                         *(D2PAW_VVOO(NA,NB,NI,NJ,RKIofKI(KI), RKIofKI(KJ),KA))&
                         /D2PAW_OOOO(NI,NJ,RKIofKI(KI), RKIofKI(KJ))

 


            END DO
            END DO
            END DO
            END DO

         END DO
         END DO
      END DO

      CALLMPI(M_sum_z(WDES%COMM, GCC_PAIR, VBMAX**2*REALNKPTS**2 ))



      IF (IO%IU0 >= 0) WRITE (*, *) 'Calculating ps-ppl correction ... '


      ETEST = zero
      DPSPPL = zero
      DMP2 = zero

      DO KJ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE
!         CALL BCAST2ALL_FTOD_IJ(WDES, KJ, 2, PW_IJ_TMP, OC_IJ_TMP)
         DO KI = 1, WDES%NKPTS !MY_NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE

!            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
!            KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
!                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
!            KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KB)) + &
!                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
            DO NI = 1, VBMAX
            DO NJ = 1, VBMAX

               DPSPPL = DPSPPL + ( EMP2_PAIR_CBS(NI,NJ,RKIofKI(KI), RKIofKI(KJ)) &
                             -  EMP2_PAIR_FBS(NI,NJ,RKIofKI(KI), RKIofKI(KJ)) ) &
                             *( GCC_PAIR(NI,NJ,RKIofKI(KI),RKIofKI(KJ))  + &
                                GMP2_PAIR(NI,NJ,RKIofKI(KI),RKIofKI(KJ)) + &
                                GMP2_PAIR(NI,NJ,RKIofKI(KI),RKIofKI(KJ))  &
                               *GCC_PAIR(NI,NJ,RKIofKI(KI),RKIofKI(KJ)) ) 
               DMP2 = DMP2 + ( EMP2_PAIR_CBS(NI,NJ,RKIofKI(KI), RKIofKI(KJ)) &
                             -  EMP2_PAIR_FBS(NI,NJ,RKIofKI(KI), RKIofKI(KJ)))
 
            END DO
            END DO

         END DO
      END DO

      IF (IO%IU0 >= 0) WRITE (*, *) 'ppl correction :', DPSPPL
      IF (IO%IU0 >= 0) WRITE (*, *) 'mp2 correction :', DMP2
      IF (IO%IU0 >= 0) WRITE (IO%IU6, *) 'Final CCSD correlation energy: ', REAL(E_CCSD,kind=q)
      IF (IO%IU0 >= 0) WRITE (IO%IU6, *) 'Final CCSD+ps-ppl correlation energy: ', REAL(E_CCSD+DPSPPL+DMP2,kind=q)
      IF (IO%IU0 >= 0) WRITE (IO%IU6, *) 'Final CCSD-FS correlation energy: ', REAL(E_CCSD_FS,kind=q)
      IF (IO%IU0 >= 0) WRITE (IO%IU6, *) 'Final CCSD-FS+ps-ppl correlation energy: ', REAL(E_CCSD_FS+DPSPPL+DMP2,kind=q)

      IF (IO%IU0 >= 0) WRITE (*, *) 'Final CCSD correlation energy: ', REAL(E_CCSD,kind=q)
      IF (IO%IU0 >= 0) WRITE (*, *) 'Final CCSD+ps-ppl correlation energy: ', REAL(E_CCSD+DPSPPL+DMP2,kind=q)
      IF (IO%IU0 >= 0) WRITE (*, *) 'Final CCSD-FS correlation energy: ', REAL(E_CCSD_FS,kind=q)
      IF (IO%IU0 >= 0) WRITE (*, *) 'Final CCSD-FS+ps-ppl correlation energy: ', REAL(E_CCSD_FS+DPSPPL+DMP2,kind=q)


   END SUBROUTINE

   SUBROUTINE CALC_GMP2_PAIR(LATT_CUR, WGW, WDES, IO, ENCUTGW, ENCUTGWSOFT)
      USE constant
      USE wave
      USE lattice
      USE full_kpoints
      IMPLICIT NONE
      TYPE(latt) LATT_CUR
      TYPE(wavedes) WGW
      TYPE(wavedes) WDES
      TYPE(in_struct) IO
      REAL(q) :: ENCUTGW, ENCUTGWSOFT
! local variables
      INTEGER :: NG, KQ, NI, NJ, KI, KJ, NA, NB, KA, KB
      COMPLEX(q) :: ETEST
      REAL(q) ::  ENERGYINT

      IF (IO%IU0 >= 0) WRITE (*, *) 'Preparing calculation of g(1)_ij (correlation depth scaling factor).'
      IF (IO%IU0 >= 0) WRITE (*, *) 'Computing Delta_ijij'
      D2_OOOO= zero

      DO KJ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE
         CALL BCAST2ALL_FTOD_IJ(WDES, KJ, 2, PW_IJ_TMP, OC_IJ_TMP)
         DO KI = 1, MY_NKPTS
            KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KI)) - &
                                     WDES%VKPT(:, KPTS_MKPTS(KI)), KPOINTS_FULL)

!            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
!            KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
!                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
!            KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KB)) + &
!                                     WDES%VKPT(:, KQ), KPOINTS_FULL)

            DO NI = 1, VBMAX
            DO NJ = 1, VBMAX

            DO NG = 1, NGVECTOR
               D2_OOOO(NI,NJ,RKIofKI(KPTS_MKPTS(KI)), RKIofKI(KJ)) = D2_OOOO(NI,NJ,RKIofKI(KPTS_MKPTS(KI)), RKIofKI(KJ)) + &
                                             (FTOD_PW_IJ_NOPOT(NG, NI, NI, RKQofKQ(KQ), KI, 1))*CONJG(PW_IJ_TMP(NG, NJ, NJ, RKQofKQ(KQ)))
            END DO


!            WRITE(*,*)'D2_OOOO(',NI,NJ,')=',D2_OOOO(NI,NJ,RKIofKI(KPTS_MKPTS(KI)),RKIofKI(KJ))
            END DO
            END DO

         END DO
      END DO

      CALLMPI(M_sum_z(WDES%COMM, D2_OOOO, VBMAX**2*REALNKPTS**2 ))

      IF (IO%IU0 >= 0) WRITE (*, *) 'Computing  g(1)_ij'
      GMP2_PAIR= zero

      DO KB = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
         CALL BCAST2ALL_FTOD_AI(WDES, KB, 2, PW_AI_TMP, OC_AI_TMP)
         DO KA = 1, MY_NKPTS
         DO KQ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
            KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
            KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, (KB)) + &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)

            DO NI = 1, VBMAX
            DO NJ = 1, VBMAX
            DO NA = 1, NUNOCC
            DO NB = 1, NUNOCC
!            DO NG = 1, NGVECTOR
!                 GMP2_PAIR(NI,NJ,RKIofKI(KI),RKIofKI(KJ)) = &
!                         GMP2_PAIR(NI,NJ,RKIofKI(KI),RKIofKI(KJ)) + &
!                         (T2(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA))*KPOINTS_FULL%WTKPT(1) &
!                         *(FTOD_PW_AI_NOPOT(NG, NI, NA, RKQofKQ(KQ), KA, 1))*CONJG(PW_AI_TMP(NG, NJ, NB, RKQofKQ(KQ)))&
!                         /D2_OOOO(NI,NJ,RKIofKI(KI), RKIofKI(KJ))
!            END DO

                 GMP2_PAIR(NI,NJ,RKIofKI(KI),RKIofKI(KJ)) = &
                         GMP2_PAIR(NI,NJ,RKIofKI(KI),RKIofKI(KJ)) + &
                         (T2(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA))*KPOINTS_FULL%WTKPT(1) &
                         *(D2PAW_VVOO(NA,NB,NI,NJ,RKIofKI(KI), RKIofKI(KJ), KA))&
                         /D2PAW_OOOO(NI,NJ,RKIofKI(KI), RKIofKI(KJ))

            END DO
            END DO
            END DO
            END DO

         END DO
         END DO
      END DO

      CALLMPI(M_sum_z(WDES%COMM, GMP2_PAIR, VBMAX**2*REALNKPTS**2 ))

      IF (IO%IU0 >= 0) WRITE (*, *) 'Done.'
      IF (IO%IU0 >= 0) WRITE (*, *) 'GMP2_PAIR:', GMP2_PAIR


   END SUBROUTINE

   SUBROUTINE TEST_D2PAW(LATT_CUR, W, WDES, T_INFO, P, IO, ENCUTGW, ENCUTGWSOFT)
      USE constant
      USE wave
      USE pseudo
      USE lattice
      USE full_kpoints
      USE asa
      USE radial
      IMPLICIT NONE
      TYPE(latt) LATT_CUR
      TYPE(wavespin) W
      TYPE(wavedes) WDES
      TYPE(type_info) T_INFO
      TYPE(potcar) P(T_INFO%NTYP)
      TYPE(in_struct) IO
      REAL(q) :: ENCUTGW, ENCUTGWSOFT
! local variables
      INTEGER :: NG, KQ, NI, NJ, KI, KJ, NA, NB, KA, KB, ISP, KI_IN_FULL_ORIG
! for D2PAW calc
      INTEGER :: NSTRIP, I
      COMPLEX(q) :: INTC
      TYPE(wavespin) WHF
      TYPE(wavedes1) WDESKI, WDESKJ, WDESKA, WDESKB
      TYPE(wavefun1), ALLOCATABLE :: WI(:), WJ(:), WA(:), WB(:)
! for D2PAW one-center terms
      INTEGER :: CH1, CH2, CH3, CH4, ION, LMMAX, SPECIES, NBI, NBJ, NT
      INTEGER :: POSCH1, POSCH2, POSCH3, POSCH4
      INTEGER :: L0P, L1P, L2P, L3P, M0P, M1P, M2P, M3P, M0PMAX, M1PMAX, M2PMAX, M3PMAX
      COMPLEX(q), ALLOCATABLE :: DELTA_PAW(:,:,:,:,:) ! DELTA_PAW(CH1,CH2,CH3,CH4,ION_TYPE) delta-kernel in PAW channel basis for each ion type 
      COMPLEX(q), ALLOCATABLE :: DELTA_PAW_RAD(:,:,:,:,:) ! DELTA_PAW(CH1,CH2,CH3,CH4,ION_TYPE) delta-kernel in PAW channel basis for radial factor contribition
      COMPLEX(q), ALLOCATABLE :: DELTA_PAW_SPH(:,:,:,:,:) ! DELTA_PAW(CH1,CH2,CH3,CH4,ION_TYPE) delta-kernel in PAW channel basis for angular factor contribition
      REAL(q) :: QRHO
      REAL(q), ALLOCATABLE :: RHO(:)
!For computing YLM3
      INTEGER :: LM01INDX, ISTART01, IEND01, LM23INDX, ISTART23, IEND23, IC01, IC23, NPRO

      ISP=1

! compute D2PAW_OOOO & D2PAW_VVOO including PAW terms properly

!***********************************************************************
! start with D2PAW_OOOO
!***********************************************************************
      IF (IO%IU0 >= 0) WRITE (*, *) 'Preparing calculation of g(1)_ij (correlation depth scaling factor).'


!***********************************************************************
!compute four-channel--one-center--delta-kernel for each atomic species
! THIS COULD GO IN A SEPARATE ROUTINE
!***********************************************************************
     !find lmmax of all ions
     LMMAX=0
     DO SPECIES=1,T_INFO%NTYP
        IF (P(SPECIES)%LMMAX > LMMAX) WRITE(*,*)'found lmmax for ion type:',SPECIES, P(SPECIES)%LMMAX
        IF (P(SPECIES)%LMMAX > LMMAX) LMMAX=P(SPECIES)%LMMAX
     ENDDO
     !lmmax should be identical to LMDIM. However, could also be organised different. Compare to CPROJ and CQIJ arrays.

     !allocate and compute DELTA_PAW
     ALLOCATE(DELTA_PAW(LMMAX,LMMAX,LMMAX,LMMAX,T_INFO%NTYP))
     ALLOCATE(DELTA_PAW_RAD(LMMAX,LMMAX,LMMAX,LMMAX,T_INFO%NTYP))
     ALLOCATE(DELTA_PAW_SPH(LMMAX,LMMAX,LMMAX,LMMAX,T_INFO%NTYP))
     DELTA_PAW=zero
     DELTA_PAW_RAD=zero
     DELTA_PAW_SPH=zero

     DO SPECIES=1,T_INFO%NTYP

        IF (P(SPECIES)%LMAX==0) EXIT

     !for each ion type compute delta(r1-r2) kernel integral of pseudo and ae partial wave contributions
     !
     !this algorithm is not efficient but will not be a bottle neck

        IF (ALLOCATED(RHO)) DEALLOCATE(RHO)
        ALLOCATE(RHO(P(SPECIES)%R%NMAX))
        RHO=0._q

        !P(I) or PP ?
        POSCH1=1
        DO CH1=1,P(SPECIES)%LMAX
        !get L0P and M0P
        L0P =P(SPECIES)%LPS(CH1)
           WRITE(*,*)'L0P', L0P

        POSCH2=1
        DO CH2=1,P(SPECIES)%LMAX
        !get L1P and M1P
        L1P =P(SPECIES)%LPS(CH2)

        POSCH3=1
        CH3=1
        L2P =P(SPECIES)%LPS(CH3)

        POSCH4=1
        CH4=1
           !get L3P and M3P
        L3P =P(SPECIES)%LPS(CH4)



              IF (L0P>0) CYCLE
              IF (L1P>0) CYCLE
              IF (L2P>0) CYCLE
              IF (L3P>0) CYCLE

           !compute factor of delta-kernel contribution coming from radial coordinate integration 
           !
           !(not efficient, work is duplicated
           !for all values of m,m',m'' and m''' quantum numbers).
           RHO=0._q
           QRHO=0._q
!           CALL SET_SIMP(P(SPECIES)%R)
           DO I=1,P(SPECIES)%R%NMAX
!              IF (P(SPECIES)%R%R(I)>T_INFO%RWIGS(NT)) THEN
!                 WRITE(*,*)'P(SPECIES)%R%R(I)>T_INFO%RWIGS(NT)',P(SPECIES)%R%R(I),T_INFO%RWIGS(NT)
!                 EXIT
!              ENDIF
              RHO(I)=RHO(I)+P(SPECIES)%WAE(I,CH1)*P(SPECIES)%WAE(I,CH2) 
!              RHO(I)=RHO(I)-P(SPECIES)%WPS(I,CH1)*P(SPECIES)%WPS(I,CH2)

           ENDDO
           CALL SIMPI(P(SPECIES)%R,RHO,QRHO)
!           WRITE(*,*)'L0P,L1P,L2P,L3P,',L0P,L1P,L2P,L3P,QRHO

           M0PMAX=2*L0P+1
           M1PMAX=2*L1P+1
           M2PMAX=2*L2P+1
           M3PMAX=2*L3P+1

           DO M0P=0,M0PMAX-1
           DO M1P=0,M1PMAX-1
           DO M2P=0,M2PMAX-1
           DO M3P=0,M3PMAX-1
              IF (L0P>0) CYCLE
              IF (L1P>0) CYCLE
              IF (L2P>0) CYCLE
              IF (L3P>0) CYCLE
              DELTA_PAW_RAD(POSCH1+M0P,POSCH2+M1P,POSCH3+M2P,POSCH4+M3P,SPECIES)=QRHO
           ENDDO
           ENDDO
           ENDDO
           ENDDO

     
           !next we need to compite DELTA_PAW_SPH =  <l,m.l',m'| L,M><L,M | l'',m'',l''',m'''>
           !
           ! For this we ńeed the  YLM3 Coefficients from the asa module
           ! The YLM3 coeffs correspond to the integrals of three spherical harmonics using l,m,l',m' and L,M
           ! (see similar routines such as SETUP_TRANS_MATRIX in fast_aug.F)
           !
           
           ! First look up <l,m.l',m'| L,M>
           !   ... then  <L,M | l'',m'',l''',m'''>
           ! Then contract over L,M and store in DELTA_PAW_SPH


           DO M0P=0,M0PMAX-1
           DO M1P=0,M1PMAX-1
           DO M2P=0,M2PMAX-1
           DO M3P=0,M3PMAX-1
              DELTA_PAW(POSCH1+M0P,POSCH2+M1P,POSCH3+M2P,POSCH4+M3P,SPECIES)=&
                      DELTA_PAW_RAD(POSCH1+M0P,POSCH2+M1P,POSCH3+M2P,POSCH4+M3P,SPECIES) 
! &
!                     *DELTA_PAW_SPH(POSCH1+M0P,POSCH2+M1P,POSCH3+M2P,POSCH4+M3P,SPECIES)

           WRITE(*,*)'CH1,CH2,CH3,CH4,DELTA_PAW_SPH',CH1,CH2,CH3,CH4,DELTA_PAW_SPH(POSCH1+M0P,POSCH2+M1P,POSCH3+M2P,POSCH4+M3P,SPECIES)
           ENDDO
           ENDDO
           ENDDO
           ENDDO




        POSCH4=POSCH4+2*L3P+1
!        ENDDO !ch4
        POSCH3=POSCH3+2*L2P+1
!        ENDDO !ch3
        POSCH2=POSCH2+2*L1P+1
        ENDDO !ch2
        POSCH1=POSCH1+2*L0P+1
        ENDDO !ch1
     

     ENDDO



!***********************************************************************
!First add pseudo part contribution on real space grid   < ij | delta(r1-r2) | ij >
!***********************************************************************
      IF (IO%IU0 >= 0) WRITE (*, *) 'Computing pseudo part of Delta_ijij'
      D2PAW_OOOO= zero


      WHF = W
      WHF%WDES => WDES_FOCK
      CALL SETWDES(WHF%WDES, WDESKI, 0)
      CALL SETWDES(WHF%WDES, WDESKJ, 0)


      ALLOCATE (WI(VBMAX) )
      ALLOCATE (WJ(VBMAX) )
      DO NBI = 1, VBMAX
         CALL NEWWAV(WI(NBI), WDESKI, .TRUE.)
         CALL NEWWAV(WJ(NBI), WDESKJ, .TRUE.)
      END DO


      DO KI = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE

         KI_IN_FULL_ORIG = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:, KI), KPOINTS_FULL_ORIG)
         CALL SETWDES(WHF%WDES, WDESKI, KI)
         CALL W1_GATHER_GLB(WHF, 1, VBMAX, ISP, WI)
               ! k_b = k_i - k_q - G
!               KB = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KQ) + WDES%VKPT(:, KI), KPOINTS_FULL)
!               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) Write (*, *) 'error in calc_ftod with shifted k-mesh'
               ! k_a = k_i + k_q - G
!               KA = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - WDES%VKPT(:, KQ), KPOINTS_FULL)
!               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA) == 0)) Write (*, *) 'error in calc_ftod with shifted k-mesh'
      WRITE(*,*)'size of wi',WI(1)%WDES1%GRID%MPLWV
      DO KJ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE

         CALL SETWDES(WHF%WDES, WDESKJ, KJ)
         CALL W1_GATHER_GLB(WHF, 1, VBMAX, ISP, WJ)


         DO NBI = 1, VBMAX
         DO NBJ = 1, VBMAX

            INTC=zero
            DO I=1,WI(NBI)%WDES1%GRID%MPLWV
!               INTC=INTC+(CONJG(WI(NBI)%CR(I))*WJ(NBJ)%CR(I))
            ENDDO
            D2PAW_OOOO(NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)), RKIofKI(KJ))=INTC/WI(NBI)%WDES1%GRID%MPLWV

            WRITE(*,*)'D2PAW_OOOO(',NBI,NBJ,')=',D2PAW_OOOO(NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)),RKIofKI(KJ))
         ENDDO
         ENDDO

     ENDDO !kj loop
     ENDDO !ki loop
 
!***********************************************************************
!compute delta-kernel-contributions of PAW terms from all atoms to D2PAW_OOOO for each state
!***********************************************************************
      IF (IO%IU0 >= 0) WRITE (*, *) 'Computing one-center part of Delta_ijij'

      DO KI = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE

         KI_IN_FULL_ORIG = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:, KI), KPOINTS_FULL_ORIG)
         CALL SETWDES(WHF%WDES, WDESKI, KI)
         CALL W1_GATHER_GLB(WHF, 1, VBMAX, ISP, WI)
      DO KJ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE

         CALL SETWDES(WHF%WDES, WDESKJ, KJ)
         CALL W1_GATHER_GLB(WHF, 1, VBMAX, ISP, WJ)


         DO NBI = 1, VBMAX
         DO NBJ = 1, VBMAX

            ! loop over atoms
            ! select corresponding species and loop over channels
            DO NI=1,W%WDES%NIONS
               SPECIES=W%WDES%ITYP(NI)
               LMMAX=W%WDES%LMMAX(SPECIES)
               IF (LMMAX==0) CYCLE
               NPRO =W%WDES%LMBASE(NI)
!               NPRO2_=WDES1%LMBASE(NI)+NPRO_


               !P(I) or PP ?
               DO CH1=1,LMMAX
                  WRITE(*,*)'ch1,nbi,npro,cproj',CH1,nbi,npro,WI(NBI)%CPROJ(CH1+NPRO)
               DO CH2=1,LMMAX
!               DO CH3=1,LMMAX
!               DO CH4=1,LMMAX

                  D2PAW_OOOO(NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)), RKIofKI(KJ))= D2PAW_OOOO(NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)), RKIofKI(KJ)) + &
                          CONJG(W%CPROJ(CH1+NPRO,NBI,1,1))*W%CPROJ(CH2+NPRO,NBJ,1,1)* &
                          DELTA_PAW(CH1,CH2,1,1,SPECIES)

!                  D2PAW_OOOO(NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)), RKIofKI(KJ))= D2PAW_OOOO(NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)), RKIofKI(KJ)) + &
!                          (CONJG(WI(NBI)%CPROJ(CH1+NPRO)*WJ(NBJ)%CPROJ(CH2+NPRO))*WI(NBI)%CPROJ(CH3+NPRO)*WJ(NBJ)%CPROJ(CH4+NPRO)) * &
!                          DELTA_PAW(CH1,CH2,CH3,CH4,SPECIES)



!               ENDDO !ch4
!               ENDDO !ch3
               ENDDO !ch2
               ENDDO !ch1
 
            ENDDO !NI

            WRITE(*,*)'D2PAW_OOOO(',NBI,NBJ,')=',D2PAW_OOOO(NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)),RKIofKI(KJ))
         ENDDO !nbi
         ENDDO !nbj

     ENDDO !kj loop
     ENDDO !ki loop
 



   END SUBROUTINE



   SUBROUTINE CALC_D2PAW(LATT_CUR, W, WDES, T_INFO, P, IO, ENCUTGW, ENCUTGWSOFT)
      USE constant
      USE wave
      USE pseudo
      USE lattice
      USE full_kpoints
      USE asa
      USE radial
      IMPLICIT NONE
      TYPE(latt) LATT_CUR
      TYPE(wavespin) W
      TYPE(wavedes) WDES
      TYPE(type_info) T_INFO
      TYPE(potcar) P(T_INFO%NTYP)
      TYPE(in_struct) IO
      REAL(q) :: ENCUTGW, ENCUTGWSOFT
! local variables
      INTEGER :: NG, KQ, NI, NJ, KI, KJ, NA, NB, KA, KB, ISP, KI_IN_FULL_ORIG
! for D2PAW calc
      INTEGER :: NSTRIP, I
      COMPLEX(q) :: INTC
      TYPE(wavespin) WHF
      TYPE(wavedes1) WDESKI, WDESKJ, WDESKA, WDESKB
      TYPE(wavefun1), ALLOCATABLE :: WI(:), WJ(:), WA(:), WB(:)
! for D2PAW one-center terms
      INTEGER :: CH1, CH2, CH3, CH4, ION, LMMAX, SPECIES, NBI, NBJ, NT
      INTEGER :: POSCH1, POSCH2, POSCH3, POSCH4
      INTEGER :: L0P, L1P, L2P, L3P, M0P, M1P, M2P, M3P, M0PMAX, M1PMAX, M2PMAX, M3PMAX
      COMPLEX(q), ALLOCATABLE :: DELTA_PAW(:,:,:,:,:) ! DELTA_PAW(CH1,CH2,CH3,CH4,ION_TYPE) delta-kernel in PAW channel basis for each ion type 
      COMPLEX(q), ALLOCATABLE :: DELTA_PAW_RAD(:,:,:,:,:) ! DELTA_PAW(CH1,CH2,CH3,CH4,ION_TYPE) delta-kernel in PAW channel basis for radial factor contribition
      COMPLEX(q), ALLOCATABLE :: DELTA_PAW_SPH(:,:,:,:,:) ! DELTA_PAW(CH1,CH2,CH3,CH4,ION_TYPE) delta-kernel in PAW channel basis for angular factor contribition
      REAL(q) :: QRHO
      REAL(q), ALLOCATABLE :: RHO(:)
!For computing YLM3
      INTEGER :: LM01INDX, ISTART01, IEND01, LM23INDX, ISTART23, IEND23, IC01, IC23, NPRO
!  D2PAW_VVOO related variables
      INTEGER :: NSTRIPA, NSTRIPB, NBA, NBB, NBAP, NBBP
      COMPLEX(q), ALLOCATABLE :: D2PAW_INT(:,:,:) ! intermediate quantity for fast calculation of one-center terms
! k-point related quantities
      COMPLEX(q), ALLOCATABLE :: CPHASE(:)
      LOGICAL :: LPHASE

      ISP=1

! compute D2PAW_OOOO & D2PAW_VVOO including PAW terms properly

!***********************************************************************
! start with D2PAW_OOOO
!***********************************************************************
      IF (IO%IU0 >= 0) WRITE (*, *) 'Preparing calculation of g(1)_ij (correlation depth scaling factor).'


!***********************************************************************
!compute four-channel--one-center--delta-kernel for each atomic species
! THIS COULD GO IN A SEPARATE ROUTINE
!***********************************************************************
     !find lmmax of all ions
     LMMAX=0
     DO SPECIES=1,T_INFO%NTYP
!        IF (P(SPECIES)%LMMAX > LMMAX) WRITE(*,*)'found lmmax for ion type:',SPECIES, P(SPECIES)%LMMAX
        IF (P(SPECIES)%LMMAX > LMMAX) LMMAX=P(SPECIES)%LMMAX
     ENDDO
     !lmmax should be identical to LMDIM. However, could also be organised different. Compare to CPROJ and CQIJ arrays.

     !allocate and compute DELTA_PAW
     ALLOCATE(DELTA_PAW(LMMAX,LMMAX,LMMAX,LMMAX,T_INFO%NTYP))
     ALLOCATE(DELTA_PAW_RAD(LMMAX,LMMAX,LMMAX,LMMAX,T_INFO%NTYP))
     ALLOCATE(DELTA_PAW_SPH(LMMAX,LMMAX,LMMAX,LMMAX,T_INFO%NTYP))
     ALLOCATE(D2PAW_INT(LMMAX,LMMAX,W%WDES%NIONS))

     DELTA_PAW=zero
     DELTA_PAW_RAD=zero
     DELTA_PAW_SPH=zero

     DO SPECIES=1,T_INFO%NTYP

        IF (P(SPECIES)%LMAX==0) EXIT

     !for each ion type compute delta(r1-r2) kernel integral of pseudo and ae partial wave contributions
     !
     !this algorithm is not efficient but will not be a bottle neck

        IF (ALLOCATED(RHO)) DEALLOCATE(RHO)
        ALLOCATE(RHO(P(SPECIES)%R%NMAX))
        RHO=0._q

        !P(I) or PP ?
        POSCH1=1
        DO CH1=1,P(SPECIES)%LMAX
        !get L0P and M0P
        L0P =P(SPECIES)%LPS(CH1)
!           WRITE(*,*)'L0P', L0P

        POSCH2=1
        DO CH2=1,P(SPECIES)%LMAX
        !get L1P and M1P
        L1P =P(SPECIES)%LPS(CH2)

        POSCH3=1
        DO CH3=1,P(SPECIES)%LMAX
        !get L2P and M2P
        L2P =P(SPECIES)%LPS(CH3)

        POSCH4=1
        DO CH4=1,P(SPECIES)%LMAX
           !get L3P and M3P
           L3P =P(SPECIES)%LPS(CH4)



!              IF (L0P>0) CYCLE
!              IF (L1P>0) CYCLE
!              IF (L2P>0) CYCLE
!              IF (L3P>0) CYCLE

           !compute factor of delta-kernel contribution coming from radial coordinate integration 
           !
           !(not efficient, work is duplicated
           !for all values of m,m',m'' and m''' quantum numbers).
           RHO=0._q
           QRHO=0._q
!           CALL SET_SIMP(P(SPECIES)%R)
           DO I=1,P(SPECIES)%R%NMAX
!              IF (P(SPECIES)%R%R(I)>T_INFO%RWIGS(NT)) THEN
!                 WRITE(*,*)'P(SPECIES)%R%R(I)>T_INFO%RWIGS(NT)',P(SPECIES)%R%R(I),T_INFO%RWIGS(NT)
!                 EXIT
!              ENDIF
              RHO(I)=RHO(I)+(P(SPECIES)%WAE(I,CH1)*P(SPECIES)%WAE(I,CH2)*P(SPECIES)%WAE(I,CH3)*P(SPECIES)%WAE(I,CH4))/P(SPECIES)%R%R(I)**2*4*PI
              RHO(I)=RHO(I)-(P(SPECIES)%WAE(I,CH1)*P(SPECIES)%WAE(I,CH2)*P(SPECIES)%WPS(I,CH3)*P(SPECIES)%WPS(I,CH4))/P(SPECIES)%R%R(I)**2*4*PI
              RHO(I)=RHO(I)-(P(SPECIES)%WPS(I,CH1)*P(SPECIES)%WPS(I,CH2)*P(SPECIES)%WAE(I,CH3)*P(SPECIES)%WAE(I,CH4))/P(SPECIES)%R%R(I)**2*4*PI
              RHO(I)=RHO(I)+(P(SPECIES)%WPS(I,CH1)*P(SPECIES)%WPS(I,CH2)*P(SPECIES)%WPS(I,CH3)*P(SPECIES)%WPS(I,CH4))/P(SPECIES)%R%R(I)**2*4*PI

           ENDDO
           CALL SIMPI(P(SPECIES)%R,RHO,QRHO)
!           WRITE(*,*)'L0P,L1P,L2P,L3P,',L0P,L1P,L2P,L3P,QRHO

           M0PMAX=2*L0P+1
           M1PMAX=2*L1P+1
           M2PMAX=2*L2P+1
           M3PMAX=2*L3P+1

           DO M0P=0,M0PMAX-1
           DO M1P=0,M1PMAX-1
           DO M2P=0,M2PMAX-1
           DO M3P=0,M3PMAX-1
!              IF (L0P>0) CYCLE
!              IF (L1P>0) CYCLE
!              IF (L2P>0) CYCLE
!              IF (L3P>0) CYCLE
              DELTA_PAW_RAD(POSCH1+M0P,POSCH2+M1P,POSCH3+M2P,POSCH4+M3P,SPECIES)=QRHO
           ENDDO
           ENDDO
           ENDDO
           ENDDO

     
           !next we need to compite DELTA_PAW_SPH =  <l,m.l',m'| L,M><L,M | l'',m'',l''',m'''>
           !
           ! For this we ńeed the  YLM3 Coefficients from the asa module
           ! The YLM3 coeffs correspond to the integrals of three spherical harmonics using l,m,l',m' and L,M
           ! (see similar routines such as SETUP_TRANS_MATRIX in fast_aug.F)
           !
           
           ! First look up <l,m.l',m'| L,M>
           !   ... then  <L,M | l'',m'',l''',m'''>
           ! Then contract over L,M and store in DELTA_PAW_SPH

           CALL YLM3LOOKUP(L0P,L1P,LM01INDX)

           DO M0P=0,M0PMAX-1
           DO M1P=0,M1PMAX-1
              LM01INDX=LM01INDX+1

              ISTART01=INDCG(LM01INDX)
              IEND01  =INDCG(LM01INDX+1)

              CALL YLM3LOOKUP(L2P,L3P,LM23INDX)

              DO M2P=0,M2PMAX-1
              DO M3P=0,M3PMAX-1
 
                 LM23INDX=LM23INDX+1

                 ISTART23=INDCG(LM23INDX)
                 IEND23  =INDCG(LM23INDX+1)

                 DO IC01=ISTART01,IEND01-1
                 DO IC23=ISTART23,IEND23-1

                    IF (JL(IC01)==JL(IC23)) THEN
                    IF (JS(IC01)==JS(IC23)) THEN
                       DELTA_PAW_SPH(POSCH1+M0P,POSCH2+M1P,POSCH3+M2P,POSCH4+M3P,SPECIES)= &
                          DELTA_PAW_SPH(POSCH1+M0P,POSCH2+M1P,POSCH3+M2P,POSCH4+M3P,SPECIES) + YLM3(IC01)*YLM3(IC23)
!*1.2732395447351628
                            
                    ENDIF
                    ENDIF
                 ENDDO
                 ENDDO


              ENDDO
              ENDDO
           ENDDO
           ENDDO

           DO M0P=0,M0PMAX-1
           DO M1P=0,M1PMAX-1
           DO M2P=0,M2PMAX-1
           DO M3P=0,M3PMAX-1
              DELTA_PAW(POSCH1+M0P,POSCH2+M1P,POSCH3+M2P,POSCH4+M3P,SPECIES)=&
                      DELTA_PAW_RAD(POSCH1+M0P,POSCH2+M1P,POSCH3+M2P,POSCH4+M3P,SPECIES) &
                     *DELTA_PAW_SPH(POSCH1+M0P,POSCH2+M1P,POSCH3+M2P,POSCH4+M3P,SPECIES)

!           WRITE(*,*)'CH1,CH2,CH3,CH4,DELTA_PAW_SPH',CH1,CH2,CH3,CH4,DELTA_PAW_SPH(POSCH1+M0P,POSCH2+M1P,POSCH3+M2P,POSCH4+M3P,SPECIES)
           ENDDO
           ENDDO
           ENDDO
           ENDDO




        POSCH4=POSCH4+2*L3P+1
        ENDDO !ch4
        POSCH3=POSCH3+2*L2P+1
        ENDDO !ch3
        POSCH2=POSCH2+2*L1P+1
        ENDDO !ch2
        POSCH1=POSCH1+2*L0P+1
        ENDDO !ch1
     

     ENDDO



!***********************************************************************
!First add pseudo part contribution on real space grid   < ij | delta(r1-r2) | ij >
!***********************************************************************
      IF (IO%IU0 >= 0) WRITE (*, *) 'Computing pseudo part of Delta_ijij'
      D2PAW_OOOO= zero


      WHF = W
      WHF%WDES => WDES_FOCK
      CALL SETWDES(WHF%WDES, WDESKI, 0)
      CALL SETWDES(WHF%WDES, WDESKJ, 0)


      ALLOCATE (WI(VBMAX) )
      ALLOCATE (WJ(VBMAX) )
      DO NBI = 1, VBMAX
         CALL NEWWAV(WI(NBI), WDESKI, .TRUE.)
         CALL NEWWAV(WJ(NBI), WDESKJ, .TRUE.)
      END DO


      DO KI = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE

         KI_IN_FULL_ORIG = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:, KI), KPOINTS_FULL_ORIG)
         CALL SETWDES(WHF%WDES, WDESKI, KI)
         CALL W1_GATHER_GLB(WHF, 1, VBMAX, ISP, WI)
               ! k_b = k_i - k_q - G
!               KB = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KQ) + WDES%VKPT(:, KI), KPOINTS_FULL)
!               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) Write (*, *) 'error in calc_ftod with shifted k-mesh'
               ! k_a = k_i + k_q - G
!               KA = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - WDES%VKPT(:, KQ), KPOINTS_FULL)
!               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA) == 0)) Write (*, *) 'error in calc_ftod with shifted k-mesh'
      WRITE(*,*)'size of wi',WI(1)%WDES1%GRID%MPLWV
      DO KJ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE

         CALL SETWDES(WHF%WDES, WDESKJ, KJ)
         CALL W1_GATHER_GLB(WHF, 1, VBMAX, ISP, WJ)


         DO NBI = 1, VBMAX
         DO NBJ = 1, VBMAX

            INTC=zero
            DO I=1,WI(NBI)%WDES1%GRID%MPLWV
               INTC=INTC+(CONJG(WI(NBI)%CR(I)*WJ(NBJ)%CR(I))*WI(NBI)%CR(I)*WJ(NBJ)%CR(I))
            ENDDO
            D2PAW_OOOO(NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)), RKIofKI(KJ))=INTC/WI(NBI)%WDES1%GRID%MPLWV

!            WRITE(*,*)'D2PAW_OOOO(',NBI,NBJ,')=',D2PAW_OOOO(NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)),RKIofKI(KJ))
         ENDDO
         ENDDO

     ENDDO !kj loop
     ENDDO !ki loop
 
!***********************************************************************
!compute delta-kernel-contributions of PAW terms from all atoms to D2PAW_OOOO for each state
!***********************************************************************
      IF (IO%IU0 >= 0) WRITE (*, *) 'Computing one-center part of Delta_ijij'

      DO KI = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE

         KI_IN_FULL_ORIG = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:, KI), KPOINTS_FULL_ORIG)
         CALL SETWDES(WHF%WDES, WDESKI, KI)
         CALL W1_GATHER_GLB(WHF, 1, VBMAX, ISP, WI)
      DO KJ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE

         CALL SETWDES(WHF%WDES, WDESKJ, KJ)
         CALL W1_GATHER_GLB(WHF, 1, VBMAX, ISP, WJ)


         DO NBI = 1, VBMAX
         DO NBJ = 1, VBMAX

            ! loop over atoms
            ! select corresponding species and loop over channels
            DO NI=1,W%WDES%NIONS
               SPECIES=W%WDES%ITYP(NI)
               LMMAX=W%WDES%LMMAX(SPECIES)
               IF (LMMAX==0) CYCLE
               NPRO =W%WDES%LMBASE(NI)
!               NPRO2_=WDES1%LMBASE(NI)+NPRO_


               !P(I) or PP ?
               DO CH1=1,LMMAX
!                  WRITE(*,*)'ch1,nbi,npro,cproj',CH1,nbi,npro,WI(NBI)%CPROJ(CH1+NPRO)
               DO CH2=1,LMMAX
               DO CH3=1,LMMAX
               DO CH4=1,LMMAX

!                  D2PAW_OOOO(NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)), RKIofKI(KJ))= D2PAW_OOOO(NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)), RKIofKI(KJ)) + &
!                          (CONJG(W%CPROJ(CH1+NPRO,NBI,1,1)*W%CPROJ(CH2+NPRO,NBJ,1,1))*W%CPROJ(CH3+NPRO,NBI,1,1)*W%CPROJ(CH4+NPRO,NBJ,1,1)) * &
!                          DELTA_PAW(CH1,CH2,CH3,CH4,SPECIES)

                  D2PAW_OOOO(NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)), RKIofKI(KJ))= D2PAW_OOOO(NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)), RKIofKI(KJ)) + &
                          (CONJG(WI(NBI)%CPROJ(CH1+NPRO)*WJ(NBJ)%CPROJ(CH2+NPRO))*WI(NBI)%CPROJ(CH3+NPRO)*WJ(NBJ)%CPROJ(CH4+NPRO)) * &
                          DELTA_PAW(CH1,CH2,CH3,CH4,SPECIES)



               ENDDO !ch4
               ENDDO !ch3
               ENDDO !ch2
               ENDDO !ch1
 
            ENDDO !NI

!            WRITE(*,*)'D2PAW_OOOO(',NBI,NBJ,')=',D2PAW_OOOO(NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)),RKIofKI(KJ))
         ENDDO !nbi
         ENDDO !nbj

     ENDDO !kj loop
     ENDDO !ki loop
 


!***********************************************************************
!Add pseudo part contribution on real space grid   < ij | delta(r1-r2) | ab >
!***********************************************************************
      IF (IO%IU0 >= 0) WRITE (*, *) 'Computing pseudo part of Delta_ijab'
      D2PAW_VVOO= zero

      NSTRIP=30

      WHF = W
      WHF%WDES => WDES_FOCK
      CALL SETWDES(WHF%WDES, WDESKI, 0)
      CALL SETWDES(WHF%WDES, WDESKJ, 0)
      CALL SETWDES(WHF%WDES, WDESKA, 0)
      CALL SETWDES(WHF%WDES, WDESKB, 0)


      IF (ALLOCATED(WI)) DEALLOCATE(WI)
      IF (ALLOCATED(WJ)) DEALLOCATE(WJ)
      IF (ALLOCATED(WA)) DEALLOCATE(WA)
      IF (ALLOCATED(WB)) DEALLOCATE(WB)

      ALLOCATE (WI(VBMAX) )
      ALLOCATE (WJ(VBMAX) )
      ALLOCATE (WA(NSTRIP) )
      ALLOCATE (WB(NSTRIP) )

      DO NBI = 1, VBMAX
         CALL NEWWAV(WI(NBI), WDESKI, .TRUE.)
         CALL NEWWAV(WJ(NBI), WDESKJ, .TRUE.)
      END DO
      DO NBI = 1, NSTRIP
         CALL NEWWAV(WA(NBI), WDESKA, .TRUE.)
         CALL NEWWAV(WB(NBI), WDESKB, .TRUE.)
      END DO



      DO KI = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE

         KI_IN_FULL_ORIG = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:, KI), KPOINTS_FULL_ORIG)
         CALL SETWDES(WHF%WDES, WDESKI, KI)
         CALL W1_GATHER_GLB(WHF, 1, VBMAX, ISP, WI)



      DO KJ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE

         CALL SETWDES(WHF%WDES, WDESKJ, KJ)
         CALL W1_GATHER_GLB(WHF, 1, VBMAX, ISP, WJ)

      DO KA = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA) == 0)) CYCLE
         !CALL SETWDES(WHF%WDES, WDESKA, KPTS_MKPTS(KA))
         CALL SETWDES(WHF%WDES, WDESKA, KA)
         !KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
         KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                     WDES%VKPT(:, KA), KPOINTS_FULL)

!               KA = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - WDES%VKPT(:, KQ), KPOINTS_FULL)
               ! k_b = k_i - k_q - G
         KB = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KQ) + WDES%VKPT(:, KJ), KPOINTS_FULL)

         CALL SETWDES(WHF%WDES, WDESKB, KB)
!               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) Write (*, *) 'error in calc_ftod with shifted k-mesh'
               ! k_a = k_i + k_q - G
!               KA = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - WDES%VKPT(:, KQ), KPOINTS_FULL)
!               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA) == 0)) Write (*, *) 'error in calc_ftod with shifted k-mesh'
!      WRITE(*,*)'size of wi',WI(1)%WDES1%GRID%MPLWV


      IF (ALLOCATED(CPHASE)) DEALLOCATE(CPHASE)
      ALLOCATE(CPHASE(WI(1)%WDES1%GRID%MPLWV))
      CPHASE=(1.0_q,0.0_q)
      CALL SETPHASE(WDES%VKPT(:,KJ)+WDES%VKPT(:,KQ)-WDES%VKPT(:,KB),WI(1)%WDES1%GRID,CPHASE,LPHASE)


      DO NBA=1,NUNOCC,NSTRIP
         NSTRIPA = MIN(NUNOCC - NBA +1, NSTRIP)
         CALL SETWDES(WHF%WDES, WDESKA, KA)
         CALL W1_GATHER_GLB(WHF, VBMAX+NBA, VBMAX+NBA+NSTRIPA-1, ISP, WA)


      DO NBB=1,NUNOCC,NSTRIP
         NSTRIPB = MIN(NUNOCC - NBB +1, NSTRIP)
         CALL SETWDES(WHF%WDES, WDESKB, KB)
         CALL W1_GATHER_GLB(WHF, VBMAX+NBB, VBMAX+NBB+NSTRIPB-1, ISP, WB)


      DO NBAP = 1, NSTRIPA
!         WRITE(*,*)'na',NBA+NBAP-1
      DO NBBP = 1, NSTRIPB

!NBAA + NBA - 1
         DO NBI = 1, VBMAX
         DO NBJ = 1, VBMAX

            INTC=zero
            DO I=1,WI(NBI)%WDES1%GRID%MPLWV
               INTC=INTC+(CONJG(WI(NBI)%CR(I)*WJ(NBJ)%CR(I))*WA(NBAP)%CR(I)*WB(NBBP)%CR(I))*CPHASE(I)
            ENDDO
            D2PAW_VVOO(NBAP + NBA - 1,NBBP + NBB - 1,NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)), RKIofKI(KJ),KA)=INTC/WI(NBI)%WDES1%GRID%MPLWV

!            WRITE(*,*)'D2PAW_OOOO(',NBI,NBJ,')=',D2PAW_OOOO(NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)),RKIofKI(KJ)),WI(NBI)%WDES1%GRID%MPLWV,WA(NBAP)%WDES1%GRID%MPLWV
         ENDDO
         ENDDO

     ENDDO !NBBP
     ENDDO !NBAP

     ENDDO !NBB
     ENDDO !NBA

     ENDDO !ka loop
     ENDDO !kj loop
     ENDDO !ki loop
 

!***********************************************************************
!compute delta-kernel-contributions of PAW terms from all atoms to D2PAW_VVOO for each state
!***********************************************************************


      IF (IO%IU0 >= 0) WRITE (*, *) 'Computing one-center part of Delta_abij'

      DO KI = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI) == 0)) CYCLE

         KI_IN_FULL_ORIG = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:, KI), KPOINTS_FULL_ORIG)
         CALL SETWDES(WHF%WDES, WDESKI, KI)
         CALL W1_GATHER_GLB(WHF, 1, VBMAX, ISP, WI)
      DO KJ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KJ) == 0)) CYCLE

         CALL SETWDES(WHF%WDES, WDESKJ, KJ)
         CALL W1_GATHER_GLB(WHF, 1, VBMAX, ISP, WJ)

      DO KA = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA) == 0)) CYCLE
         !CALL SETWDES(WHF%WDES, WDESKA, KPTS_MKPTS(KA))
         CALL SETWDES(WHF%WDES, WDESKA, KA)
         !KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
         KQ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - &
                                     WDES%VKPT(:, KA), KPOINTS_FULL)

!               KA = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - WDES%VKPT(:, KQ), KPOINTS_FULL)
               ! k_b = k_i - k_q - G
         KB = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KQ) + WDES%VKPT(:, KJ), KPOINTS_FULL)

         CALL SETWDES(WHF%WDES, WDESKB, KB)
!               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) Write (*, *) 'error in calc_ftod with shifted k-mesh'
               ! k_a = k_i + k_q - G
!               KA = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KI) - WDES%VKPT(:, KQ), KPOINTS_FULL)
!               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA) == 0)) Write (*, *) 'error in calc_ftod with shifted k-mesh'
!      WRITE(*,*)'size of wi',WI(1)%WDES1%GRID%MPLWV


      DO NBI = 1, VBMAX
      DO NBJ = 1, VBMAX

         D2PAW_INT= zero
            ! loop over atoms
            ! select corresponding species and loop over channels
         DO NI=1,W%WDES%NIONS
            SPECIES=W%WDES%ITYP(NI)
            LMMAX=W%WDES%LMMAX(SPECIES)
            IF (LMMAX==0) CYCLE
            NPRO =W%WDES%LMBASE(NI)
!               NPRO2_=WDES1%LMBASE(NI)+NPRO_


               !P(I) or PP ?
            DO CH1=1,LMMAX
!                  WRITE(*,*)'ch1,nbi,npro,cproj',CH1,nbi,npro,WI(NBI)%CPROJ(CH1+NPRO)
            DO CH2=1,LMMAX
            DO CH3=1,LMMAX
            DO CH4=1,LMMAX

!!                  D2PAW_OOOO(NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)), RKIofKI(KJ))= D2PAW_OOOO(NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)), RKIofKI(KJ)) + &
!!                          (CONJG(W%CPROJ(CH1+NPRO,NBI,1,1)*W%CPROJ(CH2+NPRO,NBJ,1,1))*W%CPROJ(CH3+NPRO,NBI,1,1)*W%CPROJ(CH4+NPRO,NBJ,1,1)) * &
!!                          DELTA_PAW(CH1,CH2,CH3,CH4,SPECIES)

!                  D2PAW_VVOO(NBA+NBAP-1,NBB+NBBP-1,NBI,NBJ,RKIofKI(KI), RKIofKI(KJ),RKIofKI(KA))= &
!                          D2PAW_VVOO(NBA+NBAP-1,NBB+NBBP-1,NBI,NBJ,RKIofKI(KI),RKIofKI(KJ),RKIofKI(KA)) + &
!                          (CONJG(WI(NBI)%CPROJ(CH1+NPRO)*WJ(NBJ)%CPROJ(CH2+NPRO))*WA(NBA+NBAP-1)%CPROJ(CH3+NPRO)*WB(NBB+NBBP-1)%CPROJ(CH4+NPRO)) * &
!                          DELTA_PAW(CH1,CH2,CH3,CH4,SPECIES)

               D2PAW_INT(CH3,CH4,NI)= &
                       D2PAW_INT(CH3,CH4,NI) + &
                       CONJG(WI(NBI)%CPROJ(CH1+NPRO)*WJ(NBJ)%CPROJ(CH2+NPRO))* &
                       DELTA_PAW(CH1,CH2,CH3,CH4,SPECIES)



            ENDDO !ch4
            ENDDO !ch3
            ENDDO !ch2
            ENDDO !ch1
 
         ENDDO !NI


      DO NBA=1,NUNOCC,NSTRIP
         NSTRIPA = MIN(NUNOCC - NBA +1, NSTRIP)
         CALL SETWDES(WHF%WDES, WDESKA, KA)
         CALL W1_GATHER_GLB(WHF, VBMAX+NBA, VBMAX+NBA+NSTRIPA-1, ISP, WA)


      DO NBB=1,NUNOCC,NSTRIP
         NSTRIPB = MIN(NUNOCC - NBB +1, NSTRIP)

         CALL SETWDES(WHF%WDES, WDESKB, KB)
         CALL W1_GATHER_GLB(WHF, VBMAX+NBB, VBMAX+NBB+NSTRIPB-1, ISP, WB)


      DO NBAP = 1, NSTRIPA
!         WRITE(*,*)'na',NBA+NBAP-1
      DO NBBP = 1, NSTRIPB
!         WRITE(*,*)'nb',NBB+NBBP-1

         DO NI=1,W%WDES%NIONS
            SPECIES=W%WDES%ITYP(NI)
            LMMAX=W%WDES%LMMAX(SPECIES)
!            WRITE(*,*)'ni',NI,LMMAX
            IF (LMMAX==0) CYCLE
            NPRO =W%WDES%LMBASE(NI)
!            WRITE(*,*)'nPRO',NPRO

            DO CH3=1,LMMAX
            DO CH4=1,LMMAX

!!               WRITE(*,*)'',D2PAW_VVOO(NBA+NBAP-1,NBB+NBBP-1,NBI,NBJ,RKIofKI(KI),RKIofKI(KJ),RKIofKI(KA))
!               D2PAW_VVOO(NBA+NBAP-1,NBB+NBBP-1,NBI,NBJ,RKIofKI(KI), RKIofKI(KJ),RKIofKI(KA))= &
!                     D2PAW_VVOO(NBA+NBAP-1,NBB+NBBP-1,NBI,NBJ,RKIofKI(KI),RKIofKI(KJ),RKIofKI(KA)) + &
!                     WA(NBAP)%CPROJ(CH3+NPRO)*WB(NBBP)%CPROJ(CH4+NPRO) * &
!                     D2PAW_INT(CH3,CH4,NI)

            ENDDO !ch4
            ENDDO !ch3
 
         ENDDO !NI


      ENDDO !NBBP
      ENDDO !NBAP
 
      ENDDO !NBB
      ENDDO !NBA

!            WRITE(*,*)'D2PAW_OOOO(',NBI,NBJ,')=',D2PAW_OOOO(NBI,NBJ,RKIofKI(KPTS_MKPTS(KI)),RKIofKI(KJ))
      ENDDO !nbi
      ENDDO !nbj


      ENDDO !ka loop

      ENDDO !kj loop
      ENDDO !ki loop
 


 
   END SUBROUTINE







!***********************************************************************
!
! This routine computes the electronic transition structure factor from
! a set of amplitudes obtained using the PPL term only.
!
!***********************************************************************

   SUBROUTINE CALC_SFACTOR_LADDER(WGW, WDES, IO, EXT)
      USE constant
      USE wave
      USE lattice
      USE full_kpoints
      USE base
      IMPLICIT NONE
      TYPE(wavedes) WGW
      TYPE(wavedes) WDES
      TYPE(in_struct) IO
      CHARACTER(LEN=3), OPTIONAL :: EXT
      INTEGER :: NG, KQ, NI, NJ, KI, KJ, NA, NB, KA, KB, IOIU
      COMPLEX(q) :: ETEST, ETESTD

      IOIU = IO%IU0
      IF (IOIU >= 0) WRITE (*, *) 'Calculating T(G)'

      SDFACTOR = zero
      SXFACTOR = zero

      DO KB = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
         CALL BCAST2ALL_FTOD_AI(WDES, KB, 2, PW_AI_TMP, OC_AI_TMP)
         DO KA = 1, MY_NKPTS
         DO KQ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
            KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
            KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, RKIofKI(KB)) + &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)

            IF (.not. LORBREAL) THEN
               DO NG = 1, NGVECTOR
               DO NI = 1, VBMAX
               DO NJ = 1, VBMAX
               DO NA = 1, NUNOCC
               DO NB = 1, NUNOCC

                  SDFACTOR(NG, RKQofKQ(KQ)) = SDFACTOR(NG, RKQofKQ(KQ)) + (T2_N(NA, NB, NI, NJ, RKIofKI(KI), RKIofKI(KJ), KA))*KPOINTS_FULL%WTKPT(1) &
                                              *KPOINTS_FULL%WTKPT(KI)*KPOINTS_FULL%WTKPT(KI)* &
                                        (FTOD_PW_AI_NOPOT(NG, NI, NA, RKQofKQ(KQ), KA, 1))*CONJG(PW_AI_TMP(NG, NJ, NB, RKQofKQ(KQ)))
               END DO
               END DO
               END DO
               END DO
               END DO
            ELSE
               DO NG = 1, NGVECTOR
               DO NI = 1, VBMAX
               DO NJ = 1, VBMAX
               DO NA = 1, NUNOCC
               DO NB = 1, NUNOCC

                  SDFACTOR(NG, KQ) = SDFACTOR(NG, KQ) + (T2_N_R(NA, NB, NI, NJ, KI, KJ, KA))*KPOINTS_FULL%WTKPT(1) &
                                     *KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1)* &
                                     (FTOD_PW_AI_NOPOT(NG, NI, NA, KQ, KA, 1))*CONJG(PW_AI_TMP(NG, NJ, NB, KQ))
               END DO
               END DO
               END DO
               END DO
               END DO
            END IF

         END DO
         END DO
      END DO

      CALLMPI(M_sum_z(WDES%COMM, SDFACTOR, SIZE(SDFACTOR, 1)*SIZE(SDFACTOR, 2)))

      IF (IOIU >= 0) WRITE (*, *) 'Calculating direct correlation energy ... '
      ETEST = zero
      DO KQ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
         DO NG = 1, NGVECTOR
            ETEST = ETEST + SDFACTOR(NG, RKQofKQ(KQ))*POTFAK_FULL(NG, RKQofKQ(KQ))
         END DO
      END DO
      IF (IOIU >= 0) WRITE (*, *) 'Direct correlation energy for ', EXT, ' channel: ', ETEST*2
      IF (IO%IU6 >= 0) WRITE (IO%IU6, *) 'Direct correlation energy for ', EXT, ' channel: ', ETEST*2
      ETESTD = ETEST*2

      IF (IOIU >= 0) WRITE (*, *) 'Calculating X(G)'

      SXFACTOR = zero

      DO KB = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
         CALL BCAST2ALL_FTOD_AI(WDES, KB, 2, PW_AI_TMP, OC_AI_TMP)
         DO KA = 1, MY_NKPTS
         DO KQ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0)) CYCLE
            KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
            KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) + &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)

            IF (.not. LORBREAL) THEN
               DO NG = 1, NGVECTOR
               DO NI = 1, VBMAX
               DO NJ = 1, VBMAX
               DO NA = 1, NUNOCC
               DO NB = 1, NUNOCC

                  SXFACTOR(NG,RKQofKQ(KQ))=SXFACTOR(NG,RKQofKQ(KQ))+(T2_N(NA,NB,NJ,NI,RKIofKI(KJ),RKIofKI(KI),KA))*KPOINTS_FULL%WTKPT(KI) &
                                            *KPOINTS_FULL%WTKPT(KI)*KPOINTS_FULL%WTKPT(KI)* &
                                        (FTOD_PW_AI_NOPOT(NG, NI, NA, RKQofKQ(KQ), KA, 1))*CONJG(PW_AI_TMP(NG, NJ, NB, RKQofKQ(KQ)))
               END DO
               END DO
               END DO
               END DO
               END DO
            ELSE
               DO NG = 1, NGVECTOR
               DO NI = 1, VBMAX
               DO NJ = 1, VBMAX
               DO NA = 1, NUNOCC
               DO NB = 1, NUNOCC

                  SXFACTOR(NG, KQ) = SXFACTOR(NG, KQ) + (T2_N_R(NA, NB, NJ, NI, KJ, KI, KA))*KPOINTS_FULL%WTKPT(1) &
                                     *KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1)* &
                                     (FTOD_PW_AI_NOPOT(NG, NI, NA, KQ, KA, 1))*CONJG(PW_AI_TMP(NG, NJ, NB, KQ))
               END DO
               END DO
               END DO
               END DO
               END DO
            END IF

         END DO
         END DO
      END DO

      CALLMPI(M_sum_z(WDES%COMM, SXFACTOR, SIZE(SXFACTOR, 1)*SIZE(SXFACTOR, 2)))

      IF (IOIU >= 0) WRITE (*, *) 'Calculatung exchange-like correlation energy ... '
      ETEST = zero
      DO KQ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
         DO NG = 1, NGVECTOR
            ETEST = ETEST + SXFACTOR(NG, RKQofKQ(KQ))*POTFAK_FULL(NG, RKQofKQ(KQ))
         END DO
      END DO
      IF (IOIU >= 0) WRITE (*, *) 'Exchange-like correlation energy for ', EXT, ' channel: ', ETEST
      IF (IO%IU6 >= 0) WRITE (IO%IU6, *) 'Exchange-like correlation energy for ', EXT, ' channel: ', ETEST
      IF (IO%IU6 >= 0) WRITE (IO%IU6, *) 'Total correlation energy for ', EXT, ' channel: ', ETESTD - ETEST

      IF (IOIU >= 0) THEN
         WRITE (*, *) 'Writing to CORRofG ...'
         OPEN (unit=7, file="CORRofG_LADDER")
         WRITE (7, '(7G20.12)') 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
         DO KQ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
            DO NG = 1, NGVECTOR
               IF (GVECLEN(NG, RKQofKQ(KQ)) == 0.0_q) CYCLE
             WRITE(7,'(7G20.12)') GVECX(NG,RKQofKQ(KQ)),GVECY(NG,RKQofKQ(KQ)),GVECZ(NG,RKQofKQ(KQ)),GVECLEN(NG,RKQofKQ(KQ)),POTFAK_FULL(NG,RKQofKQ(KQ)),REAL(SDFACTOR(NG,RKQofKQ(KQ)),kind=q),REAL(SXFACTOR(NG,RKQofKQ(KQ)),kind=q)
            END DO
         END DO
         CLOSE (7)
      END IF


   END SUBROUTINE


!***********************************************************************
!
! This is an experimental routine to manipulate the structure factor
!
!***********************************************************************

   SUBROUTINE CALC_INV_APP_QMATRIX(W, WGW, WDES, IOIU, EXT)
      USE constant
      USE wave
      USE lattice
      USE full_kpoints
      IMPLICIT NONE
      TYPE(wavespin) W
      TYPE(wavedes) WGW
      TYPE(wavedes) WDES
      CHARACTER(LEN=3), OPTIONAL :: EXT
      INTEGER :: NG, NGp, KQ, KQp, KQm, NI, NJ, KI, KJ, NA, NB, KA, KB, IOIU
      COMPLEX(q) :: ETEST
!      needed for q-matrix inversion
      COMPLEX(q), dimension(size(QDMATRIX, 1)) :: WORK  ! work array for LAPACK
      INTEGER, dimension(size(QDMATRIX, 1)) :: IPIV   ! pivot indices
      INTEGER :: N, INFO, INDI, INDJ

      IF (IOIU >= 0) WRITE (*, *) 'Calculating Q_d(G,Gp,kq)'

      QDMATRIX(1:(NGVECTOR*REALNKPTS), 1:(NGVECTOR*REALNKPTS)) = (0.0_q, 0.0_q)
      QDMATRIX(1, 1) = (1.0_q, 0.0_q)
      QXMATRIX(1, 1) = (1.0_q, 0.0_q)

      ETEST = zero

      DO KB = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
         CALL BCAST2ALL_FTOD_AI(WDES, KB, 2, PW_AI_TMP, OC_AI_TMP)
         DO KA = 1, MY_NKPTS
         DO KQ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
            KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
            KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, RKIofKI(KB)) + &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)

            KQp = KPOINT_IN_FULL_GRID(-WDES%VKPT(:, KQ), KPOINTS_FULL)

            DO NG = 1, NGVECTOR
            DO NGp = 1, NGVECTOR
            DO NI = 1, VBMAX
            DO NJ = 1, VBMAX
            DO NA = 1, NUNOCC
            DO NB = 1, NUNOCC

               INDI = NG + (KQ - 1)*NGVECTOR
               INDJ = NGp + (KQp - 1)*NGVECTOR

               QDMATRIX(INDI,INDJ)= &
                  QDMATRIX(INDI,INDJ)+KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(KI)*KPOINTS_FULL%WTKPT(KI)* &
                              (FTOD_PW_AI_NOPOT(NG,NI,NA,RKQofKQ(KQ),KA,1))*CONJG(PW_AI_TMP(NG,NJ,NB,RKQofKQ(KQ)))* &
                              CONJG((FTOD_PW_AI_NOPOT(NGp,NI,NA,RKQofKQ(KQ),KA,1))*CONJG(PW_AI_TMP(NGp,NJ,NB,RKQofKQ(KQ)))) !&

            END DO
            END DO
            END DO
            END DO
            END DO
            END DO

         END DO
         END DO
      END DO

      CALLMPI(M_sum_z(WDES%COMM, QDMATRIX, SIZE(QDMATRIX, 1)*SIZE(QDMATRIX, 2)))

      IF (IOIU >= 0) WRITE (*, *) 'Calculating QX(G,Gp,KQ,KQp)'

      ETEST = zero
      QXMATRIX = zero
      QXMATRIX(1, 1) = (1.0_q, 0.0_q)

      DO KB = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB) == 0)) CYCLE
         CALL BCAST2ALL_FTOD_AI(WDES, KB, 2, PW_AI_TMP, OC_AI_TMP)
         DO KA = 1, MY_NKPTS
         DO KQ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0)) CYCLE
            KI = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
            KJ = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KB) + &
                                     WDES%VKPT(:, KQ), KPOINTS_FULL)
            KQp = KPOINT_IN_FULL_GRID(WDES%VKPT(:, KPTS_MKPTS(KA)) - &
                                      WDES%VKPT(:, KJ), KPOINTS_FULL)

            KQm = KPOINT_IN_FULL_GRID(-WDES%VKPT(:, KQp), KPOINTS_FULL)

            DO NG = 1, NGVECTOR
            DO NGp = 1, NGVECTOR
            DO NI = 1, VBMAX
            DO NJ = 1, VBMAX
            DO NA = 1, NUNOCC
            DO NB = 1, NUNOCC

               QXMATRIX(NG+(RKQofKQ(KQ)-1)*NGVECTOR,NGp+(RKQofKQ(KQm)-1)*NGVECTOR)=QXMATRIX(NG+(RKQofKQ(KQ)-1)*NGVECTOR,NGp+(RKQofKQ(KQp)-1)*NGVECTOR)+ &
                                                             KPOINTS_FULL%WTKPT(KI)*KPOINTS_FULL%WTKPT(KI)*KPOINTS_FULL%WTKPT(KI)* &
                                     (FTOD_PW_AI_NOPOT(NG, NI, NA, RKQofKQ(KQ), KA, 1))*CONJG(PW_AI_TMP(NG, NJ, NB, RKQofKQ(KQ)))* &
                             CONJG((FTOD_PW_AI_NOPOT(NGp, NJ, NA, RKQofKQ(KQp), KA, 1))*CONJG(PW_AI_TMP(NGp, NI, NB, RKQofKQ(KQp))))

            END DO
            END DO
            END DO
            END DO
            END DO
            END DO

         END DO
         END DO
      END DO

      CALLMPI(M_sum_z(WDES%COMM, QXMATRIX, SIZE(QXMATRIX, 1)*SIZE(QXMATRIX, 2)))

      IF (IOIU >= 0) WRITE (*, *) 'Calculating singlet and triplet q-matrices'
      QSMATRIX = QDMATRIX + QXMATRIX
      QSMATRIX(1, 1) = (1.0_q, 0.0_q)
      QTMATRIX = QDMATRIX - QXMATRIX
      QTMATRIX(1, 1) = (1.0_q, 0.0_q)

      DO KQ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
         IF (NGVECTOR > WGW%NGVECTOR(KQ)) THEN
            DO NG = WGW%NGVECTOR(KQ) + 1, NGVECTOR
               QTMATRIX(NG + (RKQofKQ(KQ) - 1)*NGVECTOR, NG + (RKQofKQ(KQ) - 1)*NGVECTOR) = (1.0_q, 0.0_q)
            END DO
         END IF
      END DO

      DO KQ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
         IF (NGVECTOR > WGW%NGVECTOR(KQ)) THEN
            DO NG = WGW%NGVECTOR(KQ) + 1, NGVECTOR
               QSMATRIX(NG + (RKQofKQ(KQ) - 1)*NGVECTOR, NG + (RKQofKQ(KQ) - 1)*NGVECTOR) = (1.0_q, 0.0_q)
            END DO
         END IF
      END DO

      IF (IOIU >= 0) WRITE (*, *) 'Inverting singlet and triplet q-matrices'

      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      QSMATRIX_INV = QSMATRIX
      N = size(QSMATRIX, 1)

      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call ZGETRF(N, N, QSMATRIX_INV, N, IPIV, INFO)

      if (INFO /= 0) then
         stop 'Matrix is numerically singular!'
      end if

      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call ZGETRI(N, QSMATRIX_INV, N, IPIV, WORK, N, INFO)

      if (INFO /= 0) then
         stop 'Matrix inversion failed!'
      end if

      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      QTMATRIX_INV = QTMATRIX
      N = size(QTMATRIX, 1)

      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call ZGETRF(N, N, QTMATRIX_INV, N, IPIV, INFO)

      if (INFO /= 0) then
         stop 'Matrix is numerically singular!'
      end if

      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call ZGETRI(N, QTMATRIX_INV, N, IPIV, WORK, N, INFO)

      if (INFO /= 0) then
         stop 'Matrix inversion failed!'
      end if

      DO KQ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
         IF (NGVECTOR > WGW%NGVECTOR(KQ)) THEN
            DO NG = WGW%NGVECTOR(KQ) + 1, NGVECTOR
               QTMATRIX_INV(NG + (RKQofKQ(KQ) - 1)*NGVECTOR, NG + (RKQofKQ(KQ) - 1)*NGVECTOR) = (0.0_q, 0.0_q)
            END DO
         END IF
      END DO

      DO KQ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
         IF (NGVECTOR > WGW%NGVECTOR(KQ)) THEN
            DO NG = WGW%NGVECTOR(KQ) + 1, NGVECTOR
               QSMATRIX_INV(NG + (RKQofKQ(KQ) - 1)*NGVECTOR, NG + (RKQofKQ(KQ) - 1)*NGVECTOR) = (0.0_q, 0.0_q)
            END DO
         END IF
      END DO

      IF (IOIU >= 0) WRITE (*, *) 'Calculating singlet and triplet correlation factor by Q^-1 S'

      SSFACTOR = SDFACTOR + SXFACTOR
      STFACTOR = SDFACTOR - SXFACTOR
      F12SFACTOR(:, :) = zero
      F12TFACTOR(:, :) = zero
      QSMATRIX_INV(1, 1) = (0.0_q, 0.0_q)
      QTMATRIX_INV(1, 1) = (0.0_q, 0.0_q)

      DO KQ = 1, WDES%NKPTS
         IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
         DO KQP = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQP) /= 0.0_q)) CYCLE

            DO NG = 1, NGVECTOR
            DO NGp = 1, NGVECTOR
               F12SFACTOR(NG, RKQofKQ(KQ)) = F12SFACTOR(NG, RKQofKQ(KQ)) + &
                        QSMATRIX_INV(NG + (RKQofKQ(KQ) - 1)*NGVECTOR, NGp + (RKQofKQ(KQp) - 1)*NGVECTOR)*SSFACTOR(NGp, RKQofKQ(KQp))

               F12TFACTOR(NG, RKQofKQ(KQ)) = F12TFACTOR(NG, RKQofKQ(KQ)) + &
                        QTMATRIX_INV(NG + (RKQofKQ(KQ) - 1)*NGVECTOR, NGp + (RKQofKQ(KQp) - 1)*NGVECTOR)*STFACTOR(NGp, RKQofKQ(KQp))
            END DO
            END DO

         END DO
      END DO

      IF (IOIU >= 0) THEN
         WRITE (*, *) 'Writing to F12ofG ...'
         OPEN (unit=7, file='F12ofG.'//EXT)
         DO KQ = 1, WDES%NKPTS
            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ) /= 0.0_q)) CYCLE
            DO NG = 1, NGVECTOR
               WRITE(7,'(7G20.12)') GVECX(NG,RKQofKQ(KQ)),GVECY(NG,RKQofKQ(KQ)),GVECZ(NG,RKQofKQ(KQ)),GVECLEN(NG,RKQofKQ(KQ)),REAL(F12SFACTOR(NG,RKQofKQ(KQ)),kind=q),REAL(F12TFACTOR(NG,RKQofKQ(KQ)),kind=q),REAL(SSFACTOR(NG,RKQofKQ(KQ)),kind=q)
            END DO
         END DO
         CLOSE (7)
      END IF

   END SUBROUTINE

!******************** SUBROUTINE CCINCAR_READER ************************
!
! This is the main INCAR reader for the CCSD module to control the solver,
! various approximations and the in- and output.
!
!***********************************************************************
   SUBROUTINE CCINCAR_READER(IU5, IU6, IU0)
      USE reader_tags
      USE string, ONLY: lowercase
      INTEGER IU5, IU6, IU0
      ! local variables
      INTEGER N, IERR
      REAL(q) RDUM
      LOGICAL LOPEN
      INTEGER                    :: NUM_ELEMENTS
      CHARACTER(:), ALLOCATABLE  :: READ_STRING
      CHARACTER(:), ALLOCATABLE  :: CCALGOSTRING
      INTEGER, EXTERNAL :: LENGTH

      CCMP2 = .FALSE.
      LOPEN = .TRUE.
      READ_STRING = ''
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'CCALGO', READ_STRING, 0, IERR, WRITEXMLINCAR, FOUNDNUMBER=NUM_ELEMENTS)
      READ_STRING = REPEAT(' ', NUM_ELEMENTS)
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'CCALGO', READ_STRING, 0, IERR, WRITEXMLINCAR)
      CCALGOSTRING = TRIM(ADJUSTL(lowercase(READ_STRING)))
      N = NUM_ELEMENTS
      IF (CCALGOSTRING(1:3) == 'mp2') CCMP2 = .TRUE.
      IF (CCMP2) SINGLES = .FALSE.
      LCCD = .FALSE.
      IF (CCALGOSTRING(1:4) == 'lccd') LCCD = .TRUE.
      RING = .FALSE.
      IF (CCALGOSTRING(1:4) == 'ring') RING = .TRUE.
      LDISTING = .FALSE.
      IF (CCALGOSTRING(1:7) == 'disting') LDISTING = .TRUE.

      MAX_IT = 10
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'CCMAXIT', MAX_IT, IERR, WRITEXMLINCAR)

      NELM = 0
      LSFACTOR = .FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSFACTOR', LSFACTOR, IERR, WRITEXMLINCAR)

      LWRITET = .FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LWRITET', LWRITET, IERR, WRITEXMLINCAR)

      LPSPPL = .FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LPSPPL', LPSPPL, IERR, WRITEXMLINCAR)

      CCMIX = 0.5_q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'CCMIX', CCMIX, IERR, WRITEXMLINCAR)

      CCNINT = 25
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'CCNINT', CCNINT, IERR, WRITEXMLINCAR)

      RETURN

   END SUBROUTINE CCINCAR_READER

!***********************************************************************
!
! This routine determines the type of k-mesh employed.
!
!***********************************************************************

   SUBROUTINE CHECK_SHIFTED_KPOINTS(WDES, W, KPOINTS)
      USE constant
      USE full_kpoints
      USE mkpoints
      USE wave
      implicit NONE
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      TYPE(kpoints_struct) KPOINTS
      integer :: MKI, KI
      LOGICAL :: GAMMA_FOUND
      LOGICAL :: SECOND_GAMMA
      REAL(q) :: SECOND_GAMMA_WEIGHT
      REAL(q) :: GAMMA_WEIGHT
      INTEGER :: NRKI, NRKQ
      REAL(q) :: TINY = 1.D-8

      GAMMA_FOUND = .FALSE.
      SHIFTED_KPOINTS = .TRUE.
      SECOND_GAMMA = .FALSE.
      do MKI = 1, WDES%NKPTS
   if ((ABS(W%WDES%VKPT(1, MKI)) < (TINY)) .and. (ABS(W%WDES%VKPT(2, MKI)) < (TINY)) .and. (ABS(W%WDES%VKPT(3, MKI)) < (TINY))) then
            IF (GAMMA_FOUND) SECOND_GAMMA = .TRUE.
            IF (SECOND_GAMMA) SECOND_GAMMA_WEIGHT = W%WDES%WTKPT(MKI)
            IF (.NOT. GAMMA_FOUND) THEN
               GAMMA_WEIGHT = W%WDES%WTKPT(MKI)
               GAMMA_FOUND = .true.
            END IF
         end if
      end do
      if ((GAMMA_FOUND) .and. (GAMMA_WEIGHT /= 0.0_q)) THEN
         SHIFTED_KPOINTS = .FALSE.
      end if
      IF (SECOND_GAMMA) SHIFTED_KPOINTS = .TRUE.

      REALNKPTS = WDES%NKPTS
      IF (SHIFTED_KPOINTS) REALNKPTS = WDES%NKPTS/2

      ALLOCATE (RKQofKQ(WDES%NKPTS))
      ALLOCATE (RKIofKI(WDES%NKPTS))
      ALLOCATE (KQofRKQ(REALNKPTS))
      ALLOCATE (KIofRKI(REALNKPTS))
      NRKI = 0
      NRKQ = 0
      RKQofKQ = 0
      RKIofKI = 0
      KQofRKQ = 0
      KIofRKI = 0
      IF (SHIFTED_KPOINTS) THEN
         DO KI = 1, WDES%NKPTS
            IF ((WDES%WTKPT(KI) /= 0.00_q)) THEN
               NRKI = NRKI + 1
               RKIofKI(KI) = NRKI
               KIofRKI(NRKI) = KI
            END IF
            IF ((WDES%WTKPT(KI) == 0.00_q)) THEN
               NRKQ = NRKQ + 1
               RKQofKQ(KI) = NRKQ
               KQofRKQ(NRKQ) = KI
            END IF

         END DO
      ELSE
         DO KI = 1, WDES%NKPTS
            RKIofKI(KI) = KI
            KIofRKI(KI) = KI
            RKQofKQ(KI) = KI
            KQofRKQ(KI) = KI
         END DO
      END IF

      WTKPT = 1.0_q/REALNKPTS

   END SUBROUTINE CHECK_SHIFTED_KPOINTS

#endif // gammareal
#endif // scaLAPACK

END MODULE ccsd
EOF
    echo "Adding file  ccsd.F " 
cp .tmp_vaspcc4s_patch  ccsd.F
