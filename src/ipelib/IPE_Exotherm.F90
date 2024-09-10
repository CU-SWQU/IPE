module IPE_Exotherm
    use spher_harm
    use iso_fortran_env
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,307)
    INTEGER, PARAMETER :: prec = dp
    integer, parameter :: nspecies= 10
    REAL(prec), PARAMETER :: pi      = 3.14159265358979323844_prec
    REAL(prec), PARAMETER :: rtd     = 180.0_prec / pi   !radian-->degree
    REAL(prec), PARAMETER :: dtr     = pi / 180.0_prec   !degree-->radian
  
    contains 

subroutine exotherm_IPE_output(iyd, msis_sec, msis_alt, msis_lat, msis_lon, dens, temps)
    implicit none
    !
          real(dp) , intent(in) :: msis_alt, msis_lat , msis_lon , msis_sec
          integer, intent(in) :: iyd

          real(dp) , intent(inout) , dimension(2) :: temps
          real(dp) , intent(inout) , dimension(9) :: dens
          character(:), allocatable :: monaco_file, tgcm_file, coeffsdir_start, monaco_coeffs_dir, tgcm_coeffs_dir
          
          real(dp), allocatable ::  X(:,:,:),sh_height(:),X_tgcm(:,:,:),avg_zgmid(:)
          real(8), allocatable :: t_pmn(:), m_pmn(:)
          real(dp) :: monaco_height, hden, heden, n2den, o2den, oden, tex, tn, wgt, h1_m3, he_m3, o1_m3
          integer :: tgcm_len_z, monaco_len_z, k, m_nmax, t_nmax

        
    ! inputs:
    ! iyd,         &    ! Input, year and day as YYDDD (day of year from 1 to 365 (or 366))
    ! msis_sec,    &    ! Input, universal time ( sec )
    ! msis_alt,    &    ! Input, altitude ( km )
    ! msis_lat,    &    ! Input, geodetic latitude ( degrees )
    ! msis_lon,    &    ! Input, geodetic longitude ( degrees )
    ! outputs:
     !     D(1) - HE NUMBER DENSITY(CM-3)
     !     D(2) - O NUMBER DENSITY(CM-3)
     !     D(3) - N2 NUMBER DENSITY(CM-3)
     !     D(4) - O2 NUMBER DENSITY(CM-3)
     !     not needed: D(5) - AR NUMBER DENSITY(CM-3)
     !     not needed: D(6) - TOTAL MASS DENSITY(GM/CM3)
     !     D(7) - H NUMBER DENSITY(CM-3)
     !     D(8) - N NUMBER DENSITY(CM-3)
     !     not needed: D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
     !     T(1) - EXOSPHERIC TEMPERATURE
     !     T(2) - TEMPERATURE AT ALT
    
    
    print *, msis_sec, msis_alt, msis_lat, msis_lon
    coeffsdir_start = '/Users/sarahluettgen/Documents_local/ipe_evaluation/sh_coeffs_';
    
    monaco_coeffs_dir = coeffsdir_start//'monaco'
    tgcm_coeffs_dir= coeffsdir_start//'tgcm'

    call get_SH_coeff_file(msis_sec/3600., monaco_coeffs_dir, monaco_file)
    call get_SH_coeff_file(msis_sec/3600., tgcm_coeffs_dir, tgcm_file)
 
    call get_model_coeffs(monaco_file,m_nmax,monaco_len_z,X,     sh_height) 
    call get_model_coeffs(tgcm_file,  t_nmax,tgcm_len_z,  X_tgcm,avg_zgmid) 

    print *, X_tgcm(3,1,:)
    
    sh_height = sh_height/1e3; ! km
    monaco_height = sh_height(2); ! km
    
    if(minval(sh_height) > maxval(avg_zgmid)) then
         print *, 'gap between max(zgmid) and min(monaco shell height)!'
         stop
    end if
   
    ! we will never use these values
    dens(5) = 0;
    dens(6) = 0;
    dens(9) = 0;

    allocate(t_pmn((t_nmax+1)*(t_nmax+1)))
    allocate(m_pmn((m_nmax+1)*(m_nmax+1)))
    call SPHAR(t_nmax, msis_lon*dtr,msis_lon*dtr,t_pmn) ! TO DO

    call SPHAR(m_nmax, msis_lon*dtr,msis_lon*dtr,m_pmn) ! TO DO
    
    if(msis_alt < maxval(avg_zgmid)) then
        ! if the height we want is within TIMEGCM, get N2, O2, N from there
        ! k=find(avg_zgmid>msis_alt,1)-1; !k is shell below msis_alt, k+1 is above 
        print *, 'N2, O2 from TIME-GCM'
        do k=1,size(avg_zgmid)-1
            if(avg_zgmid(k+1)>msis_alt) EXIT
        end do ! k=lev0,lev1-1  
        wgt = (msis_alt-avg_zgmid(k))/(avg_zgmid(k+1)-avg_zgmid(k))

     
        call SH_synthesis(t_pmn,X_tgcm(3,:,:),t_nmax,k,wgt,size(avg_zgmid),n2den) ! SH Synthesis
        call SH_synthesis(t_pmn,X_tgcm(4,:,:),t_nmax,k,wgt,size(avg_zgmid),o2den) ! SH Synthesis
        dens(3) = n2den
        dens(4) = o2den
        ! dens(8,:) = (1.-wgt)*SH_fit_synthesis(msis_lon,msis_lat,X_tgcm{8,k},nmax) + (wgt)*SH_fit_synthesis(msis_lon,msis_lat,X_tgcm{8,k+1},nmax); % SH Synthesis
       
    
        if(msis_alt <= monaco_height) then
            ! if we cannot pull Tn, He, O, H from Monaco, get them from TIMEGCM too
            print *,  'Tn, He, O, H from TIME-GCM'
            call SH_synthesis(t_pmn,X_tgcm(1,:,:),t_nmax,k,wgt,size(avg_zgmid), heden) ! SH Synthesis
            call SH_synthesis(t_pmn,X_tgcm(2,:,:),t_nmax,k,wgt,size(avg_zgmid), oden) ! SH Synthesis
            call SH_synthesis(t_pmn,X_tgcm(7,:,:),t_nmax,k,wgt,size(avg_zgmid), hden) ! SH Synthesis
            
            dens(1) = heden
            dens(2) = oden
            dens(7) = hden

            call SH_synthesis(t_pmn,X_tgcm(10,:,:),t_nmax,size(avg_zgmid)-1,wgt,size(avg_zgmid),tex) ! SH Synthesis
            call SH_synthesis(t_pmn,X_tgcm(10,:,:),t_nmax,k,wgt,size(avg_zgmid),tn) ! SH Synthesis
            
            temps(1) = tex
            temps(2) = tn


        end if ! msis_alt <= monaco_height
    
    else ! out of range of TIME-GCM
        dens(3) = 0 ! N2
        dens(4) = 0 ! O2
        dens(8) = 0 ! N
    end if
    
    if(msis_alt > monaco_height) then ! if we CAN pull Tn, He, O, H from Monaco...
        print *, 'Tn, He, O, H from Monaco'
        ! Calculate weights for height interpolation
        do k=1,size(sh_height)-1
            if(sh_height(k+1)>msis_alt) EXIT
        end do ! k=lev0,lev1-1  
        wgt = (msis_alt-sh_height(k)) / (sh_height(k+1)-sh_height(k))
    
        call SH_synthesis(m_pmn,X(1,:,:),m_nmax,k,wgt,size(sh_height),he_m3) ! SH Synthesis
        call SH_synthesis(m_pmn,X(2,:,:),m_nmax,k,wgt,size(sh_height),o1_m3) ! SH Synthesis
        call SH_synthesis(m_pmn,X(7,:,:),m_nmax,k,wgt,size(sh_height),h1_m3) ! SH Synthesis
        
        dens(1) = he_m3/(100.**3);
        dens(2) = o1_m3/(100.**3);
        dens(7) = h1_m3/(100.**3);
    
        call SH_synthesis(m_pmn,X(10,:,:),m_nmax,k,wgt,size(sh_height),tn) ! SH Synthesis
        call SH_synthesis(m_pmn,X(10,:,:),m_nmax,2,wgt,size(sh_height),tex) ! assume second shell is ~ exobase
        
        temps(1) = tex;
        temps(2) = tn;
    
    end if
    print *, 'den = ', dens
    print *, 'tn = ', temps
end subroutine exotherm_IPE_output
    
    
    subroutine SH_synthesis(pmn,sh_coeffs,nmax,k,wgt,zmax,val)
       
        integer,  intent(in) :: nmax,k,zmax
        real(dp) , intent(in) :: wgt
        real(dp) , intent(in) , dimension((nmax+1)*(nmax+1)) :: pmn
        real(dp) , intent(in) , dimension(zmax,(nmax+1)*(nmax+1)) :: sh_coeffs
        real(dp), intent(out) :: val

        integer :: ifunc,m,n

        val = 0
        ! print *, 'sh_coeffs size = ',size(sh_coeffs,1),size(sh_coeffs,2)
        do ifunc=1,size(pmn)
         ! print *, 'P, SH = ',pmn(ifunc), sh_coeffs(k,ifunc)
          val = val+pmn(ifunc)*((1.-wgt)*sh_coeffs(k,ifunc)+wgt*sh_coeffs(k+1,ifunc))
        enddo 
 

    end subroutine SH_synthesis



    subroutine get_model_coeffs(filenam,nmax,len_zarr,sh_coeffs,altitudes)
        implicit none
        
        character(:), allocatable :: filenam
        integer :: nmax,len_zarr,nspecs,specnum, i , j, coeflen
        real(dp), allocatable :: altitudes(:)
        real(dp), allocatable :: sh_coeffs(:,:,:)
        real(dp), allocatable :: sh_coeffs_at_alt(:)

        
        open(unit=15,file=filenam, status='old',access='sequential', form='formatted', action='read')
        read (15,*)  nmax,len_zarr, nspecs
        ! print *, filenam, ' with nmax = ', nmax
        coeflen = (nmax+1)*(nmax+1)
        allocate(sh_coeffs_at_alt(coeflen))
        allocate(sh_coeffs(nspecies,len_zarr,coeflen))
        allocate(altitudes(len_zarr))
        read (15,*) altitudes
        do i = 1,nspecs
            read (15,*)  specnum
            ! print *, specnum
            do j = 1,len_zarr
                read (15, *)  sh_coeffs_at_alt
                
                sh_coeffs(specnum,j,:) = sh_coeffs_at_alt;
                !print *, sh_coeffs_at_alt
            end do
            
        end do
       ! print *, size(sh_coeffs,3)
        ! print *, sh_coeffs(1,1,:)
        
        close(unit=15)
    
    
    
    print *, 'get_model_coefs done!'
    end subroutine get_model_coeffs


    subroutine get_SH_coeff_file(ut_want, coeff_dir, final_best_time_file)
        ! use iso_fortran_env
        implicit none
        character(:), allocatable, intent(in) :: coeff_dir
        real(dp), intent(in) :: ut_want
        character(:), allocatable, intent(out) :: final_best_time_file

        character(len=*), parameter :: ls_file = 'coeff_file.txt'
        integer :: u, ios,last_underscore_loc
        character(100) :: filename, str_time,best_time_file
        real(dp) :: real_time, best_time_diff
        logical :: is_first

        print *, 'get_SH_coeff_file called!'
        call execute_command_line('ls -1 '//coeff_dir//' > '//ls_file, wait=.TRUE., exitstat=ios)
        if (ios /= 0) stop 'ls command did not work'

        open(newunit=u, file=ls_file, iostat=ios, status="old", action="read")
        if ( ios /= 0 ) stop 'Could not open list of SH coeff files'
        is_first = .true.
        do
            print *, 'next do'
            read(u, *, iostat=ios) filename
            print *, 'read line: ', filename
            if (is_iostat_end(ios)) exit
            if (ios /= 0) STOP "Unexpected error while reading listing file"

            if (index(filename, ".txt") > 0) then
                print *, filename
                last_underscore_loc = index(filename, '_', .TRUE.)
                str_time = filename( (last_underscore_loc+1) : (len_trim(filename)-4) )
                
                read(str_time, '(F7.5)') real_time ! 2 digits before decimal pt + 5 after 
                ! print*, filename, ' time = ',real_time, ' diff = ', abs(real_time - time_want)
                if (abs(real_time - ut_want) < best_time_diff .or. is_first) then
                    
                    is_first = .FALSE.
                    print *, is_first
                    best_time_diff = abs(real_time - ut_want)
                    print *, 'new best time: ', best_time_diff
                    best_time_file = filename
                    print *, 'new best file: ',best_time_file
                end if 
                print *, 'endif 1'
            end if
            print *, 'endif 2'
        end do
        ! print *, best_time_file
        close(u)
        if ((best_time_diff > 0.5) .or. is_first) then
            stop 'No appropriate files found!'
        else
            final_best_time_file = coeff_dir//'/'//trim(best_time_file)
            print *, 'Chose ', best_time_file
        end if

        ! call execute_command_line('rm '//ls_file, wait=.FALSE.)
        
    end subroutine get_SH_coeff_file

end module IPE_Exotherm
