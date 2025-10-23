c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Notices:
! Copyright 2025 United States Government as represented by the Administrator of the National Aeronautics and Space Administration. 
! No copyright is claimed in the United States under Title 17, U.S. Code. All Other Rights Reserved.

! The NASA Software “Progressive Release Explicit Virtual Crack Closure Technique” (LAR-20580-1) 
! calls the following third-party software, which is subject to the terms and conditions of its licensor, as applicable at the time of licensing.  
! The third-party software is not bundled or included with this software but may be available from the licensor.  
! License hyperlinks are provided here for information purposes only.

! NumPy
! https://numpy.org/doc/stable/license.html
! https://numpy.org/devdocs/license.html

! Copyright (c) 2005-2025, NumPy Developers.
! All rights reserved.

! Disclaimers
! No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, 
! ! BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, 
! ! FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY 
! ! WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN 
! ! ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER 
! ! APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING 
! ! THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS." 

! Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, 
! ! AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES 
! ! ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, 
! ! RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS, AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, 
! ! TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)
c     
      include 'aba_param.inc'
    
      dimension time(2)

c     ***************common block 1************
      double precision jelem_tipels(200000,21) !input
      double precision jelem_force(200000,13)  !output
      double precision jelem_delta(200000,13)  !output
      double precision jelem_onset(200000,14)  !output (el_id, forces, "nodal" area)
      integer ntipels
    
      double precision ncycles_prev,ncycles_current

      common /crackdata/ jelem_tipels,
     &     jelem_force,jelem_delta,jelem_onset,
     &     ntipels,ncycles_prev,ncycles_current
c    *****************************************            

c     ***common block 2***
      double precision pnewdt_global
      common /sol_contr_exp/ pnewdt_global
c     ********************

c     ***common block 3*** 
      double precision jelem_dstat(200000,4)
      integer jelem_prxvcct
      integer jelem_ghost
      logical viz
      common /dstat_common/jelem_dstat,jelem_prxvcct,
     &jelem_ghost,viz     
c     ***end common block***

      CHARACTER (LEN=256) :: DIRECTORY
      integer lenOutputDir
c      CHARACTER (LEN=*), PARAMETER :: DIRECTORY=
      logical file_exists, file_opened
      integer i,j,nnode,mcrd
      integer  max_actelements,max_elements
      integer read_file_aux
c     hardcoded maximum number of elements and activated elements
      max_actelements = 200000
      max_elements = 150000

      mcrd = 3
      nnode_side = 4

      call getoutdir(DIRECTORY, lenOutputDir)
      
      if  (lop .eq. 0) then
         u_lop = lop
c     cut back control, initiatlized to be inactive (1E6 >>> 1)
         pnewdt_global = 1.0E6 
      
c     initializes the cycle count
         ncycles_prev = 0.d0
         ncycles_current = 0.d0
         
c     initializes jelem matrices
         jelem_onset = 0.d0
         jelem_dstat = 0.d0
         jelem_force = 0.d0
         jelem_delta = 0.d0
         jelem_tipels = 0.d0
         
c         call initialize_matrix(jelem_dstat,max_actelements,nnode_side)
         
c         call initialize_matrix(jelem_force,max_actelements,
c     &        mcrd*nnode_side+1) !no, f_node_xyz
c         call initialize_matrix(jelem_delta,max_actelements,
c     &        mcrd*nnode_side+1)
c         call initialize_matrix(jelem_tipels,max_actelements,21)
c         ntipels = 0
          
c     calls the python function that loads the mesh
         call system (trim(directory)//'/'//'update_elements_nvcct.sh') !shell script that runs a local python
 
c     determine if the last file exists
         file_exists = .false. 
         do while (file_exists .eq. .false.)
            inquire(file=trim(directory)//'/'//'TIP_ELS.txt',
     &           exist=file_exists)
            write(7,*) 'file_exists',file_exits
            call sleep(1)
         end do      
         
c     reads tip elements at the interface         
         open(103, file=trim(directory)//'/'//'TIP_ELS.txt')
         read(103,*) ntipels
         do i=1,ntipels
            read(103,*) (jelem_tipels(i,j),j=1,21)
            jelem_force(i,1) = jelem_tipels(i,1)
            jelem_delta(i,1) = jelem_tipels(i,1)
         end do
         close(103)

         call read_viz_txt(directory,'viz.uel',viz,jelem_prxvcct,
     &        jelem_ghost)

      end if

      if (lop .eq. 1) then       

c     reads number of cycles
         inquire(file=trim(directory)//'/'//'ncycles.txt',
     &        exist=file_exists)
         if (file_exists .eq. .True.) then
            open(103, file=trim(directory)//'/'//'ncycles.txt')
            ncycles_prev = ncycles_current
            do 
               read(103,*,IOSTAT=read_file_aux) ncycles_current
               if (read_file_aux < 0) then !eof was reached
                  exit
               end if
            end do
            close(103)
         end if

c     reads cut-back control
         open(103, file=trim(directory)//'/'//'pnewdt.txt')
         read(103,*) pnewdt_global        
         close(103)

c     initialize info exchange arrays
         jelem_onset = 0.d0
         jelem_force = 0.d0
         jelem_delta = 0.d0
         jelem_tipels = 0.d0

         
c         call initialize_matrix(jelem_force,max_actelements,
c     &        mcrd*nnode_side+1)
c         call initialize_matrix(jelem_delta,max_actelements,
c     &        mcrd*nnode_side+1)
c         call initialize_matrix(jelem_tipels,max_actelements,21)
c         ntipels = 0
         
c     determine if the last file printed by the script exists
         file_exists = .false. 
         do while (file_exists .eq. .false.)
            inquire(file=trim(directory)//'/'//'TIP_ELS.txt',
     &           exist=file_exists)
            write(7,*) 'file_exists',file_exits
            call sleep(1)
         end do     
         
c     reads tip elements at the interface         
         open(103, file=trim(directory)//'/'//'TIP_ELS.txt')
         read(103,*) ntipels
         do i=1,ntipels
            read(103,*) (jelem_tipels(i,j),j=1,21)
            jelem_force(i,1) = jelem_tipels(i,1)
            jelem_delta(i,1) = jelem_tipels(i,1)
         end do
         close(103)

      end if
      
      if (lop .eq. 2) then
         u_lop = lop

c     write step info
         open(103,file=trim(directory)//'/'//'STEP.txt',     
     &        status='unknown',iostat=ierr)
         if ( ierr .eq. 0) then
            close(103, status='delete')
         endif
         open(103,file=trim(directory)//'/'//'STEP.txt')
         write(103,*),kstep,kinc
         write(103,*)
         close(103)

         open(103,file=trim(directory)//'/'//'ELONSET.txt',     
     &        status='unknown',iostat=ierr)
         if ( ierr .eq. 0) then
            close(103, status='delete')
         endif
         
         open(103, file=trim(directory)//'/'//'ELONSET.txt')
         do i=1, max_actelements
            if (jelem_onset(i,1) .ne. 0.d0) then
               do j=1,mcrd*nnode_side+2
                  write(103,'(e25.17$)') jelem_onset(i,j)
               end do
               write(103,*)
            end if
         end do
         close(103)
      
c      update forces computed at eltip els
         open(103,file=trim(directory)//'/'//'ELFORCE.txt',     
     &        status='unknown',iostat=ierr)
         if ( ierr .eq. 0) then
            close(103, status='delete')
         endif
         
         open(103, file=trim(directory)//'/'//'ELFORCE.txt')
         do i=1,ntipels
            do j=1,mcrd*nnode_side+1
               write(103,'(e25.17$)') jelem_force(i,j)
            end do
            write(103,*)
         end do
         close(103)

c     update displacements computed at eltip els
         open(103,file=trim(directory)//'/'//'ELDELTA.txt',     
     &        status='unknown',iostat=ierr)
         if ( ierr .eq. 0) then
            close(103, status='delete')
         endif
         
         open(103, file=trim(directory)//'/'//'ELDELTA.txt')
         do i=1,ntipels
            do j=1,mcrd*nnode_side+1
               write(103,'(e25.17$)') jelem_delta(i,j)
            end do
        
            write(103,*)
         end do
         close(103)

         open(103, file=trim(directory)//'/'//'TIP_ELS.txt',
     &        status='unknown',iostat=ierr)
         if ( ierr .eq. 0) then
            close(103, status='delete')
         endif
 
         call system (trim(directory)//'/'//'update_elements_nvcct.sh') !shell script that runs a local python
         
      end if
      
      
      return
      end 


c********************************************************************************************
c********************* subroutine uel for floating node method *******************************
c********************************************************************************************
      subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     &       props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,
     &       kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf,
     &       lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njprop,period)
c
	  include 'aba_param.inc'

      dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*)
      dimension svars(*),energy(8),coords(mcrd,nnode),u(ndofel)
      dimension du(mlvarx,*),v(ndofel),a(ndofel),time(2),params(*)
      dimension jdltyp(mdload,*),adlmag(mdload,*),ddlmag(mdload,*)
      dimension predef(2,npredf,nnode),lflags(*),jprops(*)

c     ***************common block 1************
      double precision jelem_tipels(200000,21) !input
      double precision jelem_force(200000,13)  !output
      double precision jelem_delta(200000,13)  !output
      double precision jelem_onset(200000,14)  !output (el_id, forces, "nodal" area)
      integer ntipels
    
      double precision ncycles_prev,ncycles_current

      common /crackdata/ jelem_tipels,
     &     jelem_force,jelem_delta,jelem_onset,
     &     ntipels,ncycles_prev,ncycles_current
c    *****************************************            
      
c     ***common block 2***
      double precision pnewdt_global
      common /sol_contr_exp/ pnewdt_global
c     ********************          

c     ***common block 3*** 
      double precision jelem_dstat(200000,4)
      integer jelem_prxvcct
      integer jelem_ghost
      logical viz
      common /dstat_common/jelem_dstat,jelem_prxvcct,
     &jelem_ghost,viz     
c     ***end common block***

c     variables associated with the onset methodology:
      integer stress_crit_type  !1 - quadratic stress criterion
      integer nintp !number of integration points the underlying cohesive elements (for calculation of the area)
      double precision X_st(3)  !3 strenghts
      double precision fc_thd !threshold for considering the node may have failed (final check performed by the python file
      logical f_bool            !weather the stress at one of the ones exceeded the threshold
      character(len=15):: element
      
c     remaining variables:
      integer max_actelements,no_st_props
      double precision st_props(4) !strength properties
      double precision kmatrx(ndofel,ndofel), ku(ndofel)
      double precision delta_temp
      
      integer eltip_ind      !index to the current element in the jelem_tractdelta_i matrix
      double precision area_n !avg area associated with each node
      double precision dstat(nnode/2)
      double precision m_xyz(ndofel/2)
      double precision penalty_f(nnode/2)
      double precision nforce(ndofel/2) !forces at the nodes, at the interface of a given element
      double precision delta(ndofel/2) !openings at the nodes, at the interface of a given element
      double precision rot_m(ndofel,ndofel)
      double precision coords_lc(mcrd,ndofel)
      integer nvars_status
      integer nvars_m_xyz,nvars_dstat,nvars_penalty_f
      integer nvars_ncycle_del
      CHARACTER (LEN=256) :: DIRECTORY
      integer lenOutputDir

       
c      CHARACTER (LEN=*), PARAMETER :: DIRECTORY=
      logical  nlgeom, contact
      logical file_exists, file_opened
      
c     nonlinear geometry flag (set via python script; defaul '.False')
      nlgeom=.False.
c     nodal based contact for delaminations flag (set via python script; default '.False')
      contact=.True.
 
c     state variables definition: 
c     svars(1): interface status1
      nvars_status = 1
c     svars(2,...,13): [m_z,m_x,m_y] slopes of unloading curves for interface nodes
      nvars_m_xyz = 1       ![[2,3,4,5],[6,7,8,9],[10,11,12,13]] = [m_z_node1...node_4,m_x_node_1...node_5,m_y_node_1...node_4]
c     svars(14,...,17): dstat_i dstat variables for interface nodes
      nvars_dstat = 13       ![14,15,16,17]
c     svars(18,...,21): factor multipling node integration 1/no_els sharing the node
      nvars_penalty_f = 17   ![18,19,20,21]
c     svars(22): interface, cycles at failure  (visualization)
      nvars_ncycle_del = 22 

      max_actelements = 200000

      element = 'quadrilateral'
      nintp = 4


       
      if (pnewdt_global .lt. 1.d0) then
         if (dtime .gt. 1.0E-5) then
            pnewdt = pnewdt_global
         end if
         !resets the value
         if (jelem .eq. 1) then
            open(103, file=trim(directory)//'/'//'pnewdt.txt')
            write(103,*) 100.        
            close(103)
         end if
         
      end if

c     load strengths:
      no_st_props = 4
      call load_st_props(st_props,no_st_props,
     &     props,nprops)

      
      delta_temp = predef(1,1,1)
      
      eltip_ind = 0
      do i=1,ntipels       
         if (jelem .eq. int(jelem_tipels(i,1))) then !element is a tip element
            eltip_ind = i
            call load_tip_el(svars,jelem_tipels,max_actelements,21,
     &           nsvars,i,nvars_m_xyz,nvars_penalty_f,
     &           nvars_dstat)

            exit
         end if
      end do
 
      
      do i=1,nnode/2                
         dstat(i)  = svars(nvars_dstat+i)
         penalty_f(i) = svars(nvars_penalty_f+i)
         do j=1,mcrd              
            m_xyz(i+(j-1)*nnode/2) = svars(nvars_m_xyz+i+(j-1)*nnode/2)
         end do
      end do
         
cccc  need to return rot_m
       call trans_coords(nlgeom,coords,mcrd,nnode,ndofel,u,du,
     &     rot_m,coords_lc,jelem) 
ccc 
c      call calc_rot_m_i(coords,rot_m,u,
c     &     mcrd,nnode,jelem,nlgeom)
       
      call nvcct_delta(u,ndofel,nnode,delta,rot_m)
      
      call ties_interface_nvcct_f(delta,nnode,mcrd,dstat,
     &     m_xyz,penalty_f,contact,
     &     rot_m,nforce,kmatrx,ku,ndofel,jelem)
      
      call approx_area_node(element,nintp,coords_lc,mcrd,nnode,area_n) !assumes all elements have the same area
     
      fc_thd = 0.5
      stress_crit_type = 1

      call stress_crit(nforce,st_props,no_st_props,nnode,mcrd,area_n,
     &     stress_crit_type,fc_thd,f_bool)

c      write(7,*) 'f_bool',f_bool
      if (f_bool .eq. .True.) then
         jelem_onset(jelem,1) = jelem
         do i=1,nnode/2*mcrd 
            jelem_onset(jelem,1+i) = nforce(i)
         end do
         jelem_onset(jelem,14) = area_n
        
      end if
      
      if ((eltip_ind .ne. 0)) then
         do i=1,nnode/2
            svars(nvars_dstat+i) = dstat(i)
         end do

         if (viz .eq. .True.) then
            do i=1, nnode/2
               jelem_dstat(jelem-jelem_prxvcct+1,i) = dstat(i)
            end do
         end if
        
         do i=1,nnode/2*mcrd
            jelem_force(eltip_ind,1+i) = nforce(i)
            jelem_delta(eltip_ind,1+i) = delta(i)               
         end do           
      end if
         
      do i = 1,ndofel
         do j=1,ndofel
            amatrx(i,j)=amatrx(i,j)+kmatrx(i,j)
         end do
         rhs(i,1)=rhs(i,1)-ku(i)
      end do
        
      return
      end
c********************************************************************************************
c********************* end subroutine uel ***************************************************
c********************************************************************************************

c*******************************************************************************
c********************* subroutine matrix_mul ***********************************
c*******************************************************************************
      subroutine matrix_mul(a,b,c,l,m,n)
c
      include 'aba_param.inc'
c
      dimension a(l,m), b(m,n), c(l,n)
c
      do i=1,l
        do j=1,n
          c(i,j) = 0.d0 !-initialize
          do k=1,m
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
          end do
        end do
      end do
c
      return
      end


      subroutine matrix_mul_new(a,b,c,l,m,n)
c
      implicit none
c
      integer j,i,k,l,m,n
      double precision a(l,m), b(m,n), c(l,n)
c
      do i=1,l
        do j=1,n
          c(i,j) = 0.d0 !-initialize
          do k=1,m
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
          end do
        end do
      end do
c
      return
      end

*****************************************************************************
      subroutine matrix_vect_mul(a,b,c,l,m)
      implicit none

      integer i,j,l,m
      double precision a(l,m), b(m), c(l)
c
      do i=1,l
          c(i) = 0.d0 !-initialize
          do j=1,m
            c(i) = c(i) + a(i,j)*b(j)
          end do
      end do
c
      return
      end

*****************************************************************************      
      subroutine matrix_scalar_mul(a,s,c,m,n)
      implicit none
      
!     inputs
      double precision a(m,n)
      double precision s
      integer m,n
      
!     outputs
      double precision c(m,n)
      
!     internal variables
      integer i,j
      
      do i=1,m
         do j=1,n
            c(i,j) = a(i,j)*s
         end do
      end do
     
      return
      end
*****************************************************************************      

      
      subroutine matrix_sym_a_at(a,c,m)
!     computes (a+a^t)/2.0
!     limited to squared matrices for now
      implicit none
      
!     inputs
      double precision a(m,m)
      integer m
!     internal variables
      double precision at(m,m)
      double precision ctemp(m,m)
!     outputs
      double precision c(m,m)
      
      call matrix_transpose_new(a,at,m,m)
      call matrix_sum(a,at,ctemp,m,m)
      call matrix_scalar_mul(ctemp,1.d0/2.d0,c,m,m)
      
      return
      end
      
*****************************************************************************
      subroutine matrix_sum(a,b,c,m,n)
      implicit none
      
!     inputs
      double precision a(m,n)
      double precision b(m,n)
      integer n,m

!     outputs
      double precision c(m,n)

!     internal variables
      integer i,j

      do i=1,m
         do j=1,n
            c(i,j)= a(i,j)+b(i,j)
         end do
      end do
                  

      return
      end
*****************************************************************************      
      
      
      
c*******************************************************************************
c********************* end subroutine matrix_mul *******************************
c*******************************************************************************  

c
c*******************************************************************************
c********************* subroutine invert ***************************************
c************** inverts a square matrix onto itself ****************************
c*******************************************************************************
      subroutine invert(ajacobi,detj,n)
c
      include 'aba_param.inc'
c
      dimension ajacobi(n,n),ajac(n,n)
c
      do i=1,n
        do j=1,n
          ajac(i,j)=ajacobi(i,j)
        end do
      end do
c
      if(n .eq. 2) then
        ajacobi(1,1)=ajac(2,2)
        ajacobi(2,1)=-ajac(2,1)
        ajacobi(1,2)=-ajac(1,2)
        ajacobi(2,2)=ajac(1,1)
      else if(n .eq. 3) then
        ajacobi(1,1)=ajac(2,2)*ajac(3,3)-ajac(3,2)*ajac(2,3)
        ajacobi(2,1)=ajac(3,1)*ajac(2,3)-ajac(2,1)*ajac(3,3)
        ajacobi(3,1)=ajac(2,1)*ajac(3,2)-ajac(3,1)*ajac(2,2)
        ajacobi(1,2)=ajac(3,2)*ajac(1,3)-ajac(1,2)*ajac(3,3)
        ajacobi(2,2)=ajac(1,1)*ajac(3,3)-ajac(3,1)*ajac(1,3)
        ajacobi(3,2)=ajac(3,1)*ajac(1,2)-ajac(1,1)*ajac(3,2)
        ajacobi(1,3)=ajac(1,2)*ajac(2,3)-ajac(2,2)*ajac(1,3)
        ajacobi(2,3)=ajac(2,1)*ajac(1,3)-ajac(1,1)*ajac(2,3)
        ajacobi(3,3)=ajac(1,1)*ajac(2,2)-ajac(2,1)*ajac(1,2)
      end if
c
      do i=1,n
        do j=1,n
          ajacobi(i,j)=ajacobi(i,j)/detj
        end do
      end do
c
      return
      end
c*******************************************************************************
c********************* end subroutine invert ***********************************
c*******************************************************************************

c*******************************************************************************
c********************* subroutine normalize ************************************
c********************* normalize a vector **************************************
c*******************************************************************************
      subroutine norm_vect(vect,ndim,vectnorm)
      implicit none

c     inputs
      integer ndim
      double precision vect(ndim)
c     internal variables
      integer i
      double precision vectlength
c     outputs
      double precision vectnorm(ndim)
      
      call calc_vect_length(vect,ndim,vectlength)

      if (vectlength .ne. 0.d0) then
         do i=1,ndim
            vectnorm(i) = vect(i)/vectlength
         end do
      else
         do i=1,ndim
            vectnorm(i) = 0.d0
         end do
      end if
      
      return
      end 
c*******************************************************************************
c********************* end subroutine norm_vect ********************************
c*******************************************************************************      




c*******************************************************************************
c********************* subroutine calc_vect_length *****************************
c********************* calculates the length of a vector ***********************
c*******************************************************************************
      subroutine calc_vect_length(vect,ndim,vectlength)
      implicit none

c     inputs
      integer ndim
      double precision vect(ndim)
c     internal variables
      integer i 
c     outputs
      double precision vectlength
      
      vectlength = 0.d0
      do i=1,ndim
         vectlength = vectlength + vect(i)**2.d0
      end do
      
      vectlength = sqrt(vectlength)

      return
      end 
c*******************************************************************************
c********************* end subroutine calc_vect_length *************************
c*******************************************************************************



c*******************************************************************************
c********************* subroutine matrix_transpose *****************************
c************************** matrix transpose ***********************************
c*******************************************************************************
      subroutine matrix_transpose(a,at,m,n)
c     
      include 'aba_param.inc'
      dimension a(m,n), at(n,m)
c     
      do i=1,m
         do j=1,n
            at(j,i) = a(i,j)
         end do
      end do
c     
      return
      end

      subroutine matrix_transpose_new(a,at,m,n)
      implicit none
c     inputs
      integer i,j,m,n
      double precision a(m,n), at(n,m)
c     
      do i=1,m
         do j=1,n
            at(j,i) = a(i,j)
         end do
      end do
c     
      return
      end subroutine


      subroutine initialize_matrix(a,ni,nj)
      implicit none
!     inputs
      integer ni
      integer nj
      
!     internal variables
      integer i,j

!     output
      double precision a(ni,nj)

      do i=1,ni
         do j=1,nj
            a(i,j) = 0.d0
         end do
      end do

      return
      end subroutine


c**************************************************************************
c******************* subroutine initialize_matrix_int *************************
c**************************************************************************
      subroutine initialize_matrix_int(a,ni,nj)
      implicit none
!     inputs
      integer ni
      integer nj
      
!     internal variables
      integer i,j

!     output
      integer a(ni,nj)

      do i=1,ni
         do j=1,nj
            a(i,j) = 0
         end do
      end do

      return
      end      

c**************************************************************************
c******************* subroutine initialize_vector *************************
c**************************************************************************
      subroutine initialize_vector(a,ni)
      implicit none
!     inputs
      integer ni

!     internal variables
      integer i

!     output
      double precision a(ni)
      
      do i=1,ni
            a(i) = 0.d0
      end do

      return
      end

c**************************************************************************
c******************* subroutine initialize_vector *************************
c**************************************************************************
      subroutine initialize_vector_int(a,ni)
      implicit none
!     inputs
      integer ni
!     internal variables
      integer i
!     output
      integer a(ni)
      
      do i=1,ni
            a(i) = 0
      end do

      return
      end subroutine 
   

c************************************************************************************
c*********************subroutine miner_acc_rule**************************************
c**determines the damage accumulation prior to failure, and number of cycles*********
c**required to fail at current stress level******************************************
c************************************************************************************
      subroutine miner_acc_rule(yt,s_const,maxps_prev,
     &     maxps_current,ncycles_prev,ncycles_current,
     &     c_miner,n_tofail,jelem) !included jelem to be able to assign a weak element

      implicit none
c     intput
      double precision yt,s_const !sn curve parameters (yt-static strength; s_const sn coefficient
      double precision ncycles_prev !total cycle count in the previous increment
      double precision ncycles_current !total cycle count in the current increment
      double precision maxps_prev !max principal stress at the sampling point in the previous increment
      double precision maxps_current !max principal stress at sampling point from the current increment
      integer jelem

c     internal
c     sn curve characterization:
      double precision n_inc !cycles to failure at stres level from previous increment
      double precision n_inc_current !cycles to failure at stres level from current increment

c     output
      double precision c_miner
      double precision n_tofail

c     computes cycles to failure at stress level from previous increment
      n_inc = 10**(1.d0/s_const*(1.d0-maxps_prev/yt))
      
c     computes accumulation of damage due to cycle accumulation since the previous increment
      c_miner = c_miner + (ncycles_current-ncycles_prev)/n_inc

c      computes cycles needed to fail at current stress level from the sn curve (witout accounting for damage accumulation)
      n_inc_current = 10.d0**(1/s_const*(1.d0-maxps_current/yt))
c      write(7,*) 'n_inc_current', n_inc_current
c     computes cycles needed to fail at current stress level, accounting to the miners rule for previous "damage" accumulation
      n_tofail = (1.d0-c_miner)*n_inc_current

      if (n_tofail .lt. 1.0e-14) then
         n_tofail = 1.0e-14
      end if
c      write(7,*) 'n_tofail',n_tofail

c$$$      if (n_to_fail .le. 0.d0) then
c$$$         n_to_fail = 1.d0/abs(n_to_fail)
c$$$      end if
      
c$$$      if (n_tofail .le. 1.d0) then
c$$$c$$$         write(7,*) 'warning: extrap leading to neg 
c$$$c$$$     &or 0 cycles to onset:', n_tofail
c$$$c$$$         write(7,*) 'jelem',jelem
c$$$c$$$          write(7,*) 'maxps_prev',maxps_prev
c$$$c$$$  write(7,*) 'maxps_current',maxps_prev 
c$$$c$$$          write(7,*) 's_const',s_const
c$$$c$$$          write(7,*) 'yt',yt 
c$$$c$$$          write(7,*) 'n_inc', n_inc
c$$$      
c$$$         n_tofail = 1.0d0
c$$$      end if

      return
      end


 
c*******************************************************************************
c******************* subroutine cross_product    **********************************
c*******************************************************************************

      subroutine cross_product(a,b,c,dim)
      
      implicit none
!     inputs
      integer dim
      double precision a(dim),b(dim)
!     output
      double precision c(dim)

      c(1) = a(2)*b(3)-a(3)*b(2)
      c(2) = -(a(1)*b(3)-a(3)*b(1))
      c(3) = a(1)*b(2)-a(2)*b(1)
      
      return
      end


      subroutine dot_product(a,b,dot_ab,dim)
      implicit none
      
      integer dim,i
      double precision a(dim)
      double precision b(dim)
      double precision dot_ab

      dot_ab = 0.d0
      do i=1,dim
         dot_ab = dot_ab+a(i)*b(i)
      end do
      
      return
      end subroutine
      
       
c******************************************************************************
c******************* subroutine calc_tie_contribution    **********************
c**calculates the contribution to the stiffness matrix from tieing ************
c**********************the floating nodes************************************** 
c******************************************************************************

      subroutine calc_tie_contribution(cnct,u,nnodes_tie,ndofel,
     &     ndim,kmatrxall,kuall,edget)

      implicit none
!     inputs
      integer cnct(2,3*ndim)
      integer nnodes_tie,ndofel,ndim
      double precision u(ndofel)
!     internal variables
      double precision krg
      double precision coef_v(ndim,3*ndim)
      double precision ut(3*ndim)
      double precision kmatrx(3*ndim,3*ndim), ku(3*ndim)
      double precision edget(2)
      integer i,j,k,l
      
!     output
      double precision kmatrxall(ndofel,ndofel)
      double precision kuall(ndofel)
      
    
      krg = 100000000.d0
      do k=1,nnodes_tie     
         call initialize_matrix(coef_v,ndim,3*ndim)
         
         call initialize_matrix(kmatrx,3*ndim,3*ndim)
         call initialize_vector(ku,3*ndim)
         do i=1,3
            coef_v(i,i) = 1.d0
            coef_v(i,i+3) = -(1.d0-edget(k))
            coef_v(i,i+6) = -edget(k)  
         end do

         do i=1,9
            do j=1,9
               do l=1,3
                  kmatrx(i,j) = kmatrx(i,j)+ krg*coef_v(l,i)*coef_v(l,j)
               end do
            end do
            ut(i) = u(cnct(k,i))
         end do

         call matrix_mul(kmatrx,ut,ku,9,9,1)
   
         do i=1,9
            do j=1,9
               kmatrxall(cnct(k,i),cnct(k,j)) =  
     &              kmatrxall(cnct(k,i),cnct(k,j))+
     &              kmatrx(i,j)  
            end do
            kuall(cnct(k,i))=kuall(cnct(k,i))+ku(i)
         end do
         
      end do
      
      return
      end
     
      subroutine load_tip_el(svars,jelem_tip,max_actelements,
     &     max_tip_vars,nvars,index,nvars_m_xyz,nvars_penalty_f,
     &     nvars_dstat)
      implicit none     
!     input
      double precision jelem_tip(max_actelements,max_tip_vars)
      integer max_actelements,max_tip_vars,nvars,index
      integer nvars_m_xyz,nvars_penalty_f
      integer nvars_dstat
!     output
      double precision svars(nvars)
!     internal variables
      integer j
      
      
      do j=1,4
         if (svars(nvars_dstat+j) .le. 0.001d0) then !node has not started to soften yet
            if  (jelem_tip(index,1+16+j) .gt. 0.001d0) then !node has started to soften in this increment update unloading slopes
               svars(nvars_m_xyz+j) = jelem_tip(index,1+j) !updates state variables: slopes
               svars(nvars_m_xyz+j+4) = jelem_tip(index,1+4+j) !updates state variables: slopes
               svars(nvars_m_xyz+j+8) = jelem_tip(index,1+8+j) !updates state variables: slopes
               svars(nvars_penalty_f+j) = jelem_tip(index,1+12+j) !updates state variables: shared nodes fraction
               svars(nvars_dstat+j) = jelem_tip(index,1+16+j) !updates state variables: dstat
            else                !element has not started to soften, update dstat (should be between 0 and 0.001)
               svars(nvars_dstat+j) = jelem_tip(index,1+16+j)
            end if
         else                   !element has started to soften
            if (svars(nvars_dstat+j).lt.
     &           jelem_tip(index,1+16+j)) then
               svars(nvars_dstat+j) = jelem_tip(index,1+16+j)
            end if
         end if
      end do


      
      return
      end
     

c******************** load_st_props *************************************
c**  "loads the strength properties for matrix cracks"                 **
c************************************************************************  
      subroutine load_st_props(st_props,no_st_props,
     &     props,nprops)
      
      implicit none
!     input
      integer nprops
      double precision props(nprops)
!     internal
      integer offset,i
!     output
      integer no_st_props
      double precision st_props(no_st_props)

c      st_props(1)  = yt
c      st_props(2)  = s_t !slope of the s-n curve
c      st_props(3)  = ys
c      st_props(4)  = s_s !slope of the s-n curve 

      offset = 6
      do i=1,no_st_props
         st_props(i) = props(offset+i)
      end do

      return
      end



c******************** load_e1_dir ****************************************
c**  "loads a reference vector for the ply orientations"                **
c*************************************************************************     
      subroutine load_e1_dir(ne1g,offset,
     &     props,nprops)

      implicit none
      double precision ne1g(3),props(nprops)
      integer offset,nprops,i

c     ne1g(1)  = component in global x
c     ne1g(2)  = component in global y
c     ne1g(3)  = component in global z

      do i=1,3
         ne1g(i) = props(offset+i)
      end do

      return
      end
     

      
c******************** det_fdir_vect ***********************
c**   "obtains fiber direction versor"                   **
c**********************************************************      
      subroutine det_fdir_vect(fdir_vect1,fdir_vect2,
     &     ne1g,thetam1,thetam2)

      implicit none
      
c     inputs
      double precision thetam1,thetam2
      double precision ne1g(3)
c     internal variables
      double precision cos_thetam1,sin_thetam1,cos_thetam2,sin_thetam2
      double precision pi
c     outputs
      double precision fdir_vect1(3),fdir_vect2(3)

      parameter (pi=3.14159265359d0)

      cos_thetam1=cos(pi*thetam1/180.d0)
      sin_thetam1=sin(pi*thetam1/180.d0)
      cos_thetam2=cos(pi*thetam2/180.d0)
      sin_thetam2=sin(pi*thetam2/180.d0)
      
      fdir_vect1(1) = ne1g(1)*cos_thetam1-sin_thetam1*ne1g(2)
      fdir_vect1(2) = ne1g(1)*sin_thetam1+cos_thetam1*ne1g(2)
      fdir_vect1(3) = ne1g(3)
      
      fdir_vect2(1) = ne1g(1)*cos_thetam2-sin_thetam2*ne1g(2)
      fdir_vect2(2) = ne1g(1)*sin_thetam2+cos_thetam2*ne1g(2)
      fdir_vect2(3) = ne1g(3)
  
      return
      end



c     ********************************************************
c     *********************check_perp*************************
c     ********************************************************
c     **checks if two unit vectors are perpendicular**********
c     ********************************************************
     
      subroutine check_perp(vect1,vect2,perp,ndim)
      implicit none

c     inputs
      integer ndim
      double precision vect1(ndim),vect2(ndim)

c     internal variables
      double precision inner_prod
      integer i

c     outputs
      logical perp
      
      perp = .true.
      
      inner_prod = 0.d0
      do i=1,3
         inner_prod = inner_prod + vect1(i)*vect2(i)
      end do
     
      if (abs(inner_prod) .gt. 0.1) then
         perp = .false.
      end if

      return
      end 


c     calculates a rotation matrix for the interface
c     key for nlgeom calculations
      subroutine calc_rot_m_i(coords,rot_m_i,u,
     &     ndim,nnode,jelem,nlgeom)
    
      
      implicit none
!     inputs
      integer nnode,ndim
      integer jelem
      logical nlgeom
      double precision coords(ndim,nnode),u(ndim*nnode)
!     internal variables
      integer i,j,k,nf,ndofel
      double precision e_x(3)
      double precision x(8),y(8),z(8),tri(ndim,ndim)
      double precision vx(ndim),vy(ndim),vz(ndim)
      double precision v_dummy(ndim),v_dummy2(ndim)
      double precision tri_out(ndim,ndim)
      double precision al
!     outputs
      double precision rot_m_i(ndim*nnode,ndim*nnode)
      
      ndofel = ndim*nnode 
      call initialize_matrix(rot_m_i,ndofel,ndofel)
      
c     compute middle surface to define element

c      write(7,*) 'coords',coords
      do i = 1, nnode/2     
         x(i) = coords(1,i)
         y(i) = coords(2,i)
         z(i) = coords(3,i)
      end do

      if (nlgeom == .True.) then
         do i=1,nnode/2
           x(i) = x(i) + 0.5d0*(u(3*(i-1)+1)+u(3*(i+3)+1))
           y(i) = y(i) + 0.5d0*(u(3*(i-1)+2)+u(3*(i+3)+2))
           z(i) = z(i) + 0.5d0*(u(3*(i-1)+3)+u(3*(i+3)+3))
        end do
      end if
c      write(7,*) 'u',u
      
      
c     vx' (node 2 - node 1):
      al = sqrt((x(1)-x(2))**2+(y(1)-y(2))**2+(z(1)-z(2))**2)
      vx(1) = (x(2)-x(1))/al
      vx(2) = (y(2)-y(1))/al
      vx(3) = (z(2)-z(1))/al
c      write(7,*) 'al',al
      
c     vy' (node 3 - node 2):         
      al = sqrt((x(3)-x(2))**2+(y(3)-y(2))**2+(z(3)-z(2))**2)
      vy(1) = (x(3)-x(2))/al
      vy(2) = (y(3)-y(2))/al
      vy(3) = (z(3)-z(2))/al
      
c     vz' - normal to the plane:      
      call cross_product(vx,vy,v_dummy,ndim)
      call norm_vect(v_dummy,ndim,vz)
c      write(7,*) 'vz',vz
      
c     x-direction:
      e_x(1) = 1.D0
      e_x(2) = 0.D0
      e_x(3) = 0.D0

c     compute projection of e_x onto the plane e_x_p:
c     v_x = vz'X'(e_x'X'vz)
      call cross_product(e_x,vz,v_dummy,ndim)
      call cross_product(vz,v_dummy,v_dummy2,ndim)
      call norm_vect(v_dummy2,ndim,vx)

c     calculate v_y, orthogonal to v_z and v_x
      call cross_product(vz,vx,v_dummy,ndim)
      call norm_vect(v_dummy,ndim,vy)  

c$$$      if (jelem .eq. 219) then
c$$$         write(7,*) 'vx',sqrt(vx(1)**2.0+vx(2)**2.0+vx(3)**2.0)
c$$$         write(7,*) 'vy',sqrt(vy(1)**2.0+vy(2)**2.0+vy(3)**2.0)
c$$$         write(7,*) 'vz',sqrt(vz(1)**2.0+vz(2)**2.0+vz(3)**2.0)
c$$$      end if
c$$$      
      do i=1,3
         tri(1,i) = vx(i)
         tri(2,i) = vy(i)
         tri(3,i) = vz(i)
      end do

      
c$$$      write(7,*) 'tri-JELEM:',JELEM
c$$$      do i =1,3
c$$$         do j=1,3
c$$$            write(7,'(f10.3,$)'),tri(i,j)
c$$$         end do
c$$$         write(7,*)
c$$$      end do
      
!     transformation matrix
      do i = 1, nnode              
         nf = ndim*(i-1)
         do j = 1, ndim
            do k = 1, ndim
               rot_m_i(nf+j,nf+k) = tri(j,k)
            end do
         end do
      end do

      return
      end






      
c********************** trans8n_nvc2 ****************************************
c* "calculate transformation matrices for converting                       **
c*  cz elements quantities between element and global coordinate systems   **
c*  in the presence of angled cracks such as migrated cracks"              **
c****************************************************************************       
      subroutine trans8n_nvc2(coords,ctrn,u,
     &     ndim,nnode,ndofel,jelem,fdir_vect,nlgeom)
      
      implicit none
!     inputs
      integer nnode,ndofel,ndim,i,j,k,nf
      logical nlgeom
      double precision x(8),y(8),z(8),tri(ndim,ndim)
      double precision vx(ndim),vy(ndim),vz(ndim)
      double precision v_dummy(ndim)
      double precision fdir_vect(ndim)
      double precision ctrn(ndofel,ndofel)
      double precision coords(ndim,nnode),u(ndofel)
      double precision al
      double precision pi
      double precision dev_angle,dot_p
      integer jelem,subel_int
      
      pi = acos(-1.0d0)
  
      call initialize_matrix ( ctrn, ndofel, ndofel )

      do i = 1, nnode/2     
         x(i) = coords(1,i)
         y(i) = coords(2,i)
         z(i) = coords(3,i)
      end do

      if (nlgeom .eq. .True.) then
         do i=1,nnode/2
            x(i) =x(i) + 0.5d0*(u(3*(i-1)+1)+u(3*(i+3)+1))
            y(i) =y(i) + 0.5d0*(u(3*(i-1)+2)+u(3*(i+3)+2))
            z(i) =z(i) + 0.5d0*(u(3*(i-1)+3)+u(3*(i+3)+3))
        end do
      end if   
      
c     vx' (node 2 - node 1):
      al = sqrt( (x(1)-x(2))**2+(y(1)-y(2))**2+(z(1)-z(2))**2 )
      vx(1) = (x(2)-x(1))/al
      vx(2) = (y(2)-y(1))/al
      vx(3) = (z(2)-z(1))/al
      
c     vy' (node 3 - node 2):         
      al = sqrt( (x(3)-x(2))**2+(y(3)-y(2))**2+(z(3)-z(2))**2 )
      vy(1) = (x(3)-x(2))/al
      vy(2) = (y(3)-y(2))/al
      vy(3) = (z(3)-z(2))/al
  
c     compute the orthogonal vector to the above to complete
c     the coordinate system definition. 
c     vz':
      call cross_product(vx,vy,v_dummy,ndim)
      call norm_vect(v_dummy,ndim,vz)
      
! confirm whether the cs is perpendicular to the fibers
      call dot_product(vz,fdir_vect,dot_p,ndim)
      
      dev_angle = acos(dot_p)/pi*180.d0 
      
      do i=1,3
         vx(i) = fdir_vect(i)
      end do
      
      call cross_product(vz,vx,v_dummy,ndim)
      call norm_vect(v_dummy,ndim,vy)
      call cross_product(vx,vy,v_dummy,ndim)
      call norm_vect(v_dummy,ndim,vz)

      do i=1,3
         tri(1,i) = vx(i)
         tri(2,i) = vy(i)
         tri(3,i) = vz(i)
      end do

!     transformation matrix
      do i = 1, nnode              
         nf = ndim*(i-1)
         do j = 1, ndim
            do k = 1, ndim
               ctrn(nf+j,nf+k) = tri(j,k)
            end do
         end do
      end do

      return
      end
      


  
      subroutine trans8n_i(coords,ctrn,tri,
     &     ndim,nnode,ndofel,jelem)
      
      implicit none                                       
      integer nnode,ndofel,ndim,i,j,k,nf
      logical mel
      double precision x(8),y(8),z(8),tri(ndim,ndim)
      double precision vx(ndim),vy(ndim),vz(ndim)
      double precision v_dummy(ndim)
      double precision tri_out(ndim,ndim)
      double precision fdir_vect(ndim)
      double precision ctrn(ndofel,ndofel)
      double precision ctrn_nvc(ndofel,ndofel),tri_nvc(ndim,ndim)
      double precision coords(ndim,nnode)
      double precision al
      double precision pi
      double precision mig_crack_ori,theta_i,theta_rot
      double precision theta_rot_ce
      double precision rot(ndim,ndim)
      double precision crack_type
      double precision dev_angle, m_angle, dot_p
      integer jelem,subel_int
      
      pi = acos(-1.0d0)
  
      call initialize_matrix ( ctrn, ndofel, ndofel )

c      mig_crack_ori = 0.d0

      
c     compute middle surface to define element
c     orientation

      do i = 1, nnode/2     
         x(i) = coords(1,i)
         y(i) = coords(2,i)
         z(i) = coords(3,i)
      end do

c     vx' (node 2 - node 1):
      al = sqrt( (x(1)-x(2))**2+(y(1)-y(2))**2+(z(1)-z(2))**2 )
      vx(1) = (x(2)-x(1))/al
      vx(2) = (y(2)-y(1))/al
      vx(3) = (z(2)-z(1))/al

      
c     vy' (node 3 - node 2):         
      al = sqrt( (x(3)-x(2))**2+(y(3)-y(2))**2+(z(3)-z(2))**2 )
      vy(1) = (x(3)-x(2))/al
      vy(2) = (y(3)-y(2))/al
      vy(3) = (z(3)-z(2))/al

      
c     compute the orthogonal vector to the above to complete
c     the coordinate system definition. 
c     vz':
      call cross_product(vx,vy,v_dummy,ndim)
      call norm_vect(v_dummy,ndim,vz)
      
c     compute vy such that it is perpendicular to vx
      call cross_product(vz,vx,v_dummy,ndim)
      call norm_vect(v_dummy,ndim,vy) 

      do i=1,3
         tri(1,i) = vx(i)
         tri(2,i) = vy(i)
         tri(3,i) = vz(i)
      end do
      
!     transformation matrix
      do i = 1, nnode              
         nf = ndim*(i-1)
         do j = 1, ndim
            do k = 1, ndim
               ctrn(nf+j,nf+k) = tri(j,k)
            end do
         end do
      end do

      return
      end
      

c*******************************************************************************
c********************** subroutine a_dmat **************************************
c*** assembles the d matrix which enables to pass from displacements to ********
c**  displacement jumps and vice-versa *****************************************
c*******************************************************************************
      subroutine a_dmat(dmat,ndofel,nnode)   
      implicit none
!     inputs
      integer ndofel, nnode
!     internal variables
      integer i,j
!     outputs
      double precision dmat(ndofel/2,ndofel)
    
      
      call initialize_matrix(dmat,ndofel/2,ndofel)
      
      do i=1,nnode/2
         do j=1,3
            dmat((i-1)*3+j,(i-1)*3+j) = -1.d0
            dmat((i-1)*3+j,(i-1)*3+nnode/2*3+j) = 1.d0
         end do
      end do
      
      return
      end
      
      
c**********************************************************************
c********************** subroutine a_dfddelta *************************
c***********assembles the jacobian in displacement space **************
c**********************************************************************
      subroutine k_secant(dstat,node_pair,delta_z,
     &     m_xyz,penalty_f,nnode,ndim,contact,ksec)
      
      implicit none
!     inputs
      integer node_pair,nnode,ndim
      logical contact
      double precision delta_z,dstat,penalty_f
      double precision m_xyz(ndim*nnode/2)
!     internal
      double precision k_c_user,k_pristine,d_ini
      double precision m_x,m_y,m_z
      double precision k_penalty_c,p_f
      integer i,j
!     outputs
      double precision ksec(nnode/2*ndim,nnode/2*ndim)

      k_c_user = 1.D0           !user defined factor that multiplies penalty stiffness in compression
      k_pristine = 100000000.D0
c     k_pristine = 10000000000.D0
c      k_pristine = 100000000000000.D0 !brit
      d_ini = 0.001D0
c      write(7,*) 'penalty_f',penalty_f
c      write(7,*) 'dstat',dstat
      if (penalty_f  .eq. 0.d0) then
         p_f = 1.D0
      else
         p_f = penalty_f
      end if
      m_x = m_xyz(node_pair)
      m_y = m_xyz(node_pair+nnode/2)
      m_z = m_xyz(node_pair+nnode)

      if (dstat .lt. d_ini) then
         do i=1,ndim
            ksec((node_pair-1)*ndim+i,(node_pair-1)*ndim+i) =
     &           k_pristine*p_f
         end do
         
      else if (dstat .ge. d_ini) then            
         if (m_x .eq. 0.d0) then
            m_x = -1000.D0
         end if 
         
         if (m_y .eq. 0.D0) then
            m_x = -1000.D0
         end if

         if (m_z .eq. 0.d0) then
            m_z = -1000.D0
         end if
         
         ksec((node_pair-1)*ndim+1,(node_pair-1)*ndim+1)=p_f*
     &        (m_x*(1.d0-1.d0/dstat))
         
         ksec((node_pair-1)*ndim+2,(node_pair-1)*ndim+2)=p_f*
     &        (m_y*(1.d0-1.d0/dstat))
         
         if (contact .eq. .true.) then
            if (delta_z .ge. 0.d0) then
               ksec((node_pair-1)*ndim+ndim,
     &              (node_pair-1)*ndim+ndim)=
     &              penalty_f*(m_z*(1.d0-1.d0/dstat))
            else              
               k_penalty_c = k_pristine*k_c_user
               
               ksec((node_pair-1)*ndim+ndim,(node_pair-1)*ndim+ndim) = 
     &              p_f*k_penalty_c
            end if
         else
            ksec((node_pair-1)*ndim+ndim,(node_pair-1)*ndim+ndim) = 
     &           p_f*(m_z*(1.d0-1.d0/dstat))
         end if
         
      end if


c$$$      write(7,*) '---in ksec:---'
c$$$      do i =1,nnode/2*ndim
c$$$         do j=1,nnode/2*ndim
c$$$            write(7,'(f10.3,$)'),ksec(i,j)
c$$$         end do
c$$$         write(7,*)
c$$$      end do
      
      
      return
      end                 
      
c*******************************************************************************
c********************** subroutine ties_interface ******************************
c*******************************************************************************
      subroutine ties_interface_nvcct_f(delta,nnode,ndim,dstat,
     &     m_xyz,penalty_f,contact,rot_m,nforce,kmatrx,ku,ndofel,jelem)    
       implicit none
!     inputs
       integer nnode,ndim,ndofel
       logical contact
       double precision dstat(nnode/2),delta(ndofel/2)
       double precision rot_m(ndofel,ndofel)
       double precision m_xyz(ndofel/2)
       double precision penalty_f(nnode/2)
       integer jelem
       
!     internal variables
       integer i,j
       double precision dmat(ndofel/2,ndofel)
       double precision dmatt(ndofel,ndofel/2)
       double precision tmp1(ndofel,ndofel)
       double precision ksec_dmat(ndofel/2,ndofel)
       double precision ksec(ndofel/2,ndofel/2)
       double precision nforce(ndofel/2)
       double precision lkmatrx(ndofel,ndofel),lku(ndofel)
       double precision kmatrx(ndofel,ndofel),ku(ndofel)
       double precision rot_mt(ndofel,ndofel)
       
       
       call initialize_vector(nforce,ndofel/2)
       call a_dmat(dmat,ndofel,nnode) !assemble the dmat matrix
       call initialize_matrix(ksec,ndofel/2,ndofel/2) 
       
       do i=1,nnode/2
!     assembles dfddelta - jacobian in disp jump space
          call k_secant(dstat(i),i,delta(ndim+(i-1)*ndim),
     &         m_xyz,penalty_f(i),nnode,ndim,contact,ksec) 
       end do

!     nforce:
       call matrix_vect_mul(ksec,delta,nforce,ndofel/2,
     &     ndofel/2)
       
!     determine local stiffness matrix lkmatrx      
       call matrix_mul_new(ksec,dmat,ksec_dmat,
     &      ndofel/2,ndofel/2,ndofel)

       call matrix_transpose_new(dmat,dmatt,ndofel/2,ndofel)
       call matrix_transpose_new(rot_m,rot_mt,ndofel,ndofel)
       call matrix_mul_new(dmatt,ksec_dmat,lkmatrx,ndofel,
     &      ndofel/2,ndofel)
!     determine local force vector lku          
       call matrix_vect_mul(dmatt,nforce,lku,ndofel,ndofel/2)

!     transforms local strifness matrix to global
       call matrix_mul(lkmatrx,rot_m,tmp1,
     &      ndofel,ndofel,ndofel)
       call matrix_tmul(rot_m,tmp1,kmatrx,ndofel,ndofel,
     &      ndofel)
!     transforms local force vector to global
       call matrix_tmul(rot_m,lku,ku,ndofel,ndofel,1)

       
       return
       end
      
c*******************************************************************************
c********************** subroutine nvcct_force ******************************
c*******************************************************************************

      subroutine nvcct_force(kuall,ndofel,nnodes,ndim,cnci,n_force_i_lc,
     &     rot_m_i)
      

      implicit none
c     input
      double precision rot_m_i(nnodes*ndim,nnodes*ndim)
      double precision kuall(ndofel)
      integer ndofel,ndim,nnodes,cnci(nnodes*ndim)
      
c     internal
      integer i,j,ndofel_i
      double precision rot_m_i_s(nnodes/2*ndim,nnodes/2*ndim)
      double precision n_force_i_lc(nnodes/2*ndim)
      double precision n_force_i(nnodes/2*ndim)
      
c     output

      ndofel_i = nnodes*ndim
      do i=1,ndofel_i/2
         do j=1,ndofel_i/2
            rot_m_i_s(i,j) = rot_m_i(i,j)
         end do
      end do
      
      do i=1,ndofel_i/2
         n_force_i(i) = kuall(cnci(i))
      end do

      call matrix_mul(rot_m_i_s,n_force_i,n_force_i_lc,
     &    ndofel_i/2,ndofel_i/2,1)
      

      return
      end


      subroutine nvcct_delta(u,ndofel,nnode,delta_i,rot_m_i)

      implicit none
c     input
      integer ndofel,nnode
      double precision rot_m_i(ndofel,ndofel)
      double precision u(ndofel)
c     internal
      integer i,j,ndofel_i
      double precision u_r(ndofel)
      double precision dmat(ndofel/2,ndofel)
c     output
      double precision delta_i(ndofel/2)
      


      
      call matrix_vect_mul(rot_m_i,u,u_r,ndofel,ndofel)
      call a_dmat(dmat,ndofel,nnode) !assemble the dmat matrix       
      call matrix_vect_mul(dmat,u_r,delta_i,ndofel/2,ndofel)
     
      return
      end

c*********************  matrix_tmul **********************************
c**  "c = a^t*b"                                                    **
c*********************************************************************
      subroutine matrix_tmul(a,b,c,l,m,n)

      include 'aba_param.inc'

      dimension a(l,m), b(m,n), c(l,n)

      do i=1,l
        do j=1,n
          c(i,j) = 0.d0 !-initialize
          do k=1,m
            c(i,j) = c(i,j) + a(k,i)*b(k,j)
          end do
        end do
      end do

      return
      end 
      
c*********************  tract_split_el *******************************
c**   "obtains tractions in the plane defined by xyzei"             **
c**   based on the forces nforce in an element that has been split  **
c**   but has not failed"                                           **
c*********************************************************************
      subroutine tract_split_el(nforce,xyzei,tract)

      implicit none
!     input
      double precision nforce(4,3),xyzei(3,8)
!     internal variables
      double precision avg_nforce(3),side1(3),side2(3)
      double precision side1_le, side2_le, area
      integer i,j
!     output variables
      double precision tract(3)
      
      do i=1,4
         do j=1,3
            avg_nforce(j) = avg_nforce(j) + nforce(i,j)
         end do
      end do
      
      do j=1,3
         avg_nforce(j) = avg_nforce(j)/4.d0         
      end do
  
      do j=1,3
         side1(j) = xyzei(j,1)-xyzei(j,2)
         side2(j) = xyzei(j,2)-xyzei(j,3)
      end do

      call norm_vect(side1,3,side1_le)
      call norm_vect(side2,3,side2_le)
      area = side1_le*side2_le
      
      do j=1,3
         tract(j) = avg_nforce(j)/area
      end do

      return
      end 



      subroutine override_sn_off_i(jelem,stress_bk,strength_bk,giic,e22,
     &     area_el_i,n_tofail)

      implicit none
      
!     inputs
      integer jelem
      double precision tract(3) !traction
      double precision stress_bk,strength_bk !strenght (having into account mode-mixity)
      double precision e22,giic !e22 - young's modulus; giic - mode ii fracture toughness (worst case)
      double precision giic_dummy
      double precision area_el_i
      double precision pi
!     internal
      double precision clength  !characteristic lenght
      double precision alpha
!     output
      double precision n_tofail

      pi = acos(-1.0d0)
      
      giic_dummy = 1.0

      clength = sqrt(area_el_i/4.D0)/4.0
      
      alpha = sqrt(giic_dummy*e22/(2*pi*clength))

c$$$      write(7,*) 'jelem',jelem
c$$$      write(7,*) 'clength',clength
c$$$      write(7,*) 'stress_bk',stress_bk
c$$$      write(7,*) 'giic,e22',giic,e22
c$$$      write(7,*) 'alpha',alpha
c$$$      write(7,*) 'strength_bk',strength_bk
c$$$      write(7,*) 'stress_bk - alpha',stress_bk-alpha
      
      stress_bk = stress_bk - alpha

      if (stress_bk .lt. 0.d0) then
         stress_bk = 0.0
      end if
      
      return
      end



      

      subroutine read_viz_txt(directory,file,viz,jelem_pligcoe,
     &     jelem_ghost)
!     inputs
      character(len=*) directory
      character(len=*) file 
!     outputs
      logical viz
      integer jelem_pligcoe,jelem_ghost
      
      inquire(file=trim(trim(directory)//'/'//file),
     &     exist=viz)
      if (viz .eq. .False.) then
         write(7,*) '*********vis.txt not provided'
      else
         write(7,*) '*********vis.txt provided:'
         write(7,*) 'ensure part using czes starts with underscore'
         write(7,*) 'e.g. _part; that element numbers within'
         write(7,*) 'the part are sequential starting in 1 with no'
         write(7,*) 'gaps; the the first pligcoe eleement is given' 
         write(7,*) 'in the first line and the the first ghost/dummy'
         write(7,*) 'element is given in the second'
         write(7,*) '***********************************************'
         open(103, file=trim(directory)//'/'//file)
         read(103,*) jelem_pligcoe
         read(103,*) jelem_ghost
      end if

      return
      end subroutine read_viz_txt


      
      subroutine approx_area_node(element,nintp,coords_lc,ndim,
     &     nnode,area_n)

      implicit none
!     inputs:
      integer ndim,nintp,nnode
      double precision coords_lc(ndim-1,nnode/2)
      character(len=15):: element
      
!     internal variables
      integer i
      double precision fn(nnode/2),dn(nnode/2,ndim-1)
      double precision w(nintp),intp_nat(ndim-1,nintp)
      double precision detj
      
!     outputs
      double precision area_n
      
      area_n = 0.D0
      call initialgauss_ce(element,ndim,intp_nat,w,nintp)
      do i=1,nintp
         call shape_ce(element,nnode/2,ndim,intp_nat,nintp,i,fn,dn)
         call detj_ce(ndim,nnode,coords_lc,dn,detj)
         area_n = area_n+detj
      end do

      area_n = area_n/(nnode/2)
      
      return
      end subroutine


      subroutine stress_crit(nforce,st_props,no_st_props,nnode,ndim,
     &     area_n,stress_crit_type,stress_thd,f_bool)

      implicit none
!     inputs:
      integer stress_crit_type !1 - quad_stresss
      integer nnode,ndim,no_st_props      
      double precision area_n,stress_thd
      double precision nforce(nnode/2*ndim)
      double precision st_props(no_st_props)
     
!     internal variables
      integer i,j
      double precision f_c,Z_T,Z_S,Z_st
      double precision stress
!     outputs
      logical f_bool


      Z_T = st_props(1)
      Z_S = st_props(3)
      Z_st = Z_S
      f_bool = .False.
      if (stress_crit_type .eq. 1) then !quad stress crit
         do i=1, nnode/2
            f_c = 0.d0
            do j=1,3
               stress = (nforce(j+(i-1)*ndim)/area_n)
               if (j .eq. 3) then
                  if (stress .lt. 0.d0) then
                     stress = 0.D0
                  end if
                  Z_st = Z_T 
               end if
               f_c=f_c+(stress/Z_st)**2.D0  
            end do
c$$$            write(7,*) 'nforce',nforce
c$$$            write(7,*) 'area_n',area_n
c$$$            write(7,*) 'Z_st',Z_st
            
            if (f_c .gt. stress_thd) then
c               write(7,*) 'f_c',f_c
               f_bool = .True.
            end if
         end do
      end if

      
      return
      end subroutine


      subroutine shape_ce(element,nnode,ndim,pq,no_intps,i_intp,f,df)
c     shape function/derivatives evaluation for bi-linear tri/quad cohesive elements
      implicit none
      
c     inputs
      integer nnode,ndim,i_intp,no_intps
      double precision pq(ndim-1,no_intps)
      character(len=15):: element
c     internal variables
      integer j,k
      double precision xi,eta,etam,etap,xim,xip
c     outputs
      double precision f(nnode), df(nnode,ndim-1)

      
      do k=1,nnode
         f(k)=0.d0
         do j=1,ndim-1
            df(k,j)=0.d0
         end do
      end do


      if (element.eq.'line') then
         xi = pq(1,i_intp)
         f(1) = 0.5d0*(1.d0-xi)
         f(2) = 0.5d0*(1.d0+xi)
         df(1,1)= -0.5d0
         df(2,1)= 0.5d0
      else
         xi=pq(1,i_intp)
         eta=pq(2,i_intp)
         if(element.eq.'triangle') then
            f(1) = 1.d0-xi-eta
            f(2) = xi
            f(3) = eta

            df(1,1) = -1.d0
            df(2,1) = 1.d0
            df(3,1) = 0.d0
            df(1,2) = -1.d0
            df(2,2) = 0.d0 
            df(3,2) = 1.d0
         else if(element.eq.'quadrilateral') then        
            etam= 0.25d0*(1.0-eta)
            etap= 0.25d0*(1.0+eta)
            xim = 0.25d0*(1.0-xi)
            xip = 0.25d0*(1.0+xi)

            f(1)=0.25d0*(1.d0-xi)*(1.d0-eta)
            f(2)=0.25d0*(1.d0+xi)*(1.d0-eta)
            f(3)=0.25d0*(1.d0+xi)*(1.d0+eta)
            f(4)=0.25d0*(1.d0-xi)*(1.d0+eta)

            df(1,1) = -etam
            df(2,1) =  etam
            df(3,1) =  etap
            df(4,1) = -etap
            df(1,2) = -xim
            df(2,2) = -xip
            df(3,2) =  xip
            df(4,2) =  xim
         end if
      end if
         
      return
      end

      subroutine detj_ce(ndim,nnode,xyzei,dn,detj)
c     calculate determinant for intpoints of a cohesive element
      implicit none
!     inputs
      integer ndim,nnode
      double precision xyzei(ndim-1,nnode/2)
      double precision dn(nnode/2,ndim-1)
!     internal variables
      double precision jac(ndim-1,ndim-1)
!     outputs
      double precision detj
      
      call matrix_mul(xyzei,dn,jac,ndim-1,nnode/2,ndim-1)
      if (ndim .eq. 3) then
         call determinant(jac,detj,ndim-1) 
      else
         detj = jac(1,1)
      end if
      
      return
      end subroutine

      
      
      subroutine determinant(ajacobi,detj,n)
c     determinant of a 1x1 2x2 3x3 matrix
c     inputs
      integer n
      double precision ajacobi(n,n)
c     output
      double precision detj

      if (n .eq. 1) then
        detj= 1.d0
      else if (n .eq. 2) then
        detj= ajacobi(1,1)*ajacobi(2,2) - ajacobi(1,2)*ajacobi(2,1)
      else if (n .eq. 3) then
        detj= ajacobi(1,1)*(ajacobi(2,2)*ajacobi(3,3)-
     +            ajacobi(3,2)*ajacobi(2,3))
        detj= detj- ajacobi(1,2)*(ajacobi(2,1)*ajacobi(3,3)-
     +            ajacobi(3,1)*ajacobi(2,3))
        detj= detj+ ajacobi(1,3)*(ajacobi(2,1)*ajacobi(3,2)-
     +            ajacobi(3,1)*ajacobi(2,2))
      end if
c
      return
      end subroutine

      
      
      subroutine initialgauss_ce(element,ndim,xyzg,wt,nig)
c     initialize gauss integrating points
c     returns the local coordinates of the integrating points and weights
      implicit none
c     inputs
      character(len=15):: element
      integer nig     
c     internal variables
      integer i,j,ndim
      double precision zero,one,two,three,four,six,eight
      double precision root3,one_half,one_third,two_third,one_sixth
c     outputs
      double precision xyzg(ndim-1,nig),wt(nig)
   
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0,
     &     six=6.d0,eight=8.d0)
      
      root3= one/sqrt(three)
      one_half=one/two
      one_third= one/three
      two_third= two/three
      one_sixth= one/six

      do i=1,nig
         wt(i)=zero
         do j=1,ndim-1
            xyzg(j,i)=0.D0
         end do
      end do

      if (element .eq. 'line') then
         if (nig .eq. 2) then
            xyzg(1,1) = -root3
            xyzg(1,2) = root3
            do i=1,nig
               wt(i) = one
            end do
         end if
      else if (element .eq. 'triangle') then
        if (nig .eq. 1) then
          xyzg(1,1)= one_third
          xyzg(2,1)= one_third
          wt(1) = one_half
        else if (nig .eq. 3) then
          xyzg(1,1)= one_half
          xyzg(2,1)= one_half
          xyzg(1,2)= one_half
          xyzg(2,2)= zero
          xyzg(1,3)= zero
          xyzg(2,3)= one_half
          wt(1) = one_sixth
          wt(2) = wt(1)
          wt(3) = wt(1)
        end if
      else if (element .eq. 'quadrilateral') then
        if (nig .eq. 1) then
          xyzg(1,1)= zero
          xyzg(2,1)= zero
          wt(1) = four
        else if (nig .eq. 4) then
           xyzg(1,1)= -root3
           xyzg(2,1)= -root3
           xyzg(1,2)=  root3
           xyzg(2,2)= -root3
           xyzg(1,3)=  root3
           xyzg(2,3)=  root3
           xyzg(1,4)= -root3
           xyzg(2,4)=  root3
          do i=1,nig
            wt(i) = one
          end do
        end if     
      end if

      return
      end


      subroutine trans_coords(nlgeom,coords,ndim,nnodes,ndofel,u,du,
     &     ctrn,coords_lc,jelem)
c     calculate transformation matrices for converting cz elements quantities
c     between element and global coordinate systems
      implicit none
!     inputs:
      integer nnodes,ndim,ndofel
      integer jelem
      double precision coords(ndim,nnodes),u(ndofel),du(ndofel)
      logical nlgeom
!     internal variables:
      integer k,j
      double precision x_n(ndim,nnodes/2)
      double precision tri(ndim,ndim)
!     outputs:
      double precision coords_lc(ndim-1,nnodes/2)
      double precision ctrn(ndofel,ndofel)

      call mid_surface_coords_ce(nlgeom,u,du,coords,ndim,ndofel,nnodes,
     &     x_n)
     
      call ce_rot_matrix_node(x_n,nnodes,ndim,tri)

       if ((jelem .eq. 10001) .or. (jelem .eq. 10480)) then
          write(7,*) 'jelem',jelem
          write(7,*) 'tri',tri
          write(7,*)'x_n',x_n
       end if
      
      call map_rot_matrix_ce(tri,ndim,nnodes,ndofel,ctrn)

      call lnodes_csurface_coords(tri,x_n,ndim,nnodes,coords_lc)
      
      return
      end subroutine

      subroutine lnodes_csurface_coords(tri,x_n,
     &     ndim,nnodes,xyze_l)
      implicit none
!     inputs
      integer ndim,nnodes
      double precision tri(ndim,ndim)
      double precision x_n(ndim,nnodes/2)
!     internal variables
      integer i,j
      double precision  xyz_g_n(ndim),xyz_l_n(ndim)
!     outputs
      double precision xyze_l(ndim-1,nnodes/2)
!
      do i=1,nnodes/2
       
         do j=1,ndim
            xyz_g_n(j) = x_n(j,i)
         end do
         call matrix_vect_mul(tri,xyz_g_n,xyz_l_n,ndim,ndim)
         
         do j=1,ndim-1
            xyze_l(j,i) = xyz_l_n(j)
         end do            
      end do
      
      return
      end subroutine



      
      subroutine mid_surface_coords_ce(nlgeom,u,du,coords,ndim,ndofel,
     &     nnodes,x_n)
      implicit none
!     inputs
      integer ndofel,ndim,nnodes
      double precision coords(ndim,nnodes)
      double precision u(ndofel),du(ndofel)
      logical nlgeom
!     internal variables
      integer i,j, nnodes_bots
      double precision u_t(ndofel)
!     output
      double precision x_n(ndim,nnodes/2)
      
      do i=1,ndofel
         u_t(i) = u(i) - du(i)
      end do

      nnodes_bots = nnodes/2
      if (nlgeom .eq. .True.) then 
         do i = 1, nnodes/2
            do j=1,ndim
               x_n(j,i) = coords(j,i) +
     &              0.5d0*(u_t(ndim*(i-1)+j)+
     &              u_t(ndim*(i+nnodes_bots-1)+j))
            end do
         end do
      else
         do i = 1, nnodes/2
            do j=1,ndim
               x_n(j,i) = coords(j,i)
            end do
         end do
      end if

            
      return
      end subroutine


      
      subroutine ce_rot_matrix_node(x_n,nnodes,ndim,
     &     tri)
      implicit none
!     inputs
      integer ndim,nnodes
      double precision x_n(ndim,nnodes/2)
!     internal variables
      integer i
      double precision vx(ndim),vy(ndim),vz(ndim)
      double precision vx_n(ndim),vy_n(ndim),vz_n(ndim)     
!     outputs
      double precision tri(ndim,ndim)

c     vx (node 2 - node 1):
      do i=1,ndim
         vx(i) = x_n(i,2) - x_n(i,1) 
      end do
      call norm_vect(vx,ndim,vx_n)
       
      if (ndim .eq. 2) then
c     compute normal vector 2 the crack (2D), VZ:
         vz_n(1) = -vx_n(2)
         vz_n(2) = vx_n(1)
         
         do i=1,ndim
            tri(1,i) = vx_n(i)
            tri(2,i) = vz_n(i)
         end do
         
      else if (ndim .eq. 3) then
         
c     vy (node 3 - node 2):         
         do i=1,ndim
            vy(i) = x_n(i,3) - x_n(i,2) 
         end do
         call norm_vect(vy,ndim,vy_n)
         
c     compute normal vector 2 crack (3D), VZ
         call cross_product(vx_n,vy_n,vz,ndim)
         call norm_vect(vz,ndim,vz_n)
         
c     re-obtain VY, normal to VX and VZ:
         call cross_product(vz_n,vx_n,vy,ndim)
         call norm_vect(vy,ndim,vy_n)
         
         do i=1,ndim
            tri(1,i) = vx_n(i)
            tri(2,i) = vy_n(i)
            tri(3,i) = vz_n(i)
         end do
      end if
      
      return
      end subroutine


      
      subroutine map_rot_matrix_ce(tri,ndim,nnodes,ndofel,ctrn)
      implicit none
!     inputs
      integer ndim,ndofel,nnodes
      double precision tri(ndim,ndim)
!     internal variables
      integer nf,i,j,k
!     outputs
      double precision ctrn(ndofel,ndofel)

      call initialize_matrix(ctrn,ndofel,ndofel)
      
      do i = 1, nnodes              
         nf = ndim*(i-1)
         do j = 1, ndim
            do k = 1, ndim
               ctrn(nf+j,nf+k) = tri(j,k)
            end do
         end do
      end do

      return
      end subroutine


      subroutine nintp_points_el_type(element,nintp)
c     returns number of integration points based on element type
      implicit none
c     inputs
      character(len=15):: element
c     output
      integer nintp

      nintp = 0
      if (element .eq. 'line') then
         nintp = 2
      else if (element .eq. 'triangle') then
         nintp = 3
      else if (element .eq. 'quadrilateral') then
         nintp = 4
      end if
      
      return
      end subroutine


      
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     1                drplde,drpldt,stran,dstran,time,dtime,temp,
     2                dtemp,predef,dpred,cmname,ndi,nshr,ntens,
     3                nstatv,props,nprops,coords,drot,pnewdt,
     4                celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,
     5                kstep,kinc)

      include 'aba_param.inc'


      character*80 cmname  
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1          ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),
     2          time(2),predef(1),dpred(1),props(nprops),coords(3),
     3          drot(3,3),dfgrd0(3,3),dfgrd1(3,3)         
      double precision strain(6),stiff(6,6),ammss(6,6),alpha(2),beta(2)
      double precision e1,e2,g12,v12,v23,v13,dm
      integer offset
      parameter (zero = 0.d0, one = 1.d0, two = 2.d0)

c     ***common block 3*** 
      double precision jelem_dstat(200000,4)
      integer jelem_prxvcct
      integer jelem_ghost
      logical viz
      common /dstat_common/jelem_dstat,jelem_prxvcct,
     &jelem_ghost,viz     
c     ***end common block***      
      double precision dstat(4)

      do i=1, 4
         dstat(i) = jelem_dstat(noel-jelem_ghost+1,i)
      end do

      statev(1) = maxval(dstat)
      statev(2) = minval(dstat)
      
      return
      end subroutine umat
