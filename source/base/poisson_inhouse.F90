module poisson_inhouse
  use kind_parameters
  use common_parameter
  use common_2d
  use omp_lib
  use svd_lib
  !use sparse_lsys_inhouse
  use linear_solver
  use analytic_functions
  implicit none

  !! TO DO:
  !! Parallelise Poisson building (as in my ISPH code)
  

contains
  subroutine solve_poisson_inhouse(phi,rhs_ftn)
    real(rkind),dimension(:),intent(out) :: phi
    real(rkind),dimension(:),intent(in) :: rhs_ftn
    integer(ikind) :: i,j,k,ij,ii,num,n_t
    real(rkind) :: temp,tmp
    real(rkind),dimension(dims) :: riii
    integer(ikind) :: iters,itermax
    real(rkind) :: residual,tol,x,y
    real(rkind), dimension(:,:), allocatable :: l_mat
  
    real(rkind), dimension(:), allocatable :: rhs_vec,sol_vec
   ! real(rkind), dimension(:,:), allocatable :: l_mat
    integer(ikind), dimension(:),allocatable :: ij_num
    real(rkind),dimension(:),allocatable :: p_fs 
    logical :: rebuild_lhs
  
    real(rkind) :: tt1,tt2,tt3,tt4

    allocate(rhs_vec(npfb+nb_n),sol_vec(npfb+nb_n));rhs_vec = 0.0d0;sol_vec=0.0d0
    allocate(lhs_mat((npfb+nb_n)*nplink));lhs_mat = 0.0d0
    allocate(ija((npfb+nb_n)*nplink));ija = 0
    allocate(l_mat(np,np));l_mat=0.0d0


    !! Build the RHS
    !$OMP PARALLEL DO
    do i=1,npfb
       rhs_vec(i) = rhs_ftn(i)  
       sol_vec(i)=ftn(rp(i,1),rp(i,2))
    end do
    !$OMP END PARALLEL DO
    if(nb_n.ne.0)then
       !! Homogeneous Neumann BC
       do i=npfb+1,npfb+nb_n
          rhs_vec(i)=0.0d0
       end do
    end if

    if(.true.)then  !! Conditional in case I want to avoid rebuilding LHS
       !! Make ija which is some indexing variable...
       allocate(ij_num(npfb+nb_n));ij_num = 0
       ija(1) = npfb + nb_n + 2
       num = npfb + nb_n + 1
       do i=1,npfb + nb_n
          if(i.le.npfb) then
             ii=i
          else
             ii=irelation(i)
          end if
          do k=1,ij_count(ii)
             j=ij_link(ii,k)
             if(j.le.npfb + nb_n) then
                num = num + 1
                ij_num(j) = num 
                ija(num)=j
             end if
             if(j.gt.npfb + nb_n) then
                num=num+1
                ij = irelation(j)
                ij_num(ij) = num
                ija(num) = ij 
             end if
          end do
          ija(i+1) = num + 1
          !! Build the LHS
          if(i.le.npfb)then
          do k=1,ij_count(i)
             j=ij_link(i,k)
             temp = -ij_w_lap(i,k)
             if(i.ne.j)then
                lhs_mat(i) = lhs_mat(i) + temp
                if(j.gt.npfb + nb_n)then    ! If j is mirror, then
                   ij=irelation(j)    ! look for parent of j, ij
                   if(i.ne.ij)then
                      n_t = ij_num(ij)
                      lhs_mat(n_t)=lhs_mat(n_t)-temp 
                   else
                      lhs_mat(i) = lhs_mat(i) - temp
                   end if
                elseif(j.le.npfb + nb_n)then   ! if j is not a mirror or a boundary node
                   n_t = ij_num(j)
                   lhs_mat(n_t) = lhs_mat(n_t) - temp
                end if
             end if
             j=ij_link(i,k)
          if (j>npfb) then
            j=irelation(j)
            endif
         ! print *, ij_num(j)
          l_mat(i,j)=lhs_mat(ij_num(j))
          end do
          else   !! Neumann BC - not relevant at present, as only simple
          !! domains available in this version of code.
             ii=irelation(i)    !! the boundary particle ii
             riii = rp(i,:)-rp(ii,:)
             tmp = sqrt(dot_product(riii,riii))
             riii=riii/tmp   !! The unit normal vector...
             do k=1,ij_count(ii)   !! for all particles in stencil of ii
                j=ij_link(ii,k)
                temp = dot_product(ij_w_grad(ii,k,:),riii) !! The normal derivative weighting
                if(ii.ne.j) then !! contribution of j added to ii
                   n_t=ij_num(ii)
                   lhs_mat(n_t) = lhs_mat(n_t) + temp
                   if(j.le.npfb+nb_n)then
                      n_t=ij_num(j)
                      lhs_mat(n_t)=lhs_mat(n_t) - temp
                   else  
                   !! For now, assuming no solid ghost is near a non-solid ghost
                   end if             
                end if
             end do
          end if
       end do
  !     deallocate(ij_num)
    end if
    
 !   print *, ij_num
 !   stop
 !  print *, nb
 !  stop
    !! Build Dirichlet BCs if required
    if(nb.ne.0)then
       do i=1,max(1,nb)
          lhs_mat(i)=1.0d12
          rhs_vec(i)=0.0d0
       end do
    end if
    
    do i=1, npfb
       l_mat(i,i)=lhs_mat(i)
       
       !print *,i
      ! do k=2, ij_count(i)
       
          
      ! end do
  !     print *, l_mat(i,1:npfb), rhs_vec(i),rp(i,1), rp(i,2)
    end do
    
  !  stop
!stop
    !! Set tolerance and max iterations allowed
    tol = 1d-10;itermax=3000
!print *, rhs_vec
!stop
    !! Solve the equation
    if(nb.ne.0)then !! If there are any Dirichlet boundary nodes
       call sparse_pbicgstab(lhs_mat,ija,sol_vec,rhs_vec,npfb+nb_n,(npfb+nb_n)*nplink,itermax,tol,iters,residual)
    else    !! If Neumann/periodic BCs only
    
       call sparse_pbicgstab_nonull(lhs_mat,ija,sol_vec,rhs_vec,npfb+nb_n,(npfb+nb_n)*nplink,itermax,tol,iters,residual)
    end if
    deallocate(ija,lhs_mat,rhs_vec)

    !! Pass solution into phi
    !$OMP PARALLEL DO
    do i=1,npfb
       phi(i) = sol_vec(i)
     !  print *, phi(i), ftn(rp(i,1),rp(i,2)), sqrt(rp(i,1)**2+rp(i,2)**2)

    end do
    !$OMP END PARALLEL DO
    deallocate(sol_vec)
!stop
 !   write(6,*) residual,iters

    return
  end subroutine solve_poisson_inhouse
!! ------------------------------------------------------------------------------------------------  

!! ------------------------------------------------------------------------------------------------  
!! ------------------------------------------------------------------------------------------------  
  subroutine solve_poisson_inhouse_dir(phi,rhs_ftn,tmp_vec)
    real(rkind),dimension(:),intent(out) :: phi
    real(rkind),dimension(:),intent(in) :: rhs_ftn,tmp_vec
    integer(ikind) :: i,j,k,ij,ii,num,n_t,nnz
    real(rkind) :: temp,tmp
    real(rkind),dimension(dims) :: riii
    integer(ikind) :: iters,itermax
    real(rkind) :: residual,tol
    real(rkind), dimension(:), allocatable :: rhs_vec,sol_vec
    integer(ikind), dimension(:),allocatable :: ij_num,ia
    real(rkind),dimension(:),allocatable :: p_fs 
    logical :: rebuild_lhs
  
    real(rkind) :: tt1,tt2,tt3,tt4
    
    !! Calculate the size of the PPE and allocate
    allocate(ia(npfb+1))
    ia(1)=npfb+2
    num=npfb+1
    do i=1,npfb
       ia(i+1) = ia(i) + ij_count(i)
    end do

    nnz = ia(npfb+1)
    allocate(lhs_mat(nnz),ija(nnz));ija=0;lhs_mat=0.0d0
    allocate(ij_num(npfb));ij_num=0
    allocate(rhs_vec(npfb+nb_n),sol_vec(npfb+nb_n));rhs_vec = 0.0d0;sol_vec=0.0d0

    !! Build the RHS
    !$OMP PARALLEL DO
    do i=1,npfb
       rhs_vec(i) = rhs_ftn(i)  
    end do
    !$OMP END PARALLEL DO

    if(.true.)then  !! Conditional in case I want to avoid rebuilding LHS

!       !$omp parallel do private(num,k,j,ij_num,temp,ij,n_t)
       do i=1,npfb   ! for all particles except compressible particles
          ija(i) = ia(i)
          !! First embedded loop builds ija so we can track where things go...
          num = ija(i) - 1
          do k=1,ij_count(i)   ! for all neighbours
             j=ij_link(i,k)     ! now considering particles i and j
             if(j.le.npfb)then
                num=num+1          
                ij_num(j)=num
                ija(num)=j
             endif
             if(j.gt.npfb)then
                num=num+1
                ij = irelation(j)
                ij_num(ij)=num
                ija(num) = ij
             end if
          enddo
          
          !! Second embedded loop is used to populate A
          do k=1,ij_count(i)
             j=ij_link(i,k)
             temp = -ij_w_lap(i,k)
             if(i.ne.j)then
                lhs_mat(i)=lhs_mat(i) + temp   ! Contribution to diagonal of A
                if(j.gt.npfb)then    ! If j is mirror, then
                   ij=irelation(j)    ! look for parent of j, ij
                   if(i.ne.ij)then
                      n_t = ij_num(ij)
                      lhs_mat(n_t)=lhs_mat(n_t)-temp   ! Contribution to A(i,ij)
                   else
                      lhs_mat(i)=lhs_mat(i)-temp       ! Contribution to A(i,ij) if ij=i 
                   end if
                elseif(j.le.npfb)then   ! if j is not a mirror
                   n_t = ij_num(j)   
                   lhs_mat(n_t)=lhs_mat(n_t)-temp    ! Contribution to A(i,j)
                end if
             end if
          end do

       end do
!       !$omp end parallel do
       ija(npfb+1)=ia(npfb+1)  !! Final element of ija            
      
       if(nb.ne.0) then !! Conditional for Dirichlets.
       !! This loop over neighbours begins to build Dirichlet condition
       !! By adding  -p_fs(i)*A_[i,:] to the RHS,   
       !$omp parallel do private(k,j) reduction(-:rhs_vec)
       do i=1,npfb  
          do k=ia(i),ia(i+1)-1   !! Loop over row
             j=ija(k)   !! We are considering element A_[i,j]
             if(j.le.nb) then
                rhs_vec(i)=rhs_vec(i) - lhs_mat(k)*tmp_vec(j)
             end if
          end do
       end do
       !$omp end parallel do
       
       !! Final loop will delete info in a columns and rows corresponding to a free surface node
       !! i.e. Set A_[:,j]=0,A_[j,j]=1 for all j in free surface
       !! also set RHS_[i]=p_fs[i] for all i in free surface
       !! I *think* this loop needs to be separate to previous loop.
       !$omp parallel do private(j,k)
       do i=1,npfb
          do k=ia(i),ia(i+1)-1
             j=ija(k)
             if(j.le.nb)then
                lhs_mat(k)=0.0d0
             end if
          end do
          if(i.le.nb) then
             lhs_mat(i) = 1.0_rkind   ! set A_[i,i] = 1
             rhs_vec(i) = tmp_vec(i)  ! set ith entry of RHS to p_fs(i)
             do k=ia(i),ia(i+1)-1   !! Loop over row
                j=ija(k)   !! We are considering element A_[i,j]
                lhs_mat(k)=0.0d0     !! and setting A_[i,:]=0,A_[i,i]=1
             end do
          end if     
       end do
       !$omp end parallel do

       end if
    end if
!stop
    !! Set tolerance and max iterations allowed
    tol = 1d-10;itermax=3000

    !! Solve the equation
    if(nb.ne.0)then !! If there are any Dirichlet boundary nodes
       call sparse_pbicgstab(lhs_mat,ija,sol_vec,rhs_vec,npfb,nnz,itermax,tol,iters,residual)
    else    !! If Neumann/periodic BCs only
       call sparse_pbicgstab_nonull(lhs_mat,ija,sol_vec,rhs_vec,npfb,nnz,itermax,tol,iters,residual)
    end if
    deallocate(ija,lhs_mat,rhs_vec)

    !! Pass solution into phi
    !$OMP PARALLEL DO
    do i=1,npfb
       phi(i) = sol_vec(i)
    end do
    !$OMP END PARALLEL DO
    deallocate(sol_vec)

  !  write(6,*) residual,iters

    return
  end subroutine solve_poisson_inhouse_dir
!! ------------------------------------------------------------------------------------------------  

subroutine poisson_solver_periodic(phi, rhs_ftn)
    real(rkind),dimension(:),intent(out) :: phi
    real(rkind),dimension(:),intent(in) :: rhs_ftn
    integer(ikind) :: i,j,k,n,info
    real(rkind), dimension(:), allocatable :: rhs_vec,ipiv
    real(rkind), dimension(:,:), allocatable :: l_mat
    real(rkind) :: tmp,x,y
    
    allocate(l_mat(np,np));l_mat=0.0d0
    allocate(rhs_vec(np),ipiv(np));rhs_vec=0.0d0;ipiv=0.0d0
    
    !Fill in LHS matrix  
    do i=1,np
    if (i<=npfb) then
      
       tmp=0.0d0
       do j=2,ij_count(i)
          k=ij_link(i, j)
          
            l_mat(i, k)=ij_w_lap(i,j)
            tmp=tmp-ij_w_lap(i,j)
            ! print *, k
             
          
         
       end do
       l_mat(i,i)=tmp
            
         !   print *, l_mat(i,i), sum(ij_w_lap(i,:))

    
       rhs_vec(i)=rhs_ftn(i)

   
    !stop
    elseif (i>npfb) then
      x=rp(i,1);y=rp(i,2);
      l_mat(i,i)=1.0d0
      rhs_vec(i)=ftn(x,y)
    
    
    end if
    end do
    
    call dgesv(np,1,l_mat,np,ipiv,rhs_vec,np,info)   
  
    phi=rhs_vec(1:npfb)
    do i=1,npfb
      x=rp(i,1)
      y=rp(i,2)
    !  print *, rhs_vec(i),sin(2.0*pi*x)*sin(2.0*pi*y)
    end do
   ! stop
    
    deallocate(l_mat,rhs_vec,ipiv)
    
end subroutine poisson_solver_periodic
   

end module poisson_inhouse