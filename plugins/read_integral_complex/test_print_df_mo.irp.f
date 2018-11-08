program print_stuff

  PROVIDE ezfio_filename
  call ezfio_set_integrals_bielec_disk_access_df_ao_integral_array('Read')
  call run
end


BEGIN_PROVIDER [complex*16, m_array, (mo_num_per_kpt, mo_num_per_kpt, df_num, num_kpt_pairs)]
  implicit none
  BEGIN_DOC
  !
  END_DOC
  
  call fill_m(m_array,mo_num_per_kpt,df_num,num_kpt_pairs)

END_PROVIDER

subroutine run

PROVIDE m_array
call print_mo_bielec_from_df
print*,'done'

end

subroutine print_mo_bielec_from_df
  use map_module
  implicit none
  BEGIN_DOC
  ! fill ao bielec integral map using 3-index df integrals
  ! todo: add OMP directives, do smarter things with pointers/reshaping/transposing
  END_DOC

  integer :: i,k,j,l
  integer :: ki,kk,kj,kl
  integer :: ii,ik,ij,il
  integer :: kikk2,kjkl2,jl2,ik2

  complex*16,allocatable :: ints_ik(:,:), ints_jl(:,:), ints_ikjl_tmp(:,:), ints_ikjl(:,:,:,:)
  complex*16,allocatable :: ints_tmp1(:,:,:)!, ints_tmp2(:,:)

  complex*16 :: integral
  integer(key_kind)              :: tmp_idx1,tmp_idx2
  double precision               :: tmp_re,tmp_im
  integer                        :: mo_num_kpt_2
  mo_num_kpt_2 = mo_num_per_kpt * mo_num_per_kpt

  allocate( ints_jl(mo_num_kpt_2,df_num), &
    ints_tmp1(mo_num_per_kpt,mo_num_per_kpt,df_num),&
    ints_ik(mo_num_kpt_2,df_num),&
    ints_ikjl_tmp(mo_num_kpt_2,mo_num_kpt_2),&
    ints_ikjl(mo_num_per_kpt,mo_num_per_kpt,mo_num_per_kpt,mo_num_per_kpt)&
  )
  do kl=1, num_kpts
    do kj=1, kl
      call idx2_tri_int(kj,kl,kjkl2)
      ints_jl = reshape(m_array(:,:,:,kjkl2),(/mo_num_kpt_2,df_num/))

      do kk=1,kl
        ki=kconserv(kl,kk,kj)
        if ((kl == kj) .and. (ki > kk)) cycle
        call idx2_tri_int(ki,kk,kikk2)
        if (kikk2 > kjkl2) cycle
!        if (ki >= kk) then
        if (ki > kk) then
!          ints_ik = reshape( &
!                   conjg(reshape(df_ao_integral_array(:,:,:,kikk2),(/ao_num_per_kpt,ao_num_per_kpt,df_num/),order=(/2,1,3/))),&
!                   (/ao_num_per_kpt**2,df_num/))
          ints_tmp1 = conjg(reshape(m_array(:,:,:,kikk2),(/mo_num_per_kpt,mo_num_per_kpt,df_num/),order=(/2,1,3/)))

          ints_ik = reshape(ints_tmp1, (/mo_num_kpt_2,df_num/))
        else
          ints_ik = reshape(m_array(:,:,:,kikk2),(/mo_num_kpt_2,df_num/))
        endif

        call zgemm('N','T', mo_num_kpt_2, mo_num_kpt_2, df_num, &
               (1.d0,0.d0), ints_ik, size(ints_ik,1), &
               ints_jl, size(ints_jl,1), &
               (0.d0,0.d0), ints_ikjl_tmp, size(ints_ikjl_tmp,1))

        ! this is bad
        ! use a pointer instead?
        ints_ikjl = reshape(ints_ikjl_tmp,(/mo_num_per_kpt,mo_num_per_kpt,mo_num_per_kpt,mo_num_per_kpt/))
        
        do il=1,mo_num_per_kpt
          l=il+(kl-1)*mo_num_per_kpt
          do ij=1,mo_num_per_kpt
            j=ij+(kj-1)*mo_num_per_kpt
            if (j>l) exit
            call idx2_tri_int(j,l,jl2)
            do ik=1,mo_num_per_kpt
              k=ik+(kk-1)*mo_num_per_kpt
              if (k>l) exit
              do ii=1,mo_num_per_kpt
                i=ii+(ki-1)*mo_num_per_kpt
                if ((j==l) .and. (i>k)) exit
                call idx2_tri_int(i,k,ik2)
                if (ik2 > jl2) exit
                integral = ints_ikjl(ii,ik,ij,il)
                write(*,'(4(I3,X),2(F16.7))')i,j,k,l,real(integral),imag(integral)

              enddo !ii
            enddo !ik
          enddo !ij
        enddo !il
      enddo !kk
    enddo !kj
  enddo !kl

  deallocate( &
    ints_tmp1,&
    ints_ik,&
    ints_ikjl_tmp,&
    ints_ikjl&
    )
  
  deallocate( ints_jl ) 

end subroutine print_mo_bielec_from_df

subroutine fill_m(df_mo,n_mo,n_df,n_k_pairs)
  use map_module
  implicit none
  BEGIN_DOC
  ! create 3-idx mo ints from 3-idx ao ints
  END_DOC
  integer,intent(in) :: n_mo,n_df,n_k_pairs
  complex*16,intent(out) :: df_mo(n_mo,n_mo,n_df,n_k_pairs)
  integer :: kl,kj,kjkl2,mu
  complex*16,allocatable :: coef_l(:,:), coef_j(:,:), ints_jl(:,:), ints_tmp(:,:)

  integer :: i,j
  integer :: a,b
  complex*16 :: val1,val2

  print*,'providing 3-index MO integrals from 3-index AO integrals'

  allocate( &
            coef_l(ao_num_per_kpt,mo_num_per_kpt),&
            coef_j(ao_num_per_kpt,mo_num_per_kpt),&
            ints_jl(ao_num_per_kpt,ao_num_per_kpt),&
            ints_tmp(mo_num_per_kpt,ao_num_per_kpt)&
          )

  do kl=1, num_kpts
    coef_l = mo_coef_kpts_trunc(:,:,kl)
    do kj=1, kl
      coef_j = mo_coef_kpts_trunc(:,:,kj)
      call idx2_tri_int(kj,kl,kjkl2)
      do mu=1, df_num
        ints_jl = df_ao_integral_array(:,:,mu,kjkl2)
        call zgemm('C','N',mo_num_per_kpt,ao_num_per_kpt,ao_num_per_kpt, &
              (1.d0,0.d0), coef_j, ao_num_per_kpt, &
              ints_jl, ao_num_per_kpt, &
              (0.d0,0.d0), ints_tmp, mo_num_per_kpt)

        call zgemm('N','N',mo_num_per_kpt,mo_num_per_kpt,ao_num_per_kpt, &
              (1.d0,0.d0), ints_tmp, mo_num_per_kpt, &
              coef_l, ao_num_per_kpt, &
              (0.d0,0.d0), df_mo(:,:,mu,kjkl2), mo_num_per_kpt)
!        call zgemm('C','N',mo_num_per_kpt,ao_num_per_kpt,ao_num_per_kpt, &
!              (1.d0,0.d0), mo_coef_kpts_trunc(:,:,kj), ao_num_per_kpt, &
!              df_ao_integral_array(:,:,mu,kjkl2), ao_num_per_kpt, &
!              (0.d0,0.d0), ints_tmp, mo_num_per_kpt)
!
!        call zgemm('N','N',mo_num_per_kpt,mo_num_per_kpt,ao_num_per_kpt, &
!              (1.d0,0.d0), ints_tmp, mo_num_per_kpt, &
!              mo_coef_kpts_trunc(:,:,kl), ao_num_per_kpt, &
!              (0.d0,0.d0), df_mo(:,:,mu,kjkl2), mo_num_per_kpt)
      enddo
    enddo
  enddo

!  do kl=1, num_kpts
!    do kj=1, kl
!      call idx2_tri_int(kj,kl,kjkl2)
!      do mu=1, df_num
!        do i=1,mo_num_per_kpt
!          do j=1,mo_num_per_kpt
!            val1=df_mo(i,j,mu,kjkl2)
!            write(*,'(I3,X,I3,X,I3,X,I3,X,2(F16.7))')i,j,mu,kjkl2,real(val1),imag(val1)
!          enddo
!        enddo
!      enddo
!    enddo
!  enddo
  

  deallocate( &
            coef_l, &
            coef_j, &
            ints_jl, &
            ints_tmp &
          )

end subroutine fill_m
