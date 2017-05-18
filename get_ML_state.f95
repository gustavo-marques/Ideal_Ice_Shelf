! Compute mean tracer properties in the mixed layer.
! Author: Gustavo Marques
! compile: f2py -c get_ML_state.f95 -m get_ML_state

subroutine ml_average(SST,SSS,temp,salt,hml,h,km,jm,im)
implicit none
real (kind=8), dimension(1:km,1:jm,1:im), intent(in):: temp
real (kind=8), dimension(1:km,1:jm,1:im), intent(in):: salt
real (kind=8), dimension(1:km,1:jm,1:im), intent(in):: h
real (kind=8), dimension(1:jm,1:im), intent(in):: hml
real (kind=8), dimension(1:jm,1:im), intent(out):: SST
real (kind=8), dimension(1:jm,1:im), intent(out):: SSS

integer, intent(in) :: km,jm,im
real (kind=8), dimension(1:im) :: depth
real :: dh
integer :: i,j,k

SSS(:,:) = 0.0; SST(:,:) = 0.0
do j=1,jm
      depth(:) = 0.0
      do k=1,km ; do i=1,im
        if (depth(i) + h(k,j,i) < hml(j,i)) then
          dh = h(k,j,i)
        elseif (depth(i) < hml(j,i)) then
          dh = hml(j,i) - depth(i)
        else
          dh = 0.0
        endif
        SST(j,i) = SST(j,i) + dh * temp(k,j,i)
        SSS(j,i) = SSS(j,i) + dh * salt(k,j,i)
        depth(i) = depth(i) + dh
      enddo ; enddo

! Calculate the average properties of the mixed layer depth.
      do i=1,im
          SST(j,i) = SST(j,i) / depth(i)
          SSS(j,i) = SSS(j,i) / depth(i)
          !write(*,*)'SST,SSS,depth',SST(j,i),SSS(j,i),depth(i)
      enddo
 enddo ! end of j loop


! end of program
end subroutine
