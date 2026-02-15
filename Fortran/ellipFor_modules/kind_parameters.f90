module kind_parameters
!!define portable kind parameters 
 use, intrinsic :: iso_fortran_env, only : int32,int64,real32,real64,real128
implicit none
integer, parameter :: sp = real32
integer, parameter :: dp = real64
integer, parameter :: qp = real128
integer, parameter :: isp = int32
integer, parameter :: idp = int64
end module kind_parameters
