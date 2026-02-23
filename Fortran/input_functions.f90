module input_functions
 use kind_parameters,only: dp
 implicit none
 private
 public :: f_func,f_derivative,f_integral

contains
 
 real(dp) function f_func(t)
  !!user-defined function for f
  implicit none
  
  !!inputs
  real(dp),intent(in) :: t
  
  !!internal variables
  real(dp), parameter :: pii=3.1415926535897932_dp
  
  real(dp), parameter :: a=100._dp !!first sample problem
  real(dp), parameter :: b=100._dp
  
  !real(dp), parameter :: a=52.969966897973720_dp !a=600._dp/(pii*sqrt(13._dp)) !!second sample problem
  !real(dp), parameter :: b=100._dp
  !real(dp), parameter :: c=50._dp
  !real(dp), parameter :: d=397.27475173480286_dp !d=4500._dp/(pii*sqrt(13._dp))
  
  f_func=a*sin(pii*b*t) !!first sample problem
  !f_func=a*sin(pii*b*t)*exp(-c*t)+d !!second sample problem
 end function f_func
 
 real(dp) function f_derivative(t)
  !!user-defined function for derivative of f
  implicit none
  
  !!inputs
  real(dp),intent(in) :: t
  
  !!internal variables
  real(dp), parameter :: pii=3.1415926535897932_dp
  
  real(dp), parameter :: a=100._dp !!first sample problem
  real(dp), parameter :: b=100._dp
  
  !real(dp), parameter :: a=52.969966897973720_dp !a=600._dp/(pii*sqrt(13._dp)) !!second sample problem
  !real(dp), parameter :: b=100._dp
  !real(dp), parameter :: c=50._dp
  !real(dp), parameter :: d=397.27475173480286_dp !d=4500._dp/(pii*sqrt(13._dp))
  
  f_derivative=pii*a*b*cos(pii*b*t) !!first sample problem
  !f_derivative=(pii*b*cos(pii*b*t) - c*sin(pii*b*t))*a*exp(-c*t) !!second sample problem
 end function f_derivative
 
 real(dp) function f_integral(t)
  !!user-defined function for integral of f
  implicit none
  
  !!inputs
  real(dp),intent(in) :: t
  
  !!internal variables
  real(dp), parameter :: pii=3.1415926535897932_dp
  
  real(dp), parameter :: a=100._dp !!first sample problem
  real(dp), parameter :: b=100._dp
  
  !real(dp), parameter :: a=52.969966897973720_dp !a=600._dp/(pii*sqrt(13._dp)) !!second sample problem
  !real(dp), parameter :: b=100._dp
  !real(dp), parameter :: c=50._dp
  !real(dp), parameter :: d=397.27475173480286_dp !d=4500._dp/(pii*sqrt(13._dp))
  
  f_integral=(a/(pii*b))*(1._dp-cos(pii*b*t)) !!first example from "Sample results" section using f(t)=a*sin(pi*b*t)
  !f_integral=pii*a*b/((pii*b)**2+c**2)+d*t-(pii*b*cos(pii*b*t) + c*sin(pii*b*t))*a*exp(-c*t)/((pii*b)**2+c**2) !!!second sample problem
 end function f_integral

end module input_functions
