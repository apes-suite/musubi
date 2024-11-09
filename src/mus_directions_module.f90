!> This module contains parameters to address the stencil directions. These
!! parameters are valid for all stencils, however, the direction `q00` and
!! `q000` needs to be set on each use according to `QQ`.
module mus_directions_module

  implicit none

  private

  integer, parameter, public :: qN00 = 1  !< west    x-
  integer, parameter, public :: q0N0 = 2  !< south   y-
  integer, parameter, public :: q00N = 3  !< bottom  z-
  integer, parameter, public :: q100 = 4  !< east    x+
  integer, parameter, public :: q010 = 5  !< north   y+
  integer, parameter, public :: q001 = 6  !< top     z+
  integer, parameter, public :: q0NN = 7  !<         z-,y-
  integer, parameter, public :: q0N1 = 8  !<         z+,y-
  integer, parameter, public :: q01N = 9  !<         z-,y+
  integer, parameter, public :: q011 = 10 !<         z+,y+
  integer, parameter, public :: qN0N = 11 !<         x-,z-
  integer, parameter, public :: q10N = 12 !<         x+,z-
  integer, parameter, public :: qN01 = 13 !<         x-,z+
  integer, parameter, public :: q101 = 14 !<         x+,z+
  integer, parameter, public :: qNN0 = 15 !<         y-,x-
  integer, parameter, public :: qN10 = 16 !<         y+,x-
  integer, parameter, public :: q1N0 = 17 !<         y-,x+
  integer, parameter, public :: q110 = 18 !<         y+,x+
  integer, parameter, public :: qNNN = 19 !<         x-,y-,z-
  integer, parameter, public :: qNN1 = 20 !<         x-,y-,z+
  integer, parameter, public :: qN1N = 21 !<         x-,y+,z-
  integer, parameter, public :: qN11 = 22 !<         x-,y+,z+
  integer, parameter, public :: q1NN = 23 !<         x+,y-,z-
  integer, parameter, public :: q1N1 = 24 !<         x+,y-,z+
  integer, parameter, public :: q11N = 25 !<         x+,y+,z-
  integer, parameter, public :: q111 = 26 !<         x+,y+,z+

  integer, parameter, public :: qN0 = 1 !< west   x-
  integer, parameter, public :: q0N = 2 !< south  y-
  integer, parameter, public :: q10 = 3 !< east   x+
  integer, parameter, public :: q01 = 4 !< north  y+
  integer, parameter, public :: qNN = 5 !<        y-,x-
  integer, parameter, public :: qN1 = 6 !<        y+,x-
  integer, parameter, public :: q1N = 7 !<        y-,x+
  integer, parameter, public :: q11 = 8 !<        y+,x+

  integer, parameter, public :: q__W = 1 !< west   x-
  integer, parameter, public :: q__S = 2 !< south  y-
  integer, parameter, public :: q__E = 3 !< east   x+
  integer, parameter, public :: q__N = 4 !< north  y+
  integer, parameter, public :: q_SW = 5 !<        y-,x-
  integer, parameter, public :: q_NW = 6 !<        y+,x-
  integer, parameter, public :: q_SE = 7 !<        y-,x+
  integer, parameter, public :: q_NE = 8 !<        y+,x+

  ! QQ, q00 and q000 cannot be exported, as they are different for different
  ! stencils. They need to be defined in using modules.

end module
