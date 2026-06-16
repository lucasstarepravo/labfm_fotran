module gnn_mod
  use iso_fortran_env, only: int64, real64, real32
  use kind_parameters
  use common_parameter
  use common_2d
  use torchfort
  implicit none
  private

  public :: import_gnn

  contains
      ! Create the instance of the pre-trained GNN
    subroutine import_gnn
      integer :: res

      res = torchfort_create_model('gnn_deriv', './gnn/deriv_config_cpu.yaml', -1)  
      print *, 'torchfort_create_model res =', res
      if (res /= TORCHFORT_RESULT_SUCCESS) stop 1

      res = torchfort_create_model('gnn_lap', './gnn/lap_config_cpu.yaml', -1)
      print *,  'torchfort_create_model res =', res
      if (res /= TORCHFORT_RESULT_SUCCESS) stop 1


    end subroutine import_gnn
end module gnn_mod
