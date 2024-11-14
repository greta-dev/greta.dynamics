check_state_is_2d_or_3d <- function(state_n_dim){
  state_is_2d_oe_3e <- state_n_dim %in% 2:3
  if (!state_is_2d_oe_3e) {
    cli::cli_abort(
      "State must be either two- or three-dimensional"
    )
  }
}

check_initial_state_col_vec_or_3d_dim_1 <- function(state){
  state_dim <- dim(state)
  state_n_dim <- length(state_dim)
  last_dim_of_state_is_not_1 <- state_dim[state_n_dim] != 1
  if (last_dim_of_state_is_not_1) {
    cli::cli_abort(
      "{.var initial_state} must be either a column vector, or a 3D array \\
      with final dimension 1"
    )
  }

}
