# dynamics module errors informatively

    Code
      iterate_matrix(matrix = bad_mat, initial_state = good_state)
    Condition
      Error in `iterate_matrix()`:
      ! Each matrix must be a two-dimensional square greta array

---

    Code
      iterate_matrix(matrix = bad_matrices1, initial_state = good_state)
    Condition
      Error in `iterate_matrix()`:
      ! matrix and state must be either two- or three-dimensional

---

    Code
      iterate_matrix(matrix = bad_matrices1, initial_state = good_states)
    Condition
      Error in `iterate_matrix()`:
      ! matrix and state must be either two- or three-dimensional

---

    Code
      iterate_matrix(matrix = bad_matrices2, initial_state = good_state)
    Condition
      Error in `iterate_matrix()`:
      ! Each matrix must be a two-dimensional square greta array

---

    Code
      iterate_matrix(matrix = bad_matrices2, initial_state = good_states)
    Condition
      Error in `iterate_matrix()`:
      ! Each matrix must be a two-dimensional square greta array

---

    Code
      iterate_matrix(matrix = good_mat, initial_state = bad_state)
    Condition
      Error in `check_initial_state_col_vec_or_3d_dim_1()`:
      ! `initial_state` must be either a column vector, or a 3D array with final dimension 1

---

    Code
      iterate_matrix(matrix = good_matrices, initial_state = bad_state)
    Condition
      Error in `check_initial_state_col_vec_or_3d_dim_1()`:
      ! `initial_state` must be either a column vector, or a 3D array with final dimension 1

---

    Code
      iterate_matrix(matrix = good_mat, initial_state = mismatched_state)
    Condition
      Error in `iterate_matrix()`:
      ! Length of each initial_state must match the dimension of matrix

---

    Code
      iterate_matrix(matrix = good_matrices, initial_state = mismatched_state)
    Condition
      Error in `iterate_matrix()`:
      ! Length of each initial_state must match the dimension of matrix

