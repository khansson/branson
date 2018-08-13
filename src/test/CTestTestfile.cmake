# CMake generated Testfile for 
# Source directory: /Users/khansson/Code/branson/src/test
# Build directory: /Users/khansson/Code/branson/src/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_counter_rng "-n" "1" "./test_counter_rng_exe")
add_test(test_cell "-n" "1" "./test_cell_exe")
add_test(test_photon "-n" "1" "./test_photon_exe")
add_test(test_buffer "-n" "1" "./test_buffer_exe")
add_test(test_work_packet "-n" "1" "./test_work_packet_exe")
add_test(test_mpi_types "-n" "1" "./test_mpi_types_exe")
add_test(test_completion_manager "-n" "7" "./test_completion_manager_exe")
add_test(test_remap_census "-n" "8" "./test_remap_census_exe")
add_test(test_load_balance "-n" "4" "./test_load_balance_exe")
add_test(test_tally_manager "-n" "8" "./test_tally_manager_exe")
add_test(test_parmetis_1 "-n" "1" "./test_parmetis_1_exe")
add_test(test_parmetis_4 "-n" "4" "./test_parmetis_4_exe")
add_test(test_input "-n" "1" "./test_input_exe")
add_test(test_mesh "-n" "1" "./test_mesh_exe")
add_test(test_imc_state "-n" "2" "./test_imc_state_exe")
