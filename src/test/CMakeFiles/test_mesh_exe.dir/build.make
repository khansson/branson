# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.11.3/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.11.3/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/khansson/Code/branson/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/khansson/Code/branson/src

# Include any dependencies generated for this target.
include test/CMakeFiles/test_mesh_exe.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/test_mesh_exe.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/test_mesh_exe.dir/flags.make

test/CMakeFiles/test_mesh_exe.dir/test_mesh.cc.o: test/CMakeFiles/test_mesh_exe.dir/flags.make
test/CMakeFiles/test_mesh_exe.dir/test_mesh.cc.o: test/test_mesh.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/khansson/Code/branson/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/test_mesh_exe.dir/test_mesh.cc.o"
	cd /Users/khansson/Code/branson/src/test && mpic++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_mesh_exe.dir/test_mesh.cc.o -c /Users/khansson/Code/branson/src/test/test_mesh.cc

test/CMakeFiles/test_mesh_exe.dir/test_mesh.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_mesh_exe.dir/test_mesh.cc.i"
	cd /Users/khansson/Code/branson/src/test && mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/khansson/Code/branson/src/test/test_mesh.cc > CMakeFiles/test_mesh_exe.dir/test_mesh.cc.i

test/CMakeFiles/test_mesh_exe.dir/test_mesh.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_mesh_exe.dir/test_mesh.cc.s"
	cd /Users/khansson/Code/branson/src/test && mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/khansson/Code/branson/src/test/test_mesh.cc -o CMakeFiles/test_mesh_exe.dir/test_mesh.cc.s

# Object files for target test_mesh_exe
test_mesh_exe_OBJECTS = \
"CMakeFiles/test_mesh_exe.dir/test_mesh.cc.o"

# External object files for target test_mesh_exe
test_mesh_exe_EXTERNAL_OBJECTS =

test/test_mesh_exe: test/CMakeFiles/test_mesh_exe.dir/test_mesh.cc.o
test/test_mesh_exe: test/CMakeFiles/test_mesh_exe.dir/build.make
test/test_mesh_exe: test/CMakeFiles/test_mesh_exe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/khansson/Code/branson/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_mesh_exe"
	cd /Users/khansson/Code/branson/src/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_mesh_exe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/test_mesh_exe.dir/build: test/test_mesh_exe

.PHONY : test/CMakeFiles/test_mesh_exe.dir/build

test/CMakeFiles/test_mesh_exe.dir/clean:
	cd /Users/khansson/Code/branson/src/test && $(CMAKE_COMMAND) -P CMakeFiles/test_mesh_exe.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/test_mesh_exe.dir/clean

test/CMakeFiles/test_mesh_exe.dir/depend:
	cd /Users/khansson/Code/branson/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/khansson/Code/branson/src /Users/khansson/Code/branson/src/test /Users/khansson/Code/branson/src /Users/khansson/Code/branson/src/test /Users/khansson/Code/branson/src/test/CMakeFiles/test_mesh_exe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/test_mesh_exe.dir/depend

