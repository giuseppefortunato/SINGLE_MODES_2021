# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

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
CMAKE_COMMAND = /home/giuseppe/Programmi/clion/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/giuseppe/Programmi/clion/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /media/giuseppe/DATA/codes/git-RBM-steady

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/giuseppe/DATA/codes/git-RBM-steady/cmake-build-debug

# Include any dependencies generated for this target.
include src/CMakeFiles/SingleMODES_leave_one_out.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/SingleMODES_leave_one_out.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/SingleMODES_leave_one_out.dir/flags.make

src/CMakeFiles/SingleMODES_leave_one_out.dir/SingleMODES_leave_one_out.cpp.o: src/CMakeFiles/SingleMODES_leave_one_out.dir/flags.make
src/CMakeFiles/SingleMODES_leave_one_out.dir/SingleMODES_leave_one_out.cpp.o: ../src/SingleMODES_leave_one_out.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/giuseppe/DATA/codes/git-RBM-steady/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/SingleMODES_leave_one_out.dir/SingleMODES_leave_one_out.cpp.o"
	cd /media/giuseppe/DATA/codes/git-RBM-steady/cmake-build-debug/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SingleMODES_leave_one_out.dir/SingleMODES_leave_one_out.cpp.o -c /media/giuseppe/DATA/codes/git-RBM-steady/src/SingleMODES_leave_one_out.cpp

src/CMakeFiles/SingleMODES_leave_one_out.dir/SingleMODES_leave_one_out.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SingleMODES_leave_one_out.dir/SingleMODES_leave_one_out.cpp.i"
	cd /media/giuseppe/DATA/codes/git-RBM-steady/cmake-build-debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/giuseppe/DATA/codes/git-RBM-steady/src/SingleMODES_leave_one_out.cpp > CMakeFiles/SingleMODES_leave_one_out.dir/SingleMODES_leave_one_out.cpp.i

src/CMakeFiles/SingleMODES_leave_one_out.dir/SingleMODES_leave_one_out.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SingleMODES_leave_one_out.dir/SingleMODES_leave_one_out.cpp.s"
	cd /media/giuseppe/DATA/codes/git-RBM-steady/cmake-build-debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/giuseppe/DATA/codes/git-RBM-steady/src/SingleMODES_leave_one_out.cpp -o CMakeFiles/SingleMODES_leave_one_out.dir/SingleMODES_leave_one_out.cpp.s

# Object files for target SingleMODES_leave_one_out
SingleMODES_leave_one_out_OBJECTS = \
"CMakeFiles/SingleMODES_leave_one_out.dir/SingleMODES_leave_one_out.cpp.o"

# External object files for target SingleMODES_leave_one_out
SingleMODES_leave_one_out_EXTERNAL_OBJECTS =

src/SingleMODES_leave_one_out: src/CMakeFiles/SingleMODES_leave_one_out.dir/SingleMODES_leave_one_out.cpp.o
src/SingleMODES_leave_one_out: src/CMakeFiles/SingleMODES_leave_one_out.dir/build.make
src/SingleMODES_leave_one_out: libMODES.so.0.0.1
src/SingleMODES_leave_one_out: src/CMakeFiles/SingleMODES_leave_one_out.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/giuseppe/DATA/codes/git-RBM-steady/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable SingleMODES_leave_one_out"
	cd /media/giuseppe/DATA/codes/git-RBM-steady/cmake-build-debug/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SingleMODES_leave_one_out.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/SingleMODES_leave_one_out.dir/build: src/SingleMODES_leave_one_out

.PHONY : src/CMakeFiles/SingleMODES_leave_one_out.dir/build

src/CMakeFiles/SingleMODES_leave_one_out.dir/clean:
	cd /media/giuseppe/DATA/codes/git-RBM-steady/cmake-build-debug/src && $(CMAKE_COMMAND) -P CMakeFiles/SingleMODES_leave_one_out.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/SingleMODES_leave_one_out.dir/clean

src/CMakeFiles/SingleMODES_leave_one_out.dir/depend:
	cd /media/giuseppe/DATA/codes/git-RBM-steady/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/giuseppe/DATA/codes/git-RBM-steady /media/giuseppe/DATA/codes/git-RBM-steady/src /media/giuseppe/DATA/codes/git-RBM-steady/cmake-build-debug /media/giuseppe/DATA/codes/git-RBM-steady/cmake-build-debug/src /media/giuseppe/DATA/codes/git-RBM-steady/cmake-build-debug/src/CMakeFiles/SingleMODES_leave_one_out.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/SingleMODES_leave_one_out.dir/depend
