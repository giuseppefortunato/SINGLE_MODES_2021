# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /media/giuseppe/DATA/codes/git-RBM-steady

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/giuseppe/DATA/codes/git-RBM-steady/build

# Include any dependencies generated for this target.
include src/CMakeFiles/SingleMODES.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/SingleMODES.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/SingleMODES.dir/flags.make

src/CMakeFiles/SingleMODES.dir/SingleMODES.cpp.o: src/CMakeFiles/SingleMODES.dir/flags.make
src/CMakeFiles/SingleMODES.dir/SingleMODES.cpp.o: ../src/SingleMODES.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/giuseppe/DATA/codes/git-RBM-steady/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/SingleMODES.dir/SingleMODES.cpp.o"
	cd /media/giuseppe/DATA/codes/git-RBM-steady/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SingleMODES.dir/SingleMODES.cpp.o -c /media/giuseppe/DATA/codes/git-RBM-steady/src/SingleMODES.cpp

src/CMakeFiles/SingleMODES.dir/SingleMODES.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SingleMODES.dir/SingleMODES.cpp.i"
	cd /media/giuseppe/DATA/codes/git-RBM-steady/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/giuseppe/DATA/codes/git-RBM-steady/src/SingleMODES.cpp > CMakeFiles/SingleMODES.dir/SingleMODES.cpp.i

src/CMakeFiles/SingleMODES.dir/SingleMODES.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SingleMODES.dir/SingleMODES.cpp.s"
	cd /media/giuseppe/DATA/codes/git-RBM-steady/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/giuseppe/DATA/codes/git-RBM-steady/src/SingleMODES.cpp -o CMakeFiles/SingleMODES.dir/SingleMODES.cpp.s

# Object files for target SingleMODES
SingleMODES_OBJECTS = \
"CMakeFiles/SingleMODES.dir/SingleMODES.cpp.o"

# External object files for target SingleMODES
SingleMODES_EXTERNAL_OBJECTS =

src/SingleMODES: src/CMakeFiles/SingleMODES.dir/SingleMODES.cpp.o
src/SingleMODES: src/CMakeFiles/SingleMODES.dir/build.make
src/SingleMODES: libMODES.so.0.0.1
src/SingleMODES: src/CMakeFiles/SingleMODES.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/giuseppe/DATA/codes/git-RBM-steady/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable SingleMODES"
	cd /media/giuseppe/DATA/codes/git-RBM-steady/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SingleMODES.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/SingleMODES.dir/build: src/SingleMODES

.PHONY : src/CMakeFiles/SingleMODES.dir/build

src/CMakeFiles/SingleMODES.dir/clean:
	cd /media/giuseppe/DATA/codes/git-RBM-steady/build/src && $(CMAKE_COMMAND) -P CMakeFiles/SingleMODES.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/SingleMODES.dir/clean

src/CMakeFiles/SingleMODES.dir/depend:
	cd /media/giuseppe/DATA/codes/git-RBM-steady/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/giuseppe/DATA/codes/git-RBM-steady /media/giuseppe/DATA/codes/git-RBM-steady/src /media/giuseppe/DATA/codes/git-RBM-steady/build /media/giuseppe/DATA/codes/git-RBM-steady/build/src /media/giuseppe/DATA/codes/git-RBM-steady/build/src/CMakeFiles/SingleMODES.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/SingleMODES.dir/depend

