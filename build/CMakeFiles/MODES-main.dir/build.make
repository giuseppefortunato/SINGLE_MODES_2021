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
include CMakeFiles/MODES-main.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/MODES-main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MODES-main.dir/flags.make

CMakeFiles/MODES-main.dir/src/main.cpp.o: CMakeFiles/MODES-main.dir/flags.make
CMakeFiles/MODES-main.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/giuseppe/DATA/codes/git-RBM-steady/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/MODES-main.dir/src/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MODES-main.dir/src/main.cpp.o -c /media/giuseppe/DATA/codes/git-RBM-steady/src/main.cpp

CMakeFiles/MODES-main.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MODES-main.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/giuseppe/DATA/codes/git-RBM-steady/src/main.cpp > CMakeFiles/MODES-main.dir/src/main.cpp.i

CMakeFiles/MODES-main.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MODES-main.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/giuseppe/DATA/codes/git-RBM-steady/src/main.cpp -o CMakeFiles/MODES-main.dir/src/main.cpp.s

# Object files for target MODES-main
MODES__main_OBJECTS = \
"CMakeFiles/MODES-main.dir/src/main.cpp.o"

# External object files for target MODES-main
MODES__main_EXTERNAL_OBJECTS =

MODES-main: CMakeFiles/MODES-main.dir/src/main.cpp.o
MODES-main: CMakeFiles/MODES-main.dir/build.make
MODES-main: libMODES.so.0.0.1
MODES-main: CMakeFiles/MODES-main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/giuseppe/DATA/codes/git-RBM-steady/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable MODES-main"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MODES-main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MODES-main.dir/build: MODES-main

.PHONY : CMakeFiles/MODES-main.dir/build

CMakeFiles/MODES-main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/MODES-main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/MODES-main.dir/clean

CMakeFiles/MODES-main.dir/depend:
	cd /media/giuseppe/DATA/codes/git-RBM-steady/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/giuseppe/DATA/codes/git-RBM-steady /media/giuseppe/DATA/codes/git-RBM-steady /media/giuseppe/DATA/codes/git-RBM-steady/build /media/giuseppe/DATA/codes/git-RBM-steady/build /media/giuseppe/DATA/codes/git-RBM-steady/build/CMakeFiles/MODES-main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/MODES-main.dir/depend
