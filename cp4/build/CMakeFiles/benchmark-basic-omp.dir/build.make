# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

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

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.30.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.30.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/satvikverma/Workspace/CSC-746/cp4

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/satvikverma/Workspace/CSC-746/cp4/build

# Include any dependencies generated for this target.
include CMakeFiles/benchmark-basic-omp.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/benchmark-basic-omp.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/benchmark-basic-omp.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/benchmark-basic-omp.dir/flags.make

CMakeFiles/benchmark-basic-omp.dir/dgemm-basic-omp.cpp.o: CMakeFiles/benchmark-basic-omp.dir/flags.make
CMakeFiles/benchmark-basic-omp.dir/dgemm-basic-omp.cpp.o: /Users/satvikverma/Workspace/CSC-746/cp4/dgemm-basic-omp.cpp
CMakeFiles/benchmark-basic-omp.dir/dgemm-basic-omp.cpp.o: CMakeFiles/benchmark-basic-omp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/satvikverma/Workspace/CSC-746/cp4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/benchmark-basic-omp.dir/dgemm-basic-omp.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/benchmark-basic-omp.dir/dgemm-basic-omp.cpp.o -MF CMakeFiles/benchmark-basic-omp.dir/dgemm-basic-omp.cpp.o.d -o CMakeFiles/benchmark-basic-omp.dir/dgemm-basic-omp.cpp.o -c /Users/satvikverma/Workspace/CSC-746/cp4/dgemm-basic-omp.cpp

CMakeFiles/benchmark-basic-omp.dir/dgemm-basic-omp.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/benchmark-basic-omp.dir/dgemm-basic-omp.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/satvikverma/Workspace/CSC-746/cp4/dgemm-basic-omp.cpp > CMakeFiles/benchmark-basic-omp.dir/dgemm-basic-omp.cpp.i

CMakeFiles/benchmark-basic-omp.dir/dgemm-basic-omp.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/benchmark-basic-omp.dir/dgemm-basic-omp.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/satvikverma/Workspace/CSC-746/cp4/dgemm-basic-omp.cpp -o CMakeFiles/benchmark-basic-omp.dir/dgemm-basic-omp.cpp.s

# Object files for target benchmark-basic-omp
benchmark__basic__omp_OBJECTS = \
"CMakeFiles/benchmark-basic-omp.dir/dgemm-basic-omp.cpp.o"

# External object files for target benchmark-basic-omp
benchmark__basic__omp_EXTERNAL_OBJECTS = \
"/Users/satvikverma/Workspace/CSC-746/cp4/build/CMakeFiles/benchmark.dir/benchmark.cpp.o"

benchmark-basic-omp: CMakeFiles/benchmark-basic-omp.dir/dgemm-basic-omp.cpp.o
benchmark-basic-omp: CMakeFiles/benchmark.dir/benchmark.cpp.o
benchmark-basic-omp: CMakeFiles/benchmark-basic-omp.dir/build.make
benchmark-basic-omp: CMakeFiles/benchmark-basic-omp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/satvikverma/Workspace/CSC-746/cp4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable benchmark-basic-omp"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/benchmark-basic-omp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/benchmark-basic-omp.dir/build: benchmark-basic-omp
.PHONY : CMakeFiles/benchmark-basic-omp.dir/build

CMakeFiles/benchmark-basic-omp.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/benchmark-basic-omp.dir/cmake_clean.cmake
.PHONY : CMakeFiles/benchmark-basic-omp.dir/clean

CMakeFiles/benchmark-basic-omp.dir/depend:
	cd /Users/satvikverma/Workspace/CSC-746/cp4/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/satvikverma/Workspace/CSC-746/cp4 /Users/satvikverma/Workspace/CSC-746/cp4 /Users/satvikverma/Workspace/CSC-746/cp4/build /Users/satvikverma/Workspace/CSC-746/cp4/build /Users/satvikverma/Workspace/CSC-746/cp4/build/CMakeFiles/benchmark-basic-omp.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/benchmark-basic-omp.dir/depend
