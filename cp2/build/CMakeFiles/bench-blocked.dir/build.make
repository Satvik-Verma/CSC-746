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
CMAKE_SOURCE_DIR = /Users/satvikverma/Workspace/CSC-746/cp2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/satvikverma/Workspace/CSC-746/cp2/build

# Include any dependencies generated for this target.
include CMakeFiles/bench-blocked.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/bench-blocked.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/bench-blocked.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/bench-blocked.dir/flags.make

CMakeFiles/bench-blocked.dir/benchmark.cpp.o: CMakeFiles/bench-blocked.dir/flags.make
CMakeFiles/bench-blocked.dir/benchmark.cpp.o: /Users/satvikverma/Workspace/CSC-746/cp2/benchmark.cpp
CMakeFiles/bench-blocked.dir/benchmark.cpp.o: CMakeFiles/bench-blocked.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/satvikverma/Workspace/CSC-746/cp2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/bench-blocked.dir/benchmark.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/bench-blocked.dir/benchmark.cpp.o -MF CMakeFiles/bench-blocked.dir/benchmark.cpp.o.d -o CMakeFiles/bench-blocked.dir/benchmark.cpp.o -c /Users/satvikverma/Workspace/CSC-746/cp2/benchmark.cpp

CMakeFiles/bench-blocked.dir/benchmark.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/bench-blocked.dir/benchmark.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/satvikverma/Workspace/CSC-746/cp2/benchmark.cpp > CMakeFiles/bench-blocked.dir/benchmark.cpp.i

CMakeFiles/bench-blocked.dir/benchmark.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/bench-blocked.dir/benchmark.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/satvikverma/Workspace/CSC-746/cp2/benchmark.cpp -o CMakeFiles/bench-blocked.dir/benchmark.cpp.s

bench-blocked: CMakeFiles/bench-blocked.dir/benchmark.cpp.o
bench-blocked: CMakeFiles/bench-blocked.dir/build.make
.PHONY : bench-blocked

# Rule to build all files generated by this target.
CMakeFiles/bench-blocked.dir/build: bench-blocked
.PHONY : CMakeFiles/bench-blocked.dir/build

CMakeFiles/bench-blocked.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/bench-blocked.dir/cmake_clean.cmake
.PHONY : CMakeFiles/bench-blocked.dir/clean

CMakeFiles/bench-blocked.dir/depend:
	cd /Users/satvikverma/Workspace/CSC-746/cp2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/satvikverma/Workspace/CSC-746/cp2 /Users/satvikverma/Workspace/CSC-746/cp2 /Users/satvikverma/Workspace/CSC-746/cp2/build /Users/satvikverma/Workspace/CSC-746/cp2/build /Users/satvikverma/Workspace/CSC-746/cp2/build/CMakeFiles/bench-blocked.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/bench-blocked.dir/depend

