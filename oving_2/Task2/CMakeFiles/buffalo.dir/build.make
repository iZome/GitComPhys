# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_SOURCE_DIR = /home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2

# Include any dependencies generated for this target.
include CMakeFiles/buffalo.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/buffalo.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/buffalo.dir/flags.make

CMakeFiles/buffalo.dir/BuffaloMain.cpp.o: CMakeFiles/buffalo.dir/flags.make
CMakeFiles/buffalo.dir/BuffaloMain.cpp.o: BuffaloMain.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/buffalo.dir/BuffaloMain.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/buffalo.dir/BuffaloMain.cpp.o -c /home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2/BuffaloMain.cpp

CMakeFiles/buffalo.dir/BuffaloMain.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/buffalo.dir/BuffaloMain.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2/BuffaloMain.cpp > CMakeFiles/buffalo.dir/BuffaloMain.cpp.i

CMakeFiles/buffalo.dir/BuffaloMain.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/buffalo.dir/BuffaloMain.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2/BuffaloMain.cpp -o CMakeFiles/buffalo.dir/BuffaloMain.cpp.s

CMakeFiles/buffalo.dir/BuffaloMain.cpp.o.requires:

.PHONY : CMakeFiles/buffalo.dir/BuffaloMain.cpp.o.requires

CMakeFiles/buffalo.dir/BuffaloMain.cpp.o.provides: CMakeFiles/buffalo.dir/BuffaloMain.cpp.o.requires
	$(MAKE) -f CMakeFiles/buffalo.dir/build.make CMakeFiles/buffalo.dir/BuffaloMain.cpp.o.provides.build
.PHONY : CMakeFiles/buffalo.dir/BuffaloMain.cpp.o.provides

CMakeFiles/buffalo.dir/BuffaloMain.cpp.o.provides.build: CMakeFiles/buffalo.dir/BuffaloMain.cpp.o


CMakeFiles/buffalo.dir/Buffalo.cpp.o: CMakeFiles/buffalo.dir/flags.make
CMakeFiles/buffalo.dir/Buffalo.cpp.o: Buffalo.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/buffalo.dir/Buffalo.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/buffalo.dir/Buffalo.cpp.o -c /home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2/Buffalo.cpp

CMakeFiles/buffalo.dir/Buffalo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/buffalo.dir/Buffalo.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2/Buffalo.cpp > CMakeFiles/buffalo.dir/Buffalo.cpp.i

CMakeFiles/buffalo.dir/Buffalo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/buffalo.dir/Buffalo.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2/Buffalo.cpp -o CMakeFiles/buffalo.dir/Buffalo.cpp.s

CMakeFiles/buffalo.dir/Buffalo.cpp.o.requires:

.PHONY : CMakeFiles/buffalo.dir/Buffalo.cpp.o.requires

CMakeFiles/buffalo.dir/Buffalo.cpp.o.provides: CMakeFiles/buffalo.dir/Buffalo.cpp.o.requires
	$(MAKE) -f CMakeFiles/buffalo.dir/build.make CMakeFiles/buffalo.dir/Buffalo.cpp.o.provides.build
.PHONY : CMakeFiles/buffalo.dir/Buffalo.cpp.o.provides

CMakeFiles/buffalo.dir/Buffalo.cpp.o.provides.build: CMakeFiles/buffalo.dir/Buffalo.cpp.o


# Object files for target buffalo
buffalo_OBJECTS = \
"CMakeFiles/buffalo.dir/BuffaloMain.cpp.o" \
"CMakeFiles/buffalo.dir/Buffalo.cpp.o"

# External object files for target buffalo
buffalo_EXTERNAL_OBJECTS =

buffalo: CMakeFiles/buffalo.dir/BuffaloMain.cpp.o
buffalo: CMakeFiles/buffalo.dir/Buffalo.cpp.o
buffalo: CMakeFiles/buffalo.dir/build.make
buffalo: /usr/lib/libarmadillo.so
buffalo: /usr/lib/liblapack.so
buffalo: /usr/lib/x86_64-linux-gnu/libgsl.so
buffalo: /usr/lib/x86_64-linux-gnu/libgslcblas.so
buffalo: CMakeFiles/buffalo.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable buffalo"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/buffalo.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/buffalo.dir/build: buffalo

.PHONY : CMakeFiles/buffalo.dir/build

CMakeFiles/buffalo.dir/requires: CMakeFiles/buffalo.dir/BuffaloMain.cpp.o.requires
CMakeFiles/buffalo.dir/requires: CMakeFiles/buffalo.dir/Buffalo.cpp.o.requires

.PHONY : CMakeFiles/buffalo.dir/requires

CMakeFiles/buffalo.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/buffalo.dir/cmake_clean.cmake
.PHONY : CMakeFiles/buffalo.dir/clean

CMakeFiles/buffalo.dir/depend:
	cd /home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2 /home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2 /home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2 /home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2 /home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2/CMakeFiles/buffalo.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/buffalo.dir/depend

