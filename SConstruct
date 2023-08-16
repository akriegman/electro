#!/usr/bin/env python
import os
import sys

godot_cpp_path = "cpp/godot-cpp"
src_path = "cpp"
godot_project_path = "godot"

env = SConscript(godot_cpp_path + "/SConstruct")

# For reference:
# - CCFLAGS are compilation flags shared between C and C++
# - CFLAGS are for C-specific compilation flags
# - CXXFLAGS are for C++-specific compilation flags
# - CPPFLAGS are for pre-processor flags
# - CPPDEFINES are for pre-processor defines
# - LINKFLAGS are for linking flags

# tweak this if you want to use different folders, or more folders, to store your source code in.
env.Append(CPPPATH=[src_path])
sources = Glob(src_path + "/*.cpp")

if env["platform"] == "macos":
    library = env.SharedLibrary(
        godot_project_path + "/bin/libgdexample.{}.{}.framework/libgdexample.{}.{}".format(
            env["platform"], env["target"], env["platform"], env["target"]
        ),
        source=sources,
    )
else:
    library = env.SharedLibrary(
        godot_project_path + "/bin/libgdexample{}{}".format(env["suffix"], env["SHLIBSUFFIX"]),
        source=sources,
    )

Default(library)
