# Electro

This is part research project and part game. We're using the Godot game engine
for the purposes of 3D graphics, moving the camera, interacting with the fields,
and potentially making a game out of this one day. We're using C++ because speed
matters, Julia doesn't have Godot bindings, and Rust (a) doesn't have a library
with nd convolutions and (b) while cleaner, is slower to write.

Anyone wishing to do computations with discreet exterior calculus might find
this code to be a good starting point. DGtal also has a DEC module, but you may
find this smaller code base easier to work with (and the plots low key more aesthetic).

## running

There's a Linux executable in `bin`, and there will be Mac and Windows
executables as soon as I get that working. In the meantime, you can install
Godot and run the project from the editor. Dependencies:

- Godot
- scons
- Optional: [run](https://github.com/akriegman/run)
- Windows: mingw-g++

I haven't tried running on Mac, there may be more dependencies. The code for
Windows is built, so you only need Godot if you just want to run it.

## controls

`w a s d shift space` to move the camera, and use the mouse to look around.
