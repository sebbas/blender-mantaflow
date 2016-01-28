# Blender Mantaflow integration

<del>This is the current state of the Blender Mantaflow integration. It builds on previous integration steps which can be found in the Blender branch `soc-2014-fluid`.</del>

## DEPRECATED
No more development will be taking place in the branch `soc-2014-fluid`. The Blender Mantaflow integration can now be found in the `fluid-mantaflow` branch which is part of the official Blender repository. You should build using that code instead.

However, an up-to-date version of `fluid-mantaflow` will also be available in this Github repository for you to fork. Building works just as described below.

## Building (for the impatient)

### Mac OSX

 1. Get the official Blender code:
        
        $ mkdir blender-build
        $ cd blender-build
        $ git clone http://git.blender.org/blender.git
        $ cd blender
        $ git submodule update --init --recursive
        $ git submodule foreach git checkout master
        $ git submodule foreach git pull --rebase origin master
        $ cd ..
        $ mkdir lib
        $ cd lib
        $ svn checkout https://svn.blender.org/svnroot/bf-blender/trunk/lib/darwin-9.x.universal

 2. Check out the soc-2014-fluid repo and play the new Mantaflow integration (this github repo) on top:
        
        $ cd ../blender
        $ git checkout soc-2014-fluid
        $ git remote add github https://github.com/sebbas/BlenderMantaflow.git
        $ git pull github soc-2014-fluid

 3. Build the code!
        
        $ cd ..
        $ mkdir build_darwin
        $ cd build_darwin
        $ cmake ../blender/ -DCMAKE_C_COMPILER=your_path/blender-build/lib/darwin-9.x.universal/clang-omp-3.5/bin/clang -DCMAKE_CXX_COMPILER=your_path/blender-build/lib/darwin-9.x.universal/clang-omp-3.5/bin/clang++
        $ make
        $ make install

Run Blender from *your_path/build_darwin/bin/Blender.app*!

### Linux

 1. Get the official Blender code:

        $ mkdir blender-git
        $ cd blender-git
        $ git clone http://git.blender.org/blender.git
        $ cd blender
        $ git submodule update --init --recursive
        $ git submodule foreach git checkout master
        $ git submodule foreach git pull --rebase origin master

 2. Check out the soc-2014-fluid repo and play the new Mantaflow integration (this github repo) on top:
        
        $ git checkout soc-2014-fluid
        $ git remote add github https://github.com/sebbas/BlenderMantaflow.git
        $ git pull github soc-2014-fluid

 3. Build the code!
        
        $ make

## Building (with explanations)

Since this repository relies on a working Blender building environment, the easiest way to set it up is to play it on top of the official Blender code. You can do so by following these steps:

 1. If you haven't already, download the latest Blender [code](http://wiki.blender.org/index.php/Dev:Doc/Building_Blender). If you have, then make sure that it is *up-to-date*! Make sure to *only get* the code, we will build it later.  

 2. Navigate to the Blender source code directory and check out the *soc-2014-fluid* repository which contains the old Mantaflow integration:  

         $ cd build-blender/blender
         $ git checkout soc-2014-fluid


 3. Play this github repository, which contains the updated Mantaflow integration, on top. So still in the Blender source directory, you can do this by simply adding a new remote and pulling from there:  
        
        $ git remote add github https://github.com/sebbas/BlenderMantaflow.git
        $ git pull github soc-2014-fluid

 3. Now build the code, preferably with CMake. You only need to make sure that your compiler supports OpenMP as the internal Mantaflow version relies on that.  

 3.1  **Linux**: Under Linux you should normally not not have any problems with the compiler as the default GNU compiler is readily available and supports OpenMP.

 3.2  **Mac OSX**: For Mac OSX you need to specify a compiler different from the default Clang compiler since that one unfortunately does not support OpenMP. Luckily, the *lib* directory from the SVN which you downloaded along with the source code, contains an OpenMP enabled compiler! Use it by telling CMake where to find it:

        $ cmake ../blender/ -DCMAKE_C_COMPILER=your_path/blender-build/lib/darwin-9.x.universal/clang-omp-3.5/bin/clang -DCMAKE_CXX_COMPILER=your_path/blender-build/lib/darwin-9.x.universal/clang-omp-3.5/bin/clang++ 

 4. Finalise your build as you would when building the official Blender source code.


## Examples

Some fire renderings can be found on [Vimeo](https://vimeo.com/sebbas/videos). For the current development state, take a look at the latest videos!

## Known issues
- Every time Blender starts a new animation cycle, Mantaflow allocates new grids as the Blender constructor `FLUID_3D` (and `WTURBULENCE` if using the high resolution option) is called. This results in higher memory consumption than in the Vanilla Blender as allocated grids are not immediately garbage collected by Python.

So, for now, don't let the animation loop to often ;)

- The Mantaflow script that you obtain when exporting in Blender is *not* compatible with the current, official Mantaflow [distribution](http://mantaflow.com/download.html). We still need to merge it with the Blender development branch in Mantaflow.

- Wavelet turbulence flickers at intersection "cached <-> uncached" when animation was not cached entirely (Probably because of uncached uvGrids ...) 

## Troubleshooting

### Mac OSX

If you have problems building under Mac OS 10.11 "El Capitan" and with Xcode 7 (I did) then below are two possible workarounds

The first is to edit CMakeLists.txt, and change the OSX_SYSTEM setting manually, by changing:

        elseif(${MAC_SYS} MATCHES 14)
    -       set(OSX_SYSTEM 10.10)
    +       set(OSX_SYSTEM 10.11)

The second workaround is to use an older Mac OSX SDK, e.g. the SDK from 10.10. 

 1. Get a copy of the Mac OSX 10.10 SDK. It can be found in Xcode 6 under `Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.10.sdk`. You might find Xcode 6 itself on your Time Machine Backup or you can download it from [Apple](https://developer.apple.com/downloads/).
 
 2. Copy `MacOSX10.10.sdk` into your current Xcode 7 app, alongside the 10.11 SDK which is located in `/ApplicationsXcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/`.
 
 3. Edit CMakeLists.txt like this

            if(${MAC_SYS} MATCHES 15)
        -       set(OSX_SYSTEM 10.11)
        +       set(OSX_SYSTEM 10.10)


Hope this works for you, if you know a better fix then please let me know!
