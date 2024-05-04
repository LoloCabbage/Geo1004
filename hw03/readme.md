# GEO1004 homework 3: BIM to Geo conversion using voxels


## Prerequisites

Ensure the following package is installed and can be found by CMake.

- CGAL
- Eigen3 (at least 3.10 version)

Also, ensure your C++ version is not lower than 11.

## File structure

- `cpp`: stores all the code related files.

    - `include`
  
      - `json.hpp`: needed for CityJSON process.
  
    - `src`
  
        - `main.cpp`: contains all the codes of this assignment.

  - `CMakeLists.txt`: The CMake file for compilation.
  
- `data`: stores all the input IFC and OBJ files as well as the IfcConvert software.

    - `Wellness_center_Sama.ifc`: this IFC file is modifed using BlenderBIM plugin, 
        removing all the furniture and most of the street-level elements such as roads 
        and parking lots. 
    - `Wellness_center_Sama.obj`: converted OBJ file from the IFC file, using IfcConvert software.
    - `IfcConvert.exe`: the software used for converting IFC file to OBJ file.

- `output`: contains the ouput CityJSON files.

    - `output_01.json`: the output CityJSON file with voxel size of 0.1m.
    - `output_02.json`: the output CityJSON file with voxel size of 0.2m.
    - `output_03.json`: the output CityJSON file with voxel size of 0.3m.

## Convert IFC to OBJ file
To convert the IFC file to OBJ file, locate to the data directory and run the `IfcConvert.exe` in
your terminal:

        $ cd data
        $ ./IfcConvert Wellness_center_Sama.ifc Wellness_center_Sama.obj
The conversion will generate both OBJ and MTL files and output to the `data` directory.

## Running the code

### Input files

- To compile and run (without specify the input files):

      $ cd cpp
      $ mkdir build
      $ cd build
      $ cmake ..
      $ make
      $ ./geo1004_hw03

    Running the code in this way will use the default input settings:
    
  - `../../data/Wellness_center_Sama.obj`
  - `../../data/Wellness_center_Sama.mtl`

- To run the code with your own specified input files:

      $ cd cpp
      $ mkdir build
      $ cd build
      $ cmake ..
      $ make
      $ ./geo1004_hw03 input.obj input.mtl
- Note here, if run with specified input files, you must ensure your input 
  contains both OBJ and MTL files, and the OBJ file must be the first input.

### Output file
The output file is called `output.json`, and it will be output to the `build`
directory (**NOT the `output` directory**).

### Change parameters

- `voxel size`: if you want to change the voxel size, go to the `main.cpp`,
    find line **800**:

        VoxelGrid voxel_grid = create_voxel_grid(obj_info, 0.2);
    the `0.2` means the voxel size is 0.2 meter, change it to the size you want.
- `room validation threshold`: we use the *length*, *width* and *height* of the 
    oriented bounding box of a room to check whether it is a valid room. You can
    find it in `main.cpp`, line **493**:

      if (bbox_x < 1 || bbox_y < 1 || bbox_z < 1) return false;
  
  by default, we define a valid room should be greater than 1 meter in *length*, 
  *width* and *height*, you can tune this threshold if you want.